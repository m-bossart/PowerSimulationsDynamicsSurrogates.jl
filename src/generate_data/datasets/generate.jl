abstract type SurrogateDataset end

struct GenerateDataParams
    solver::String
    solver_tols::NamedTuple{(:reltol, :abstol), Tuple{Float64, Float64}}
    dtmax::Union{Nothing, Float64}
    tspan::Tuple{Float64, Float64}
    tstops::Vector{Float64}
    tsave::Vector{Float64}
    formulation::String
    all_branches_dynamic::Bool
    all_lines_dynamic::Bool
    frequency_reference::String
    seed::Int64
end

function GenerateDataParams(;
    solver = "IDA",
    solver_tols = (reltol = 1e-6, abstol = 1e-6),
    dtmax = nothing,
    tspan = (0.0, 1.0),
    tstops = [],
    tsave = [],
    formulation = "Residual",
    all_branches_dynamic = false,
    all_lines_dynamic = false,
    frequency_reference = "ReferenceBus",
    seed = 1,
)
    GenerateDataParams(
        solver,
        solver_tols,
        dtmax,
        tspan,
        tstops,
        tsave,
        formulation,
        all_branches_dynamic,
        all_lines_dynamic,
        frequency_reference,
        seed,
    )
end

"""
    function generate_surrogate_data(
        dataset_type::Type{T},
        sys_main::PSY.System,
        sys_aux::PSY.System, 
        ics::Vector{Vector{Float64}},
        data_params::D,
        data_collection_params::GenerateDataParams,
        dataset_aux = nothing,
    ) where {O <: SurrogateOperatingPoint, T <: SurrogateDataset}

This function creates a new deepcopy of `sys_main` for each initial condition and then generates a dataset acording to `data_params` and `data_collection_params`.
- `sys_main` is the system that is changed and perturbed and used to generate data. \\
- `ics` Vector of the initial conditions to be used for the simulation. Used inplace of combition of perturbations and operating points. \\
- `data_params` are the parameters for generating a specific type of dataset. \\
- `data_collection_params` are the generic parameters for generating time series data.
- `dataset_aux` - An auxiliary dataset. If provided, is passed to `fill_surrogate_data!`. One possible use of this fucntionality it to avoid re-running unstable cases.  
- `surrogate_params` - Parameters of a surrogate model. Passed in order to match up the operating point appropriately. 
"""
function generate_surrogate_data(
    dataset_type::Type{T},
    sys_main::PSY.System,
    sys_aux::PSY.System,
    ics::Vector{Vector{Float64}},
    operating_points::Vector{O},
    device_details::Dict{String, Dict{Symbol, Symbol}},
    data_collection_params::GenerateDataParams;
    dataset_aux = nothing,
    surrogate_params = nothing,
    run_precompilation_sim = false,
) where {O <: SurrogateOperatingPoint, T <: SurrogateDataset}
    Random.seed!(data_collection_params.seed)
    train_data = dataset_type[]
    ########################################################################################
    #Run full simulation for first operating point and perturbation in order to precompile
    #all code and avoid including in timing.
    ########################################################################################
    if run_precompilation_sim
        sys = deepcopy(sys_main)
        update_operating_point!(sys, operating_points[1], sys_aux)
        dummy_data = EmptyTrainDataSet(dataset_type)
        sim_full =
            _build_run_simulation_initial_conditions(sys, data_collection_params, ics[1])
        if !issubset(
            data_collection_params.tsave,
            union(data_collection_params.tstops, sim_full.tstops),
        )
            @warn "tsave not subset of tstops"
        end
        fill_surrogate_data!(dummy_data, device_details, data_collection_params, sim_full)
    end
    ########################################################################################
    ########################################################################################
    @assert length(operating_points) == length(ics)
    for (ix, ic) in enumerate(ics)
        sys = deepcopy(sys_main)
        update_operating_point!(sys, operating_points[ix], sys_aux)
        #PowerFlows.solve_powerflow!(PowerFlows.ACPowerFlow(),sys) 
        data = EmptyTrainDataSet(dataset_type)
        sim_full = _build_run_simulation_initial_conditions(sys, data_collection_params, ic)
        if !issubset(
            data_collection_params.tsave,
            union(data_collection_params.tstops, sim_full.tstops),
        )
            @warn "tsave not subset of tstops"
        end
        fill_surrogate_data!(data, device_details, data_collection_params, sim_full)
        @show PSY.get_internal_voltage(collect(PSY.get_components(PSY.Source, sys))[1])
        @show PSY.get_internal_angle(collect(PSY.get_components(PSY.Source, sys))[1])
        push!(train_data, data)
    end
    return train_data
end

"""
    function generate_surrogate_data(
        dataset_type::Type{T},
        sys_main::PSY.System,
        sys_aux::PSY.System,          
        perturbations::Vector{Vector{Union{SurrogatePerturbation, PSID.Perturbation}}},
        operating_points::Vector{O},
        data_params::D,
        data_collection_params::GenerateDataParams,
        dataset_aux = nothing,
    ) where {O <: SurrogateOperatingPoint, T <: SurrogateDataset}

This function creates a new deepcopy of `sys_main` for each pair of (`perturbations`, `operating_points`) and then generates a dataset acording to `data_params` and `data_collection_params`.
- `sys_main` is the system that is changed and perturbed and used to generate data. \\
- `sys_aux` is used by some operating points and perturbations which require randomly choosing components to perturb. \\
- `perturbations` is a vector of perturbations to be implemented combinatorially with `operating_points` \\
- `operating_points` is a vector of operating points to be implemented combinatorially with `perturbations` \\
- `data_params` are the parameters for generating a specific type of dataset. \\
- `data_collection_params` are the generic parameters for generating time series data.
- `dataset_aux` - An auxiliary dataset. If provided, is passed to `fill_surrogate_data!`. One possible use of this fucntionality it to avoid re-running unstable cases.  
- `surrogate_params` - Parameters of a surrogate model. Passed in order to match up the operating point appropriately. 
"""
function generate_surrogate_data(
    dataset_type::Type{T},
    sys_main::PSY.System,
    sys_aux::PSY.System,    #TODO make this optional
    perturbations::Vector{Vector{Union{SurrogatePerturbation, PSID.Perturbation}}},
    operating_points::Vector{O},
    device_details::Dict{String, Dict{Symbol, Symbol}},
    data_collection_params::GenerateDataParams;
    dataset_aux = nothing,
    surrogate_params = nothing,
    run_precompilation_sim = false,
) where {O <: SurrogateOperatingPoint, T <: SurrogateDataset}
    Random.seed!(data_collection_params.seed)
    train_data = dataset_type[]
    ########################################################################################
    #Run full simulation for first operating point and perturbation in order to precompile
    #all code and avoid including in timing.
    ########################################################################################
    if run_precompilation_sim
        sys = deepcopy(sys_main)
        update_operating_point!(sys, operating_points[1], sys_aux)
        PSID.Simulation!(PSID.MassMatrixModel, sys, pwd(), (0.0, 0.0))  #Run power flow and re-init devices before defining perturbation
        psid_perturbations = PSID.Perturbation[]
        for p_single in perturbations[1]
            add_surrogate_perturbation!(sys, psid_perturbations, p_single, sys_aux)
        end
        rng_state = copy(Random.default_rng())
        dummy_data = EmptyTrainDataSet(dataset_type)
        if dataset_aux !== nothing && dataset_aux[1].built == true
            match_operating_point(sys, dataset_aux[1], surrogate_params)    #TODO - not tested
            sim_full = _build_run_simulation_perturbations(
                sys,
                data_collection_params,
                psid_perturbations,
            )
            if !issubset(
                data_collection_params.tsave,
                union(data_collection_params.tstops, sim_full.tstops),
            )
                @warn "tsave not subset of tstops"
            end
            fill_surrogate_data!(
                dummy_data,
                device_details,
                data_collection_params,
                sim_full,
            )
        elseif dataset_aux === nothing
            sim_full = _build_run_simulation_perturbations(
                sys,
                data_collection_params,
                psid_perturbations,
            )
            if !issubset(
                data_collection_params.tsave,
                union(data_collection_params.tstops, sim_full.tstops),
            )
                @warn "tsave not subset of tstops"
            end
            fill_surrogate_data!(
                dummy_data,
                device_details,
                data_collection_params,
                sim_full,
            )
        else
            @warn "Simulation not attempted because the ground truth scenario could not be built"
        end
        copy!(Random.default_rng(), rng_state)
    end
    ########################################################################################
    ########################################################################################

    for (ix_o, o) in enumerate(operating_points)
        for (ix_p, p) in enumerate(perturbations)
            sys = deepcopy(sys_main)
            update_operating_point!(sys, o, sys_aux)
            PSID.Simulation!(PSID.MassMatrixModel, sys, pwd(), (0.0, 0.0))  #Run power flow and re-init devices before defining perturbation
            psid_perturbations = PSID.Perturbation[]
            for p_single in p
                add_surrogate_perturbation!(sys, psid_perturbations, p_single, sys_aux)
            end
            rng_state = copy(Random.default_rng())
            data = EmptyTrainDataSet(dataset_type)
            if dataset_aux !== nothing &&
               dataset_aux[(ix_o - 1) * size(perturbations)[1] + ix_p].built == true
                match_operating_point(
                    sys,
                    dataset_aux[(ix_o - 1) * size(perturbations)[1] + ix_p],
                    surrogate_params,
                )    #TODO - not tested
                sim_full = _build_run_simulation_perturbations(
                    sys,
                    data_collection_params,
                    psid_perturbations,
                )
                if !issubset(
                    data_collection_params.tsave,
                    union(data_collection_params.tstops, sim_full.tstops),
                )
                    @warn "tsave not subset of tstops"
                end
                fill_surrogate_data!(data, device_details, data_collection_params, sim_full)
            elseif dataset_aux === nothing
                #Instead of running the simulation, just compare the error in the model 
                sim_full = _build_run_simulation_perturbations(
                    sys,
                    data_collection_params,
                    psid_perturbations,
                )
                if !issubset(
                    data_collection_params.tsave,
                    union(data_collection_params.tstops, sim_full.tstops),
                )
                    @warn "tsave not subset of tstops"
                end
                fill_surrogate_data!(data, device_details, data_collection_params, sim_full)
            else
                @warn "Simulation not attempted because the ground truth scenario could not be built"
            end
            push!(train_data, data)
            copy!(Random.default_rng(), rng_state)
        end
    end
    return train_data
end

function _build_run_simulation_perturbations(sys, data_collection, psid_perturbations)
    tspan = data_collection.tspan
    solver = instantiate_solver(data_collection)
    abstol = data_collection.solver_tols.abstol
    reltol = data_collection.solver_tols.reltol

    if data_collection.formulation == "MassMatrix"
        if data_collection.frequency_reference == "ReferenceBus"
            sim_full = PSID.Simulation!(
                PSID.MassMatrixModel,
                sys,
                pwd(),
                tspan,
                psid_perturbations;
                all_branches_dynamic = data_collection.all_branches_dynamic,
                all_lines_dynamic = data_collection.all_lines_dynamic,
                frequency_reference = PSID.ReferenceBus(),
            )
        elseif data_collection.frequency_reference == "ConstantFrequency"
            sim_full = PSID.Simulation!(
                PSID.MassMatrixModel,
                sys,
                pwd(),
                tspan,
                psid_perturbations;
                all_branches_dynamic = data_collection.all_branches_dynamic,
                all_lines_dynamic = data_collection.all_lines_dynamic,
                frequency_reference = PSID.ConstantFrequency(),
            )
        end
    elseif data_collection.formulation == "Residual"
        if data_collection.frequency_reference == "ReferenceBus"
            sim_full = PSID.Simulation!(
                PSID.ResidualModel,
                sys,
                pwd(),
                tspan,
                psid_perturbations;
                all_branches_dynamic = data_collection.all_branches_dynamic,
                all_lines_dynamic = data_collection.all_lines_dynamic,
                frequency_reference = PSID.ReferenceBus(),
            )
        elseif data_collection.frequency_reference == "ConstantFrequency"
            sim_full = PSID.Simulation!(
                PSID.ResidualModel,
                sys,
                pwd(),
                tspan,
                psid_perturbations;
                all_branches_dynamic = data_collection.all_branches_dynamic,
                all_lines_dynamic = data_collection.all_lines_dynamic,
                frequency_reference = PSID.ConstantFrequency(),
            )
        end
    end
    if sim_full.status == PSID.BUILT
        if data_collection.dtmax === nothing
            #PSID.show_states_initial_value(sim_full)   #For debugging of SS error for surrogate
            PSID.execute!(
                sim_full,
                solver,
                abstol = abstol,
                reltol = reltol,
                maxiters = 1e4,
                tstops = union(data_collection.tstops, sim_full.tstops),    #sim_full.tstops can have tstops that are required for re-initialization after a perturbation.
                save_everystep = true,
                saveat = data_collection.tsave,
                reset_simulation = false,
                enable_progress_bar = true,
            )
        else
            PSID.execute!(
                sim_full,
                solver,
                abstol = abstol,
                reltol = reltol,
                maxiters = 1e4,
                tstops = union(data_collection.tstops, sim_full.tstops),    #sim_full.tstops can have tstops that are required for re-initialization after a perturbation.
                save_everystep = true,
                saveat = data_collection.tsave,
                reset_simulation = false,
                enable_progress_bar = true,
                dtmax = data_collection.dtmax,
            )
        end
    end

    return sim_full
end

function _build_run_simulation_initial_conditions(sys, data_collection, initial_conditions)
    tspan = data_collection.tspan
    solver = instantiate_solver(data_collection)
    abstol = data_collection.solver_tols.abstol
    reltol = data_collection.solver_tols.reltol

    if data_collection.formulation == "MassMatrix"
        if data_collection.frequency_reference == "ReferenceBus"
            sim_full = PSID.Simulation!(
                PSID.MassMatrixModel,
                sys,
                pwd(),
                tspan,
                initialize_simulation = false,
                initial_conditions = initial_conditions;
                all_branches_dynamic = data_collection.all_branches_dynamic,
                all_lines_dynamic = data_collection.all_lines_dynamic,
                frequency_reference = PSID.ReferenceBus(),
            )
        elseif data_collection.frequency_reference == "ConstantFrequency"
            sim_full = PSID.Simulation!(
                PSID.MassMatrixModel,
                sys,
                pwd(),
                tspan,
                initialize_simulation = false,
                initial_conditions = initial_conditions;
                all_branches_dynamic = data_collection.all_branches_dynamic,
                all_lines_dynamic = data_collection.all_lines_dynamic,
                frequency_reference = PSID.ConstantFrequency(),
            )
        end
    elseif data_collection.formulation == "Residual"
        if data_collection.frequency_reference == "ReferenceBus"
            sim_full = PSID.Simulation!(
                PSID.ResidualModel,
                sys,
                pwd(),
                tspan,
                initialize_simulation = false,
                initial_conditions = initial_conditions;
                all_branches_dynamic = data_collection.all_branches_dynamic,
                all_lines_dynamic = data_collection.all_lines_dynamic,
                frequency_reference = PSID.ReferenceBus(),
            )
        elseif data_collection.frequency_reference == "ConstantFrequency"
            sim_full = PSID.Simulation!(
                PSID.ResidualModel,
                sys,
                pwd(),
                tspan,
                initialize_simulation = false,
                initial_conditions = initial_conditions;
                all_branches_dynamic = data_collection.all_branches_dynamic,
                all_lines_dynamic = data_collection.all_lines_dynamic,
                frequency_reference = PSID.ConstantFrequency(),
            )
        end
    end
    if data_collection.dtmax === nothing
        PSID.execute!(
            sim_full,
            solver,
            abstol = abstol,
            reltol = reltol,
            tstops = union(data_collection.tstops, sim_full.tstops),    #sim_full.tstops can have tstops that are required for re-initialization after a perturbation.
            save_everystep = true,
            saveat = data_collection.tsave,
            reset_simulation = false,
            enable_progress_bar = true,
        )
    else
        PSID.execute!(
            sim_full,
            solver,
            abstol = abstol,
            reltol = reltol,
            tstops = union(data_collection.tstops, sim_full.tstops),    #sim_full.tstops can have tstops that are required for re-initialization after a perturbation.
            save_everystep = true,
            saveat = data_collection.tsave,
            reset_simulation = false,
            enable_progress_bar = true,
            dtmax = data_collection.dtmax,
        )
    end
    return sim_full
end

function fill_surrogate_data!(
    data,
    device_details,
    data_collection_params,
    sim_full,
    results,
    save_indices,
)
    @warn "collect_data not implemented for this type of SurrogateDataSet: $(typeof(data))"
end

function _match_operating_point(sys, P0, Q0, Vm0, θ0, surrogate_params::Nothing)
    @error "Passed auxiliary dataset without passing surrogate parameters. Should only be used for tests."
end
