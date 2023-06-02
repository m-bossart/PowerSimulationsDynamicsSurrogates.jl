abstract type SurrogateDataset end

struct GenerateDataParams
    solver::String
    solver_tols::NamedTuple{(:reltol, :abstol), Tuple{Float64, Float64}}
    tspan::Tuple{Float64, Float64}
    tstops::Vector{Float64}
    tsave::Vector{Float64}
    formulation::String
    all_branches_dynamic::Bool
    all_lines_dynamic::Bool
    seed::Int64
end

function GenerateDataParams(;
    solver = "IDA",
    solver_tols = (reltol = 1e-6, abstol = 1e-6),
    tspan = (0.0, 1.0),
    tstops = [],
    tsave = [],
    formulation = "Residual",
    all_branches_dynamic = false,
    all_lines_dynamic = false,
    seed = 1,
)
    GenerateDataParams(
        solver,
        solver_tols,
        tspan,
        tstops,
        tsave,
        formulation,
        all_branches_dynamic,
        all_lines_dynamic,
        seed,
    )
end

"""
    function generate_surrogate_data(
        sys_main::PSY.System,
        ics::Vector{Vector{Float64}},
        data_params::D,
        data_collection_params::GenerateDataParams,
        dataset_aux = nothing,
    ) where {O <: SurrogateOperatingPoint, D <: SurrogateDatasetParams}

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
    ics::Vector{Vector{Float64}},
    device_details::Dict{String, Dict{Symbol, Symbol}},
    data_collection_params::GenerateDataParams;
    dataset_aux = nothing,
    surrogate_params = nothing,
)  where { T <: SurrogateDataset}
    Random.seed!(data_collection_params.seed)
    train_data = dataset_type[]
    ########################################################################################
    #Run full simulation for first operating point and perturbation in order to precompile
    #all code and avoid including in timing.
    ########################################################################################
    sys = deepcopy(sys_main)
    dummy_data = EmptyTrainDataSet(dataset_type)
    sim_full = _build_run_simulation_initial_conditions(sys, data_collection_params, ics[1])
    @assert issubset(data_collection_params.tsave, union(data_collection_params.tstops, sim_full.tstops))
    if sim_full.status == PSID.SIMULATION_FINALIZED
        results = PSID.read_results(sim_full)
        if length(data_collection_params.tsave) == 0
            save_indices = 1:length(unique(results.solution.t))
        else
            save_indices = indexin(data_collection_params.tsave, unique(results.solution.t))
        end
        fill_surrogate_data!(
            dummy_data,
            device_details,
            data_collection_params, 
            sim_full,
            results,
            save_indices,
        )
    end 
    ########################################################################################
    ########################################################################################
    for ic in ics
        sys = deepcopy(sys_main)
        #PowerFlows.run_powerflow!(sys)
        data = EmptyTrainDataSet(dataset_type)
        sim_full = _build_run_simulation_initial_conditions(sys, data_collection_params, ic)
        @assert issubset(data_collection_params.tsave, union(data_collection_params.tstops, sim_full.tstops))
        if sim_full.status == PSID.SIMULATION_FINALIZED
            results = PSID.read_results(sim_full)
            if length(data_collection_params.tsave) == 0
                save_indices = 1:length(unique(results.solution.t))
            else
                save_indices = indexin(data_collection_params.tsave, unique(results.solution.t))
            end
            fill_surrogate_data!(
                data,
                device_details,
                data_collection_params, 
                sim_full,
                results,
                save_indices,
            )
        end 
        push!(train_data, data)
    end
    return train_data
 end 


"""
    function generate_surrogate_data(
        sys_main::PSY.System,
        sys_aux::PSY.System,          
        perturbations::Vector{Vector{Union{SurrogatePerturbation, PSID.Perturbation}}},
        operating_points::Vector{O},
        data_params::D,
        data_collection_params::GenerateDataParams,
        dataset_aux = nothing,
    ) where {O <: SurrogateOperatingPoint, D <: SurrogateDatasetParams}

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
) where {O <: SurrogateOperatingPoint, T <: SurrogateDataset}
    Random.seed!(data_collection_params.seed)
    train_data = dataset_type[]
    ########################################################################################
    #Run full simulation for first operating point and perturbation in order to precompile
    #all code and avoid including in timing.
    ########################################################################################
    sys = deepcopy(sys_main)
    update_operating_point!(sys, operating_points[1], sys_aux)
    psid_perturbations = PSID.Perturbation[]
    for p_single in perturbations[1]
        add_surrogate_perturbation!(sys, psid_perturbations, p_single, sys_aux)
    end
    dummy_data = EmptyTrainDataSet(dataset_type)
    if dataset_aux !== nothing
        match_operating_point(sys, dataset_aux[1], surrogate_params)    #TODO - not tested
    end
    sim_full = _build_run_simulation_perturbations(sys, data_collection_params, psid_perturbations)
    @assert issubset(data_collection_params.tsave, union(data_collection_params.tstops, sim_full.tstops))
    if sim_full.status == PSID.SIMULATION_FINALIZED
        results = PSID.read_results(sim_full)
        if length(data_collection_params.tsave) == 0
            save_indices = 1:length(unique(results.solution.t))
        else
            save_indices = indexin(data_collection_params.tsave, unique(results.solution.t))
        end
        fill_surrogate_data!(
            dummy_data,
            device_details,
            data_collection_params, 
            sim_full,
            results,
            save_indices,
        )
    end 
    ########################################################################################
    ########################################################################################

    for (ix_o, o) in enumerate(operating_points)
        for (ix_p, p) in enumerate(perturbations)
            sys = deepcopy(sys_main)
            update_operating_point!(sys, o, sys_aux)
            psid_perturbations = PSID.Perturbation[]
            for p_single in p
                add_surrogate_perturbation!(sys, psid_perturbations, p_single, sys_aux)
            end
            data = EmptyTrainDataSet(dataset_type)
            if dataset_aux !== nothing
                match_operating_point(sys, dataset_aux[(ix_o -1) * size(perturbations)[1] + ix_p], surrogate_params)    #TODO - not tested
            end
            sim_full = _build_run_simulation_perturbations(sys, data_collection_params, psid_perturbations)
            @assert issubset(data_collection_params.tsave, union(data_collection_params.tstops, sim_full.tstops))
            if sim_full.status == PSID.SIMULATION_FINALIZED
                results = PSID.read_results(sim_full)
                if length(data_collection_params.tsave) == 0
                    save_indices = 1:length(unique(results.solution.t))
                else
                    save_indices = indexin(data_collection_params.tsave, unique(results.solution.t))
                end
                fill_surrogate_data!(
                    data,
                    device_details,
                    data_collection_params, 
                    sim_full,
                    results,
                    save_indices,
                )
            end 
            push!(train_data, data)
        end
    end
    return train_data
end


function _build_run_simulation_perturbations(
    sys,
    data_collection,
    psid_perturbations,
)
    tspan = data_collection.tspan
    solver = instantiate_solver(data_collection)
    abstol = data_collection.solver_tols.abstol
    reltol = data_collection.solver_tols.reltol

    if data_collection.formulation == "MassMatrix"
        sim_full = PSID.Simulation(
            PSID.MassMatrixModel,
            sys,
            pwd(),
            tspan,
            psid_perturbations;
            all_branches_dynamic = data_collection.all_branches_dynamic,
            all_lines_dynamic = data_collection.all_lines_dynamic,
        )
    elseif data_collection.formulation == "Residual"
        sim_full = PSID.Simulation(
            PSID.ResidualModel,
            sys,
            pwd(),
            tspan,
            psid_perturbations;
            all_branches_dynamic = data_collection.all_branches_dynamic,
            all_lines_dynamic = data_collection.all_lines_dynamic,
        )
    end
    PSID.execute!(
        sim_full,
        solver,
        abstol = abstol,
        reltol = reltol,
        tstops = union(data_collection.tstops, sim_full.tstops),    #sim_full.tstops can have tstops that are required for re-initialization after a perturbation.
        save_everystep = true,
        saveat = data_collection.tsave,
        reset_simulation = false,
        enable_progress_bar = false,
    )
    return sim_full 
end 

function _build_run_simulation_initial_conditions(
    sys,
    data_collection,
    initial_conditions,
)
    tspan = data_collection.tspan
    solver = instantiate_solver(data_collection)
    abstol = data_collection.solver_tols.abstol
    reltol = data_collection.solver_tols.reltol

    if data_collection.formulation == "MassMatrix"
        sim_full = PSID.Simulation!(
            PSID.MassMatrixModel,
            sys,
            pwd(),
            tspan,
            initialize_simulation = false,
            initial_conditions = initial_conditions;
            all_branches_dynamic = data_collection.all_branches_dynamic,
            all_lines_dynamic = data_collection.all_lines_dynamic,
        )
        ic = PSID.get_initial_conditions(sim_full)
    elseif data_collection.formulation == "Residual"
        sim_full = PSID.Simulation!(
            PSID.ResidualModel,
            sys,
            pwd(),
            tspan,
            initialize_simulation = false,
            initial_conditions = initial_conditions;
            all_branches_dynamic = data_collection.all_branches_dynamic,
            all_lines_dynamic = data_collection.all_lines_dynamic,
        )
        ic = PSID.get_initial_conditions(sim_full)
    end
    PSID.execute!(
        sim_full,
        solver,
        abstol = abstol,
        reltol = reltol,
        tstops = union(data_collection.tstops, sim_full.tstops),    #sim_full.tstops can have tstops that are required for re-initialization after a perturbation.
        save_everystep = true,
        saveat = data_collection.tsave,
        reset_simulation = false,
        enable_progress_bar = false,
    )
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
