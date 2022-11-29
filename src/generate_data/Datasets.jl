abstract type SurrogateDataset end
abstract type SurrogateDatasetParams end

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
"""
function generate_surrogate_data(
    sys_main::PSY.System,
    sys_aux::PSY.System,    #TODO make this optional
    perturbations::Vector{Vector{Union{SurrogatePerturbation, PSID.Perturbation}}},
    operating_points::Vector{O},
    data_params::D,
    data_collection_params::GenerateDataParams;
    dataset_aux = nothing,
) where {O <: SurrogateOperatingPoint, D <: SurrogateDatasetParams}
    Random.seed!(data_collection_params.seed)
    #@assert length(stable_trajectories) == size(perturbations)[1] * length(operating_points)
    train_data = typeof(EmptyTrainDataSet(data_params))[]
    for (ix_o, o) in enumerate(operating_points)
        for (ix_p, p) in enumerate(perturbations)
            sys = deepcopy(sys_main)
            update_operating_point!(sys, o, sys_aux)
            psid_perturbations = PSID.Perturbation[]
            for p_single in p
                add_surrogate_perturbation!(sys, psid_perturbations, p_single, sys_aux)
            end
            data = EmptyTrainDataSet(data_params)
            if dataset_aux === nothing
                fill_surrogate_data!(
                    data,
                    data_params,
                    sys,
                    psid_perturbations,
                    data_collection_params,
                    nothing,
                )
            else
                fill_surrogate_data!(
                    data,
                    data_params,
                    sys,
                    psid_perturbations,
                    data_collection_params,
                    dataset_aux[(ix_o - 1) * size(perturbations)[1] + ix_p],
                )
            end
            push!(train_data, data)
        end
    end
    return train_data
end

mutable struct SteadyStateNODEDataParams <: SurrogateDatasetParams
    type::String
    location_of_data_collection::Array{Tuple{String, Symbol}}
end

function SteadyStateNODEDataParams(;
    type = "SteadyStateNODEDataParams",
    location_of_data_collection = [],           #TODO - This is [("branchname", :to)] [("branchname", :from)]  or  --> extend to allow for [("sourcename", :source)]   
)
    return SteadyStateNODEDataParams(type, location_of_data_collection)
end

mutable struct SteadyStateNODEData <: SurrogateDataset
    type::String
    tsteps::AbstractArray
    real_current::AbstractArray
    imag_current::AbstractArray
    connecting_resistance::AbstractArray        #Not used? Leave for now 
    connecting_reactance::AbstractArray
    surrogate_real_voltage::AbstractArray
    surrogate_imag_voltage::AbstractArray
    tstops::AbstractArray
    stable::Bool
end

function SteadyStateNODEData(;
    type = "SteadyStateNODEData",
    tsteps = [],
    real_current = [],
    imag_current = [],
    connecting_resistance = [],
    connecting_reactance = [],
    surrogate_real_voltage = [],
    surrogate_imag_voltage = [],
    tstops = [],
    stable = false,
)
    return SteadyStateNODEData(
        type,
        tsteps,
        real_current,
        imag_current,
        connecting_resistance,
        connecting_reactance,
        surrogate_real_voltage,
        surrogate_imag_voltage,
        tstops,
        stable,
    )
end

function fill_surrogate_data!(
    data::T,
    params::P,
    sys::PSY.System,
    psid_perturbations,
    data_collection::GenerateDataParams,
) where {T <: SurrogateDataset, P <: SurrogateDatasetParams}
    @warn "collect_data not implemented for this type of SurrogateDataSet"
end

function fill_surrogate_data!(
    data::SteadyStateNODEData,
    params::SteadyStateNODEDataParams,
    sys_train::PSY.System,
    psid_perturbations,
    data_collection::GenerateDataParams,
    data_aux::Union{SteadyStateNODEData, Nothing},
)
    tspan = data_collection.tspan
    solver = instantiate_solver(data_collection)
    abstol = data_collection.solver_tols.abstol
    reltol = data_collection.solver_tols.reltol
    #display(PSY.solve_powerflow(sys_train)["bus_results"])
    #display(PSY.solve_powerflow(sys_train)["flow_results"])
    if data_aux !== nothing
        for s in PSY.get_components(
            PSY.Source,
            sys_train,
            x -> typeof(PSY.get_dynamic_injector(x)) == SteadyStateNODE,
        )
            Vr0 = data_aux.surrogate_real_voltage[1]
            Vi0 = data_aux.surrogate_imag_voltage[1]
            Ir0 = data_aux.real_current[1]
            Ii0 = data_aux.imag_current[1]
            P0 = Vr0 * Ir0 + Vi0 * Ii0
            Q0 = Vi0 * Ir0 - Vr0 * Ii0
            Vm0 = sqrt(Vr0^2 + Vi0^2)
            θ0 = atan(Vi0, Vr0)
            PSY.set_active_power!(s, P0)
            PSY.set_reactive_power!(s, Q0)
            PSY.set_internal_voltage!(s, Vm0)
            PSY.set_internal_angle!(s, θ0)
        end
    end
    #display(PSY.solve_powerflow(sys_train)["bus_results"])
    #display(PSY.solve_powerflow(sys_train)["flow_results"])
    if data_collection.formulation == "MassMatrix"
        sim_full = PSID.Simulation(
            PSID.MassMatrixModel,
            sys_train,
            pwd(),
            tspan,
            psid_perturbations;
            all_branches_dynamic = data_collection.all_branches_dynamic,
            all_lines_dynamic = data_collection.all_lines_dynamic,
        )
    elseif data_collection.formulation == "Residual"
        sim_full = PSID.Simulation(
            PSID.ResidualModel,
            sys_train,
            pwd(),
            tspan,
            psid_perturbations;
            all_branches_dynamic = data_collection.all_branches_dynamic,
            all_lines_dynamic = data_collection.all_lines_dynamic,
        )
    end
    #=     ss = PSID.small_signal_analysis(sim_full)       #state with index x not found in the global index 
        if !(ss.stable)
            display(ss.eigenvalues)
            @error "system not small signal stable" 
        end   =#
    #Check if the simulation was built properly? What if it wasn't able to initialize? 
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
    @assert issubset(data_collection.tsave, union(data_collection.tstops, sim_full.tstops))

    if sim_full.status == PSID.SIMULATION_FINALIZED
        results = PSID.read_results(sim_full)
        if length(data_collection.tsave) == 0
            save_indices = 1:length(unique(results.solution.t))
        else
            save_indices = indexin(data_collection.tsave, unique(results.solution.t))
        end
        results = PSID.read_results(sim_full)
        location_data_collection = params.location_of_data_collection

        if location_data_collection[1][2] == :from || location_data_collection[1][2] == :to
            _fill_data_branch!(
                data,
                results,
                location_data_collection,
                save_indices,
                sys_train,
                data_collection,
            )
        elseif location_data_collection[1][2] == :source
            _fill_data_source!(
                data,
                results,
                location_data_collection,
                save_indices,
                sys_train,
            )   #TODO - test this function which collects data from a source! 
        else
            @error "Invaled entry for data collection location"
        end
    end
end

function _fill_data_source!(data, results, connecting_sources, save_indices, sys_train)
    n_save_points = length(save_indices)
    real_current = zeros(length(connecting_sources), n_save_points)
    imag_current = zeros(length(connecting_sources), n_save_points)
    connecting_resistance = zeros(length(connecting_sources))
    connecting_reactance = zeros(length(connecting_sources))
    surrogate_real_voltage = zeros(length(connecting_sources), n_save_points)
    surrogate_imag_voltage = zeros(length(connecting_sources), n_save_points)
    for (i, source_tuple) in enumerate(connecting_sources)
        source_name = source_tuple[1]
        source = PSY.get_component(PSY.Source, sys_train, source_name)
        dyn_source = PSY.get_dynamic_injector(source)
        #Note: currents multiplied by -1 to get current out of the surrogate (into the perturbing source)
        if dyn_source === nothing
            real_current[i, :] =
                -1 .*
                PSID.get_source_real_current_series(results, source_name)[2][save_indices]
            imag_current[i, :] =
                -1 .*
                PSID.get_source_imaginary_current_series(results, source_name)[2][save_indices]
        else
            real_current[i, :] =
                -1 .* PSID.get_real_current_series(results, source_name)[2][save_indices]
            imag_current[i, :] =
                -1 .*
                PSID.get_imaginary_current_series(results, source_name)[2][save_indices]
        end
        bus_number_source = PSY.get_number(PSY.get_bus(source))
        V_surrogate =
            PSID.get_voltage_magnitude_series(results, bus_number_source)[2][save_indices]
        θ_surrogate =
            PSID.get_voltage_angle_series(results, bus_number_source)[2][save_indices]
        surrogate_real_voltage[i, :] = V_surrogate .* cos.(θ_surrogate)
        surrogate_imag_voltage[i, :] = V_surrogate .* sin.(θ_surrogate)
        connecting_reactance[i] = 0.0   #TODO - should this include the source resistance??
        connecting_resistance[i] = 0.0   #TODO - should this include the source resistance??--> not used anyways? 
    end
    data.tstops = unique(results.solution.t)
    data.tsteps = unique(results.solution.t)[save_indices]
    data.real_current = real_current
    data.imag_current = imag_current
    data.connecting_resistance = connecting_resistance
    data.connecting_reactance = connecting_reactance
    data.surrogate_real_voltage = surrogate_real_voltage
    data.surrogate_imag_voltage = surrogate_imag_voltage
    data.stable = true
end

function _fill_data_branch!(
    data,
    results,
    connecting_branches,
    save_indices,
    sys_train,
    data_collection,
)
    n_save_points = length(save_indices)
    real_current = zeros(length(connecting_branches), n_save_points)
    imag_current = zeros(length(connecting_branches), n_save_points)
    connecting_resistance = zeros(length(connecting_branches))
    connecting_reactance = zeros(length(connecting_branches))
    surrogate_real_voltage = zeros(length(connecting_branches), n_save_points)
    surrogate_imag_voltage = zeros(length(connecting_branches), n_save_points)
    for (i, branch_tuple) in enumerate(connecting_branches)
        branch_name = branch_tuple[1]
        surrogate_location = branch_tuple[2]
        #TODO - get rid of the if/else logic below after this PSID issue is resolve: https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/issues/283
        if data_collection.all_lines_dynamic
            Ir_from_to =
                PSID.get_state_series(results, (branch_name, :Il_R))[2][save_indices]
            Ii_from_to =
                PSID.get_state_series(results, (branch_name, :Il_I))[2][save_indices]
        else
            Ir_from_to =
                PSID.get_real_current_branch_flow(results, branch_name)[2][save_indices]
            Ii_from_to =
                PSID.get_imaginary_current_branch_flow(results, branch_name)[2][save_indices]
        end
        branch = PSY.get_component(PSY.ACBranch, sys_train, branch_name)
        if surrogate_location == :from
            real_current[i, :] = Ir_from_to
            imag_current[i, :] = Ii_from_to
            bus_number_surrogate = PSY.get_number(PSY.get_from(PSY.get_arc(branch)))
            bus_number_opposite = PSY.get_number(PSY.get_to(PSY.get_arc(branch)))
            V_surrogate =
                PSID.get_voltage_magnitude_series(results, bus_number_surrogate)[2][save_indices]
            θ_surrogate =
                PSID.get_voltage_angle_series(results, bus_number_surrogate)[2][save_indices]
            Vr_surrogate = V_surrogate .* cos.(θ_surrogate)
            Vi_surrogate = V_surrogate .* sin.(θ_surrogate)

        elseif surrogate_location == :to
            real_current[i, :] = Ir_from_to * -1
            imag_current[i, :] = Ii_from_to * -1
            bus_number_surrogate = PSY.get_number(PSY.get_to(PSY.get_arc(branch)))
            bus_number_opposite = PSY.get_number(PSY.get_from(PSY.get_arc(branch)))
            V_surrogate =
                PSID.get_voltage_magnitude_series(results, bus_number_surrogate)[2][save_indices]
            θ_surrogate =
                PSID.get_voltage_angle_series(results, bus_number_surrogate)[2][save_indices]
            Vr_surrogate = V_surrogate .* cos.(θ_surrogate)
            Vi_surrogate = V_surrogate .* sin.(θ_surrogate)
        end
        surrogate_real_voltage[i, :] = Vr_surrogate
        surrogate_imag_voltage[i, :] = Vi_surrogate
        connecting_resistance[i] =
            _get_branch_plus_source_impedance(sys_train, branch_tuple[1])[2]
        connecting_reactance[i] =
            _get_branch_plus_source_impedance(sys_train, branch_tuple[1])[1]
    end
    data.tstops = unique(results.solution.t)
    data.tsteps = unique(results.solution.t)[save_indices]
    data.real_current = real_current
    data.imag_current = imag_current
    data.connecting_resistance = connecting_resistance
    data.connecting_reactance = connecting_reactance
    data.surrogate_real_voltage = surrogate_real_voltage
    data.surrogate_imag_voltage = surrogate_imag_voltage
    data.stable = true
end

function _get_branch_plus_source_impedance(sys_train, branch_name)
    @assert length(PSY.get_components_by_name(PSY.ACBranch, sys_train, branch_name)) == 1
    ac_branch = PSY.get_components_by_name(PSY.ACBranch, sys_train, branch_name)[1]
    bus_from = PSY.get_from(PSY.get_arc(ac_branch))
    bus_to = PSY.get_to(PSY.get_arc(ac_branch))
    source_active = collect(
        PSY.get_components(
            PSY.Source,
            sys_train,
            x -> PSY.get_available(x) && PSY.get_bus(x) in [bus_from, bus_to],
        ),
    )
    if length(source_active) == 1
        source_active = source_active[1]
        return [
            PSY.get_x(ac_branch) + PSY.get_X_th(source_active),
            PSY.get_r(ac_branch) + PSY.get_R_th(source_active),
        ]
    else
        return [PSY.get_x(ac_branch), PSY.get_r(ac_branch)]
    end
end

function instantiate_solver(inputs)
    return solver_map(inputs.solver)()
end

function solver_map(key)
    d = Dict(
        "Rodas4" => OrdinaryDiffEq.Rodas4,
        "Rodas5" => OrdinaryDiffEq.Rodas5,
        "TRBDF2" => OrdinaryDiffEq.TRBDF2,
        "Tsit5" => OrdinaryDiffEq.Tsit5,
        "IDA" => Sundials.IDA,
    )
    return d[key]
end

function EmptyTrainDataSet(key::SteadyStateNODEDataParams)
    return SteadyStateNODEData()
end
