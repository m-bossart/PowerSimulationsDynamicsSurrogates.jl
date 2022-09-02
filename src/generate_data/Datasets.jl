abstract type SurrogateDataset end
abstract type SurrogateDatasetParams end

struct GenerateDataParams
    solver::String
    solver_tols::Tuple{Float64, Float64}
    tspan::Tuple{Float64, Float64}
    steps::Int64
    tsteps_spacing::String
    formulation::String
    all_branches_dynamic::Bool
    all_lines_dynamic::Bool
    seed::Int64
end

function GenerateDataParams(;
    solver = "IDA",
    solver_tols = (1e-6, 1e-6),
    tspan = (0.0, 1.0),
    steps = 100,
    tsteps_spacing = "linear",
    formulation = "Residual",
    all_branches_dynamic = false,
    all_lines_dynamic = false,
    seed = 1,
)
    GenerateDataParams(
        solver,
        solver_tols,
        tspan,
        steps,
        tsteps_spacing,
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
    sys_aux::PSY.System,
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
    connecting_branch_names::Array{Tuple{String, Symbol}}
end

function SteadyStateNODEDataParams(;
    type = "SteadyStateNODEDataParams",
    connecting_branch_names = [],
)
    return SteadyStateNODEDataParams(type, connecting_branch_names)
end

mutable struct SteadyStateNODEData <: SurrogateDataset
    type::String
    tsteps::AbstractArray
    groundtruth_current::AbstractArray
    groundtruth_voltage::AbstractArray
    connecting_impedance::AbstractArray
    powerflow::AbstractArray
    stable::Bool
end

function SteadyStateNODEData(;
    type = "SteadyStateNODEData",
    tsteps = [],
    groundtruth_current = [],
    groundtruth_voltage = [],
    connecting_impedance = [],
    powerflow = [],
    stable = false,
)
    return SteadyStateNODEData(
        type,
        tsteps,
        groundtruth_current,
        groundtruth_voltage,
        connecting_impedance,
        powerflow,
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
    steps = data_collection.steps
    if data_collection.tsteps_spacing == "linear"
        tsteps = tspan[1]:((tspan[2] - tspan[1]) / steps):tspan[2]
    end
    solver = instantiate_solver(data_collection)
    abstol = data_collection.solver_tols[2]
    reltol = data_collection.solver_tols[1]
    #display(PSY.solve_powerflow(sys_train)["bus_results"])
    #display(PSY.solve_powerflow(sys_train)["flow_results"])
    if data_aux !== nothing
        @warn "Updating the operating point based on previously run dataset"
        for s in PSY.get_components(
            PSY.Source,
            sys_train,
            x -> typeof(PSY.get_dynamic_injector(x)) == SteadyStateNODE,
        )
            PSY.set_active_power!(s, data_aux.powerflow[1])
            PSY.set_reactive_power!(s, data_aux.powerflow[2])
            PSY.set_internal_voltage!(s, data_aux.powerflow[3])
            PSY.set_internal_angle!(s, data_aux.powerflow[4])
        end
    end
    #display(PSY.solve_powerflow(sys_train)["bus_results"])
    #display(PSY.solve_powerflow(sys_train)["flow_results"])
    connecting_branches = params.connecting_branch_names
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
    PSID.execute!(
        sim_full,
        solver,
        abstol = abstol,
        reltol = reltol,
        reset_simulation = false,
        saveat = tsteps,
        enable_progress_bar = false,
    )
    results = PSID.read_results(sim_full)
    if sim_full.status == PSID.SIMULATION_FINALIZED
        ground_truth_current = zeros(length(connecting_branches) * 2, length(tsteps))
        ground_truth_voltage = zeros(length(connecting_branches) * 2, length(tsteps))
        connecting_impedance = zeros(length(connecting_branches), 2)
        powerflow = zeros(length(connecting_branches) * 4)
        for (i, branch_tuple) in enumerate(connecting_branches)
            branch_name = branch_tuple[1]
            location = branch_tuple[2]
            if location == :from
                ground_truth_current[2 * i - 1, :] =
                    PSID.get_real_current_branch_flow(results, branch_name)[2]
                ground_truth_current[2 * i, :] =
                    PSID.get_imaginary_current_branch_flow(results, branch_name)[2]
                P0 = PSID.get_activepower_branch_flow(results, branch_name, location)[2][1]
                Q0 =
                    PSID.get_reactivepower_branch_flow(results, branch_name, location)[2][1]
                branch = PSY.get_component(PSY.ACBranch, sys_train, branch_name)
                bus_number = PSY.get_number(PSY.get_from(PSY.get_arc(branch)))
                V = PSID.get_voltage_magnitude_series(results, bus_number)[2]
                ground_truth_voltage[2 * i - 1, :] = V
                V0 = V[1]
                θ = PSID.get_voltage_angle_series(results, bus_number)[2]
                ground_truth_voltage[2 * i, :] = θ
                θ0 = θ[1]
            elseif location == :to
                ground_truth_current[2 * i - 1, :] =
                    PSID.get_real_current_branch_flow(results, branch_name)[2] * -1
                ground_truth_current[2 * i, :] =
                    PSID.get_imaginary_current_branch_flow(results, branch_name)[2] * -1
                P0 =
                    PSID.get_activepower_branch_flow(results, branch_name, location)[2][1] *
                    -1
                Q0 =
                    PSID.get_reactivepower_branch_flow(results, branch_name, location)[2][1] *
                    -1
                branch = PSY.get_component(PSY.ACBranch, sys_train, branch_name)
                bus_number = PSY.get_number(PSY.get_to(PSY.get_arc(branch)))
                V = PSID.get_voltage_magnitude_series(results, bus_number)[2]
                ground_truth_voltage[2 * i - 1, :] = V
                V0 = V[1]
                θ = PSID.get_voltage_angle_series(results, bus_number)[2]
                ground_truth_voltage[2 * i, :] = θ
                θ0 = θ[1]
            end
            connecting_impedance[i, :] =
                _get_branch_plus_source_impedance(sys_train, branch_tuple[1])
            powerflow[(i * 4 - 3):(i * 4)] = [P0, Q0, V0, θ0]
        end
        data.tsteps = tsteps
        data.groundtruth_current = ground_truth_current
        data.groundtruth_voltage = ground_truth_voltage
        data.connecting_impedance = connecting_impedance
        data.powerflow = powerflow
        data.stable = true
    end
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
