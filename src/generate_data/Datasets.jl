abstract type SurrogateDataset end
abstract type SurrogateDatasetParams end

struct GenerateDataParams
    solver::String
    solver_tols::Tuple{Float64, Float64}
    tspan::Tuple{Float64, Float64}
    steps::Int64
    tsteps_spacing::String
    formulation::String
end

function GenerateDataParams(;
    solver = "Rodas4",
    solver_tols = (1e-6, 1e-6),
    tspan = (0.0, 1.0),
    steps = 100,
    tsteps_spacing = "linear",
    formulation = "MassMatrix",
)
    GenerateDataParams(solver, solver_tols, tspan, steps, tsteps_spacing, formulation)
end

"""
    function generate_surrogate_data(
        sys_main::PSY.System,
        sys_aux::PSY.System,          
        perturbations::Vector{Vector{Union{SurrogatePerturbation, PSID.Perturbation}}},
        operating_points::Vector{O},
        data_params::D,
        data_collection_params::GenerateDataParams,
    ) where {O <: SurrogateOperatingPoint, D <: SurrogateDatasetParams}

- `sys_main` is the system that is changed and perturbed and used to generate data. \\
- `sys_aux` is used by some operating points and perturbations which require randomly choosing components to perturb. \\
- This function creates a new deepcopy of `sys_main` for each pair of (perturbations, operating_points)
"""
function generate_surrogate_data(
    sys_main::PSY.System,
    sys_aux::PSY.System,
    perturbations::Vector{Vector{Union{SurrogatePerturbation, PSID.Perturbation}}},
    operating_points::Vector{O},
    data_params::D,
    data_collection_params::GenerateDataParams;
    seed::Int64 = 1,
    stable_trajectories::BitVector = trues(
        size(perturbations)[1] * length(operating_points),
    ),
) where {O <: SurrogateOperatingPoint, D <: SurrogateDatasetParams}
    Random.seed!(seed)
    stable_trajectories_out = copy(stable_trajectories)
    @assert length(stable_trajectories) == size(perturbations)[1] * length(operating_points)
    train_data = SurrogateDataset[]
    for (ix_o, o) in enumerate(operating_points)
        for (ix_p, p) in enumerate(perturbations)
            sys = deepcopy(sys_main)
            update_operating_point!(sys, o, sys_aux)
            psid_perturbations = PSID.Perturbation[]
            for p_single in p
                add_surrogate_perturbation!(sys, psid_perturbations, p_single, sys_aux)
            end
            if stable_trajectories[(ix_o - 1) * size(perturbations)[1] + ix_p] == true
                data = EmptyTrainDataSet(data_params)
                retcode = fill_surrogate_data!(
                    data,
                    data_params,
                    sys,
                    psid_perturbations,
                    data_collection_params,
                )
                if retcode == :Success
                    push!(train_data, data)
                else
                    stable_trajectories_out[(ix_o - 1) * size(perturbations)[1] + ix_p] =
                        false
                end
            end
        end
    end
    return train_data, stable_trajectories_out
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
    connecting_impedance::AbstractArray
    powerflow::AbstractArray
end

function SteadyStateNODEData(;
    type = "SteadyStateNODEData",
    tsteps = [],
    groundtruth_current = [],
    connecting_impedance = [],
    powerflow = [],
)
    return SteadyStateNODEData(
        type,
        tsteps,
        groundtruth_current,
        connecting_impedance,
        powerflow,
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
)
    tspan = data_collection.tspan
    steps = data_collection.steps
    if data_collection.tsteps_spacing == "linear"
        tsteps = tspan[1]:((tspan[2] - tspan[1]) / steps):tspan[2]
    end
    solver = instantiate_solver(data_collection)
    abstol = data_collection.solver_tols[2]
    reltol = data_collection.solver_tols[1]

    connecting_branches = params.connecting_branch_names
    if data_collection.formulation == "MassMatrix"
        sim_full = PSID.Simulation!(
            PSID.MassMatrixModel,
            sys_train,
            pwd(),
            tspan,
            psid_perturbations,
        )
    elseif data_collection.formulation == "Residual"
        sim_full = PSID.Simulation!(
            PSID.ResidualModel,
            sys_train,
            pwd(),
            tspan,
            psid_perturbations,
        )
    end
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
                V0 = PSID.get_voltage_magnitude_series(results, bus_number)[2][1]
                θ0 = PSID.get_voltage_angle_series(results, bus_number)[2][1]
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
                V0 = PSID.get_voltage_magnitude_series(results, bus_number)[2][1]
                θ0 = PSID.get_voltage_angle_series(results, bus_number)[2][1]
            end
            connecting_impedance[i, :] =
                _get_branch_plus_source_impedance(sys_train, branch_tuple[1])
            powerflow[(i * 4 - 3):(i * 4)] = [P0, Q0, V0, θ0]
        end
        data.tsteps = tsteps
        data.groundtruth_current = ground_truth_current
        data.connecting_impedance = connecting_impedance
        data.powerflow = powerflow
        return results.solution.retcode
    else
        return :Failed
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
        "TRBDF2" => OrdinaryDiffEq.TRBDF2,
        "Tsit5" => OrdinaryDiffEq.Tsit5,
        "IDA" => Sundials.IDA,
    )
    return d[key]
end

function EmptyTrainDataSet(key::SteadyStateNODEDataParams)
    return SteadyStateNODEData()
end
