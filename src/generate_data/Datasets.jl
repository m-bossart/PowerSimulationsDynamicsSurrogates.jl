abstract type SurrogateTrainDataset end

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
    function generate_train_data(
        sys_train::PSY.System,
        perturbations::Vector{SurrogatePerturbation},
        operating_points::Vector{SurrogateOperatingPoint},
        data_name::String,   
        data_collection_params::GenerateDataParams,
    )
- Options for `data_name`: `"SteadyStateNODEData"`
"""
function generate_train_data(
    sys_train::PSY.System,
    perturbations,  #TODO - type info
    operating_points::Vector{O},
    data_name::String,
    data_collection_params::GenerateDataParams,
) where {O <: SurrogateOperatingPoint}
    train_data = SurrogateTrainDataset[]
    for o in operating_points
        for p in perturbations
            sys = deepcopy(sys_train)           #new copy of original system for each (o,p) combination
            update_operating_point!(sys, o)
            psid_perturbations = PSID.Perturbation[]
            for p_single in p
                add_surrogate_perturbation!(sys, psid_perturbations, p_single)
            end
            data = EmptyTrainDataSet(data_name)
            fill_surrogate_data!(data, sys, psid_perturbations, data_collection_params)
            push!(train_data, data)
        end
    end
    return train_data
end

mutable struct SteadyStateNODEData <: SurrogateTrainDataset
    tsteps::AbstractArray
    groundtruth_current::AbstractArray
    connecting_impedance::AbstractArray
    powerflow::AbstractArray
    branch_order::Array{String}
end

function SteadyStateNODEData()
    return SteadyStateNODEData([], [], [], [], [])
end

function fill_surrogate_data!(
    data::T,
    sys::PSY.System,
    psid_perturbations,
    data_collection::GenerateDataParams,
) where {T <: SurrogateTrainDataset}
    @warn "collect_data not implemented for this type of SurrogateDataSet"
end

function fill_surrogate_data!(
    data::SteadyStateNODEData,
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

    sources = collect(PSY.get_components(PSY.Source, sys_train))
    buses_with_sources = PSY.get_bus.(sources)

    #find connecting branches
    connecting_branches = String[]
    for b in PSY.get_components(PSY.ACBranch, sys_train)
        to_bus = PSY.get_to(PSY.get_arc(b))
        from_bus = PSY.get_from(PSY.get_arc(b))
        if (from_bus in buses_with_sources) || (to_bus in buses_with_sources)
            push!(connecting_branches, PSY.get_name(b))
        end
    end

    sim_full = PSID.Simulation!(
        PSID.MassMatrixModel,   #todo -set from parameters 
        sys_train,
        pwd(),
        tspan,
        psid_perturbations,
        #console_level = PSID_CONSOLE_LEVEL,
        # file_level = PSID_FILE_LEVEL,
    )
    display(sim_full)

    PSID.execute!(
        sim_full,
        solver,
        abstol = abstol,
        reltol = reltol,
        reset_simulation = false,
        saveat = tsteps,
        enable_progress_bar = false,
    )

    ground_truth_current = zeros(2 * length(connecting_branches), length(tsteps))
    connecting_impedance = zeros(length(connecting_branches), 2)
    powerflow = zeros(length(connecting_branches) * 4)

    for (i, branch_name) in enumerate(connecting_branches)
        ground_truth_current[2 * i - 1, :] = get_total_current_series(sim_full)[1, :]
        ground_truth_current[2 * i, :] = get_total_current_series(sim_full)[2, :]
        #ground_truth_ir[i, :] = get_branch_current(branch_name, :Ir) #TODO:  Get branch current instead when this issue closes: https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/issues/224
        #ground_truth_ii[i, :] = get_branch_current(branch_name, :Ii)
        connecting_impedance[i, :] =
            _get_branch_plus_source_impedance(sys_train, branch_name)
        powerflow[(i * 4 - 3):(i * 4)] =
            _get_powerflow_opposite_source(sys_train, branch_name)
    end

    data.tsteps = tsteps
    data.groundtruth_current = ground_truth_current
    data.connecting_impedance = connecting_impedance
    data.powerflow = powerflow
    data.branch_order = connecting_branches
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
    @assert length(source_active) == 1
    source_active = source_active[1]
    return [
        PSY.get_x(ac_branch) + PSY.get_X_th(source_active),
        PSY.get_r(ac_branch) + PSY.get_R_th(source_active),
    ]
end

#TODO - don't use run_powerflow, just use the first value from the dynamic simualtion. Will be much easier to keep track of 
function _get_powerflow_opposite_source(sys_train, branch_name)
    flow_results = PowerFlows.run_powerflow(sys_train)["flow_results"]
    bus_results = PowerFlows.run_powerflow(sys_train)["bus_results"]
    display(flow_results)
    display(bus_results)

    connecting_branch_results = filter(row -> row.line_name == branch_name, flow_results)
    @assert size(connecting_branch_results)[1] == 1

    bus_from = connecting_branch_results.bus_from[1]
    bus_to = connecting_branch_results.bus_to[1]

    source_active = collect(
        PSY.get_components(
            PSY.Source,
            sys_train,
            x ->
                PSY.get_available(x) &&
                    PSY.get_number(PSY.get_bus(x)) in [bus_from, bus_to],
        ),
    )
    @assert length(source_active) == 1
    source_active = source_active[1]

    if PSY.get_number(PSY.get_bus(source_active)) == bus_from
        P_pf = connecting_branch_results.P_to_from[1] / PSY.get_base_power(sys_train)
        Q_pf = connecting_branch_results.Q_to_from[1] / PSY.get_base_power(sys_train)
        opposite_source_bus_results = filter(row -> row.bus_number == bus_to, bus_results)
        V_pf = opposite_source_bus_results.Vm[1]
        θ_pf = opposite_source_bus_results.θ[1]
        return [P_pf, Q_pf, V_pf, θ_pf]
    elseif PSY.get_number(PSY.get_bus(source_active)) == bus_to
        P_pf = connecting_branch_results.P_from_to[1] / PSY.get_base_power(sys_train)
        Q_pf = connecting_branch_results.Q_from_to[1] / PSY.get_base_power(sys_train)
        opposite_source_bus_results = filter(row -> row.bus_number == bus_from, bus_results)
        V_pf = opposite_source_bus_results.Vm[1]
        θ_pf = opposite_source_bus_results.θ[1]
        return [P_pf, Q_pf, V_pf, θ_pf]
    else
        @error "Didn't find a source at the correct bus"
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
    )
    return d[key]
end

#TODO - replace with the branch or source function for current! 
function get_total_current_series(sim)
    ir_total = []
    ii_total = []
    for (i, g) in enumerate(
        PSY.get_components(
            PSY.DynamicInjection,
            sim.sys,
            x -> typeof(x) !== PSY.PeriodicVariableSource,
        ),
    )
        results = PSID.read_results(sim)
        if i == 1
            ir_total = PSID.get_real_current_series(results, PSY.get_name(g))[2]
            ii_total = PSID.get_imaginary_current_series(results, PSY.get_name(g))[2]
        else
            ir_total .+= PSID.get_real_current_series(results, PSY.get_name(g))[2]
            ii_total .+= PSID.get_imaginary_current_series(results, PSY.get_name(g))[2]
        end
    end
    data_array = zeros(Float64, (2, length(ir_total)))
    data_array[1, :] .= ir_total
    data_array[2, :] .= ii_total
    return data_array
end

function EmptyTrainDataSet(key)
    d = Dict("SteadyStateNODEData" => SteadyStateNODEData)
    return d[key]()
end
