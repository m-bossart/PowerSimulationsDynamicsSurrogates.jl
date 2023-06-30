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
    surrogate_params = nothing,
) where {O <: SurrogateOperatingPoint, D <: SurrogateDatasetParams}
    Random.seed!(data_collection_params.seed)
    #@assert length(stable_trajectories) == size(perturbations)[1] * length(operating_points)
    train_data = typeof(EmptyTrainDataSet(data_params))[]
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
    dummy_data = EmptyTrainDataSet(data_params)
    if dataset_aux === nothing
        fill_surrogate_data!(
            dummy_data,
            data_params,
            sys,
            psid_perturbations,
            data_collection_params,
            nothing,
            nothing,
        )
    else
        fill_surrogate_data!(
            dummy_data,
            data_params,
            sys,
            psid_perturbations,
            data_collection_params,
            dataset_aux[1],
            surrogate_params,
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
            data = EmptyTrainDataSet(data_params)
            if dataset_aux === nothing
                fill_surrogate_data!(
                    data,
                    data_params,
                    sys,
                    psid_perturbations,
                    data_collection_params,
                    nothing,
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
                    surrogate_params,
                )
            end
            push!(train_data, data)
        end
    end
    return train_data
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

##################################################################
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
    ic::Dict{Symbol, Float64}
    real_current::AbstractArray
    imag_current::AbstractArray
    surrogate_real_voltage::AbstractArray
    surrogate_imag_voltage::AbstractArray
    tstops::AbstractArray
    stable::Bool
    solve_time::Float64
end

function SteadyStateNODEData(;
    type = "SteadyStateNODEData",
    tsteps = [],
    ic = Dict{Symbol, Float64}(),
    real_current = [],
    imag_current = [],
    surrogate_real_voltage = [],
    surrogate_imag_voltage = [],
    tstops = [],
    stable = false,
    solve_time = 0.0,
)
    return SteadyStateNODEData(
        type,
        tsteps,
        ic,
        real_current,
        imag_current,
        surrogate_real_voltage,
        surrogate_imag_voltage,
        tstops,
        stable,
        solve_time,
    )
end

function fill_surrogate_data!(
    data::SteadyStateNODEData,
    params::SteadyStateNODEDataParams,
    sys_train::PSY.System,
    psid_perturbations,
    data_collection::GenerateDataParams,
    data_aux::Union{SteadyStateNODEData, Nothing},
    surrogate_params,
)
    tspan = data_collection.tspan
    solver = instantiate_solver(data_collection)
    abstol = data_collection.solver_tols.abstol
    reltol = data_collection.solver_tols.reltol
    #display(PSY.solve_powerflow(sys_train)["bus_results"])
    #display(PSY.solve_powerflow(sys_train)["flow_results"])
    if data_aux !== nothing
        match_operating_point(sys_train, data_aux, surrogate_params)    #TODO - not tested
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
    ic = PSID.read_initial_conditions(sim_full)
    location_data_collection = params.location_of_data_collection
    if location_data_collection[1][2] == :from || location_data_collection[1][2] == :to
        _fill_ic_branch!(data, ic, location_data_collection, sys_train, data_collection)
    elseif location_data_collection[1][2] == :source
        _fill_ic_source!(data, ic, location_data_collection, sys_train)   #TODO - test this function which collects data from a source! 
    else
        @error "Invaled entry for initial condition collection location"
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
        if location_data_collection[1][2] !== :source
            @assert isapprox(data.ic[:Ir0], data.real_current[1]; atol = 1e-14)
            @assert isapprox(data.ic[:Ii0], data.imag_current[1]; atol = 1e-14)
            @assert isapprox(data.ic[:Vr0], data.surrogate_real_voltage[1]; atol = 1e-14)
            @assert isapprox(data.ic[:Vi0], data.surrogate_imag_voltage[1]; atol = 1e-14)
        end
    end
end

"""
Matches the operating point from the ground truth dataset when generating the dataset for a surrogate model. 
"""
function match_operating_point(sys, data_aux, surrogate_params)
    settings_unit_cache = deepcopy(sys.units_settings.unit_system)
    PSY.set_units_base_system!(sys, "SYSTEM_BASE")
    Vr0 = data_aux.ic[:Vr0]
    Vi0 = data_aux.ic[:Vi0]
    Ir0 = data_aux.ic[:Ir0]
    Ii0 = data_aux.ic[:Ii0]
    P0 = Vr0 * Ir0 + Vi0 * Ii0
    Q0 = Vi0 * Ir0 - Vr0 * Ii0
    Vm0 = sqrt(Vr0^2 + Vi0^2)
    θ0 = atan(Vi0, Vr0)
    @info "operating point to match with surrogate model:  ",
    "Vm0: ",
    Vm0,
    "θ0: ",
    θ0,
    "P0: ",
    P0,
    "Q0: ",
    Q0
    _match_operating_point(sys, P0, Q0, Vm0, θ0, surrogate_params)
    PSY.set_units_base_system!(sys, settings_unit_cache)
end

function _match_operating_point(
    sys,
    P0,
    Q0,
    Vm0,
    θ0,
    surrogate_params::Union{SteadyStateNODEObsParams, SteadyStateNODEParams},
)
    for s in PSY.get_components(
        x -> typeof(PSY.get_dynamic_injector(x)) == SteadyStateNODE,
        PSY.Source,
        sys,
    )
        PSY.set_active_power!(s, P0)
        PSY.set_reactive_power!(s, Q0)
        PSY.set_internal_voltage!(s, Vm0)
        PSY.set_internal_angle!(s, θ0)
    end
end

function _match_operating_point(
    sys,
    P0,
    Q0,
    Vm0,
    θ0,
    surrogate_params::Union{ClassicGenParams, GFLParams, GFMParams},
)
    for s in PSY.get_components(
        x -> PSY.get_name(x) == surrogate_params.name,
        PSY.StaticInjection,
        sys,
    )
        base_power = PSY.get_base_power(s)
        base_power_dyanmic = PSY.get_base_power(PSY.get_dynamic_injector(s))
        @assert base_power == base_power_dyanmic
        PSY.set_active_power!(s, P0)
        PSY.set_reactive_power!(s, Q0)
    end
end

function _match_operating_point(sys, P0, Q0, Vm0, θ0, surrogate_params::ZIPParams)
    orientation_scale = -1.0
    load = PSY.get_component(PSY.StandardLoad, sys, surrogate_params.name)
    total_P =
        PSY.get_max_impedance_active_power(load) +
        PSY.get_max_current_active_power(load) +
        PSY.get_max_constant_active_power(load)
    total_Q =
        PSY.get_max_impedance_reactive_power(load) +
        PSY.get_max_current_reactive_power(load) +
        PSY.get_max_constant_reactive_power(load)
    PSY.set_impedance_active_power!(
        load,
        orientation_scale * P0 * PSY.get_max_impedance_active_power(load) / total_P,
    )
    PSY.set_impedance_reactive_power!(
        load,
        orientation_scale * Q0 * PSY.get_max_impedance_reactive_power(load) / total_Q,
    )

    PSY.set_current_active_power!(
        load,
        orientation_scale * P0 * PSY.get_max_current_active_power(load) / total_P,
    )
    PSY.set_current_reactive_power!(
        load,
        orientation_scale * Q0 * PSY.get_max_current_reactive_power(load) / total_Q,
    )

    PSY.set_constant_active_power!(
        load,
        orientation_scale * P0 * PSY.get_max_constant_active_power(load) / total_P,
    )
    PSY.set_constant_reactive_power!(
        load,
        orientation_scale * Q0 * PSY.get_max_constant_reactive_power(load) / total_Q,
    )
end

#TODO - Assumes static devices are loads and dynamic devices are generators 
function _match_operating_point(sys, P0, Q0, Vm0, θ0, surrogate_params::MultiDeviceParams)
    P_device_available = []
    Q_load = []
    for s in surrogate_params.static_devices
        if typeof(s) == ZIPParams
            orientation_scale = -1.0
            device = PSY.get_component(PSY.StandardLoad, sys, s.name)
            active_power_available =
                orientation_scale * (
                    PSY.get_max_impedance_active_power(device) +
                    PSY.get_max_current_active_power(device) +
                    PSY.get_max_constant_active_power(device)
                )
            reactive_power =
                orientation_scale * (
                    PSY.get_impedance_reactive_power(device) +
                    PSY.get_current_reactive_power(device) +
                    PSY.get_constant_reactive_power(device)
                )
            push!(P_device_available, active_power_available)
            push!(Q_load, reactive_power)
        else
            device = PSY.get_component(PSY.Component, sys, s.name)
            push!(P_device_available, PSY.get_max_active_power(device))
        end
    end
    for s in surrogate_params.dynamic_devices
        device = PSY.get_component(PSY.StaticInjection, sys, s.name)
        push!(P_device_available, PSY.get_max_active_power(device))
    end
    P_net = sum(P_device_available)
    n_gens = length(surrogate_params.dynamic_devices)
    P_device = P_device_available ./ P_net .* P0
    Q_remaining = Q0 - sum(Q_load)
    Q_split = Q_remaining / n_gens
    Q_device = vcat(Q_load, ones(n_gens) .* Q_split)
    ix = 1
    for s in surrogate_params.static_devices
        _match_operating_point(sys, P_device[ix], Q_device[ix], Vm0, θ0, s)
        ix += 1
    end
    for s in surrogate_params.dynamic_devices
        _match_operating_point(sys, P_device[ix], Q_device[ix], Vm0, θ0, s)
        ix += 1
    end
end

function _fill_ic_source!(data, ic, connecting_sources, sys_train)
    @error "Methods not yet implemented to collect initial conditions from source"
    return
    Ir0 = zero(Float64)
    Vr0 = zero(Float64)
    Ii0 = zero(Float64)
    Vi0 = zero(Float64)
    for (i, source_tuple) in enumerate(connecting_sources)
        source_name = source_tuple[1]
        source = PSY.get_component(PSY.Source, sys_train, source_name)
        dyn_source = PSY.get_dynamic_injector(source)
        display(ic)
        #Note: currents multiplied by -1 to get current out of the surrogate (into the perturbing source)
        if dyn_source === nothing
            #Ir0 = -1 .* ic[source_name][:Ir]
            #Ii0 = -1 .* ic[source_name][:Ii]
        else
            #Ir0 = -1 .* ic[source_name][:Ir]
            #Ii0 = -1 .* ic[source_name][:Ii]
        end
        bus_number_source = PSY.get_number(PSY.get_bus(source))
        Vr0 = ic["V_R"][bus_number_source]
        Vi0 = ic["V_I"][bus_number_source]
    end
    data.ic =
        Dict{Symbol, AbstractArray}(:Vr0 => Vr0, :Vi0 => Vi0, :Ir0 => Ir0, :Ii0 => Ii0)
end

#NOTE: this function will fail if you have double lines with reversed to/from orientation.
#NOTE: need to fix current for tripped line.
function _fill_ic_branch!(data, ic, connecting_branch_data, sys_train, data_collection)
    connecting_branches = [
        PSY.get_component(PSY.Branch, sys_train, branch_tuple[1]) for
        branch_tuple in connecting_branch_data
    ]
    connecting_arcs = unique([PSY.get_arc(branch) for branch in connecting_branches])
    Ir0 = zero(Float64)
    Vr0 = zero(Float64)
    Ii0 = zero(Float64)
    Vi0 = zero(Float64)
    Ir0_from_to = zero(Float64)
    Ii0_from_to = zero(Float64)
    for (i, arc) in enumerate(connecting_arcs)
        corresponding_branches =
            collect(PSY.get_components(x -> PSY.get_arc(x) == arc, PSY.Branch, sys_train))
        for b in corresponding_branches
            branch_name = PSY.get_name(b)
            #TODO - get rid of the if/else logic below after this PSID issue is resolve: https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/issues/283
            if data_collection.all_lines_dynamic
                Ir0_from_to += ic[string("Line ", branch_name)][:Il_R]
                Ii0_from_to += ic[string("Line ", branch_name)][:Il_I]
            else
                bus_number_from = PSY.get_number(PSY.get_from(PSY.get_arc(b)))
                bus_number_to = PSY.get_number(PSY.get_to(PSY.get_arc(b)))
                Vr0_from = ic["V_R"][bus_number_from]
                Vi0_from = ic["V_I"][bus_number_from]
                Vr0_to = ic["V_R"][bus_number_to]
                Vi0_to = ic["V_I"][bus_number_to]
                r = PSY.get_r(b)
                x = PSY.get_x(b)
                I_flow =
                    ((Vr0_from + Vi0_from * 1im) - (Vr0_to + Vi0_to * 1im)) ./ (r + x * 1im)
                Ir0_from_to += real(I_flow)
                Ii0_from_to += imag(I_flow)
            end
        end
        branch_name = PSY.get_name(corresponding_branches[1])   #assume data same for all corresponding_branches
        surrogate_location = connecting_branch_data[1][2]
        branch = PSY.get_component(PSY.ACBranch, sys_train, branch_name)
        if surrogate_location == :from
            Ir0 = Ir0_from_to
            Ii0 = Ii0_from_to
            bus_number_surrogate = PSY.get_number(PSY.get_from(PSY.get_arc(branch)))
            bus_number_opposite = PSY.get_number(PSY.get_to(PSY.get_arc(branch)))
            Vr0 = ic["V_R"][bus_number_surrogate]
            Vi0 = ic["V_I"][bus_number_surrogate]
        elseif surrogate_location == :to
            Ir0 = Ir0_from_to * -1
            Ii0 = Ii0_from_to * -1
            bus_number_surrogate = PSY.get_number(PSY.get_to(PSY.get_arc(branch)))
            bus_number_opposite = PSY.get_number(PSY.get_from(PSY.get_arc(branch)))
            Vr0 = ic["V_R"][bus_number_surrogate]
            Vi0 = ic["V_I"][bus_number_surrogate]
        end
    end
    data.ic = Dict{Symbol, Float64}(:Vr0 => Vr0, :Vi0 => Vi0, :Ir0 => Ir0, :Ii0 => Ii0)
end

function _fill_data_source!(data, results, connecting_sources, save_indices, sys_train)
    n_save_points = length(save_indices)
    real_current = zeros(length(connecting_sources), n_save_points)
    imag_current = zeros(length(connecting_sources), n_save_points)
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
    end
    data.tstops = unique(results.solution.t)
    data.tsteps = unique(results.solution.t)[save_indices]
    data.real_current = real_current
    data.imag_current = imag_current
    data.surrogate_real_voltage = surrogate_real_voltage
    data.surrogate_imag_voltage = surrogate_imag_voltage
    data.stable = true
    data.solve_time = results.time_log[:timed_solve_time]
end

function _fill_data_states(data, results, dynamic_device_name, save_indices, sys_train)
    state_trajectories_dict = Dict{Symbol, AbstractArray}()
    dynamic_injector =
        PSY.get_component(PSY.DynamicInjection, sys_train, dynamic_device_name)
    for s in PSY.get_states(dynamic_injector)
        state_trajectories_dict[s] =
            PSID.get_state_series(results, (dynamic_device_name, s))[2][save_indices]
    end
    data.states = state_trajectories_dict
end

#NOTE: this function will fail if you have double lines with reversed to/from orientation.
#NOTE: need to fix current for tripped line.
function _fill_data_branch!(
    data,
    results,
    connecting_branch_data,
    save_indices,
    sys_train,
    data_collection,
)
    connecting_branches = [
        PSY.get_component(PSY.Branch, sys_train, branch_tuple[1]) for
        branch_tuple in connecting_branch_data
    ]
    connecting_arcs = unique([PSY.get_arc(branch) for branch in connecting_branches])
    n_save_points = length(save_indices)
    real_current = zeros(length(connecting_arcs), n_save_points)
    imag_current = zeros(length(connecting_arcs), n_save_points)
    surrogate_real_voltage = zeros(length(connecting_arcs), n_save_points)
    surrogate_imag_voltage = zeros(length(connecting_arcs), n_save_points)
    for (i, arc) in enumerate(connecting_arcs)
        corresponding_branches =
            collect(PSY.get_components(x -> PSY.get_arc(x) == arc, PSY.Branch, sys_train))
        Ir_from_to = zeros(n_save_points)
        Ii_from_to = zeros(n_save_points)
        for b in corresponding_branches
            branch_name = PSY.get_name(b)
            #TODO - get rid of the if/else logic below after this PSID issue is resolve: https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/issues/283
            if data_collection.all_lines_dynamic
                Ir_from_to +=
                    PSID.get_state_series(results, (branch_name, :Il_R))[2][save_indices]
                Ii_from_to +=
                    PSID.get_state_series(results, (branch_name, :Il_I))[2][save_indices]
            else
                Ir_from_to +=
                    PSID.get_real_current_branch_flow(results, branch_name)[2][save_indices]
                Ii_from_to +=
                    PSID.get_imaginary_current_branch_flow(results, branch_name)[2][save_indices]
            end
        end

        branch_name = PSY.get_name(corresponding_branches[1])   #assume data same for all corresponding_branches
        surrogate_location = connecting_branch_data[1][2]
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
    end
    data.tstops = unique(results.solution.t)
    data.tsteps = unique(results.solution.t)[save_indices]
    data.real_current = real_current
    data.imag_current = imag_current
    data.surrogate_real_voltage = surrogate_real_voltage
    data.surrogate_imag_voltage = surrogate_imag_voltage
    data.stable = true
    data.solve_time = results.time_log[:timed_solve_time]
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

#################################################################
#make a different type of dataset which captures all of the states of a device, along with the input and output voltage as above. 
#just add -> states::Dict{Symbol, AbstractArray}, then loop through all of the states and get the values for the device. 
mutable struct AllStatesDataParams <: SurrogateDatasetParams
    type::String
    dynamic_device_name::String
    location_of_data_collection::Array{Tuple{String, Symbol}}
end

function AllStatesDataParams(;
    type = "AllStatesDataParams",
    dynamic_device_name = "",
    location_of_data_collection = [],           #TODO - This is [("branchname", :to)] [("branchname", :from)]  or  --> extend to allow for [("sourcename", :source)]   
)
    return AllStatesDataParams(type, dynamic_device_name, location_of_data_collection)
end

mutable struct AllStatesData <: SurrogateDataset
    type::String
    tsteps::AbstractArray
    real_current::AbstractArray
    imag_current::AbstractArray
    surrogate_real_voltage::AbstractArray
    surrogate_imag_voltage::AbstractArray
    states::Dict{Symbol, AbstractArray}
    tstops::AbstractArray
    stable::Bool
    solve_time::Float64
end

function AllStatesData(;
    type = "AllStatesData",
    tsteps = [],
    real_current = [],
    imag_current = [],
    surrogate_real_voltage = [],
    surrogate_imag_voltage = [],
    states = Dict{Symbol, AbstractArray}(),
    tstops = [],
    stable = false,
    solve_time = 0.0,
)
    return AllStatesData(
        type,
        tsteps,
        real_current,
        imag_current,
        surrogate_real_voltage,
        surrogate_imag_voltage,
        states,
        tstops,
        stable,
        solve_time,
    )
end

function fill_surrogate_data!(
    data::AllStatesData,
    params::AllStatesDataParams,
    sys_train::PSY.System,
    psid_perturbations,
    data_collection::GenerateDataParams,
    data_aux::Union{AllStatesData, Nothing},
    surrogate_params,
)
    tspan = data_collection.tspan
    solver = instantiate_solver(data_collection)
    abstol = data_collection.solver_tols.abstol
    reltol = data_collection.solver_tols.reltol
    #display(PSY.solve_powerflow(sys_train)["bus_results"])
    #display(PSY.solve_powerflow(sys_train)["flow_results"])
    if data_aux !== nothing
        match_operating_point(sys_train, data_aux, surrogate_params)    #TODO - not tested
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
            #fill the states dictionary with data. 
            _fill_data_states(
                data,
                results,
                params.dynamic_device_name,
                save_indices,
                sys_train,
            )
        elseif location_data_collection[1][2] == :source
            _fill_data_source!(
                data,
                results,
                location_data_collection,
                save_indices,
                sys_train,
            )   #TODO - test this function which collects data from a source! 
            _fill_data_states(
                data,
                results,
                params.dynamic_device_name,
                save_indices,
                sys_train,
            )
        else
            @error "Invaled entry for data collection location"
        end
    end
end

function EmptyTrainDataSet(key::AllStatesDataParams)
    return AllStatesData()
end
