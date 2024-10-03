mutable struct TerminalData <: SurrogateDataset
    type::String
    tsteps::AbstractArray
    tstops::AbstractArray
    stable::Bool
    built::Bool
    solve_time::Float64
    device_terminal_data::Dict{String, Dict{Symbol, AbstractArray}} #string is device name, symbol is :vr, :vi, :ir, :ii, :p, :q, arrays are the data
end

function TerminalData(;
    type = "TerminalData",
    tsteps = [],
    tstops = [],
    stable = false,
    built = false,
    solve_time = 0.0,
    device_terminal_data = Dict{String, Dict{Symbol, AbstractArray}}(),
)
    return TerminalData(
        type,
        tsteps,
        tstops,
        stable,
        built,
        solve_time,
        device_terminal_data,
    )
end

function fill_surrogate_data!(
    data::TerminalData,
    device_details,
    data_collection_params,
    sim_full,
)
    terminal_data_dict = Dict{String, Dict{Symbol, AbstractArray}}()
    sys = sim_full.sys
    for (device_name, orientation_details) in device_details
        all_components_with_name =
            PSY.get_components_by_name(PSY.Component, sys, device_name)
        exclude_dynamic_injectors =
            filter!(x -> !(typeof(x) <: PSY.DynamicInjection), all_components_with_name)
        @assert length(exclude_dynamic_injectors) == 1
        device = exclude_dynamic_injectors[1]
        _fill_terminal_data!(
            terminal_data_dict,
            sim_full,
            device,
            orientation_details,
            data_collection_params,
            sys,
        )
    end
    if sim_full.status == PSID.SIMULATION_FINALIZED
        results = PSID.read_results(sim_full)
        if length(data_collection_params.tsave) == 0
            save_indices = 1:length(unique(results.solution.t))
        else
            save_indices = indexin(data_collection_params.tsave, unique(results.solution.t))
        end
        data.tstops = unique(results.solution.t)
        data.tsteps = unique(results.solution.t)[save_indices]
        data.stable = true
        data.built = true
        data.solve_time = results.time_log[:timed_solve_time]
        data.device_terminal_data = terminal_data_dict
    elseif sim_full.status == PSID.SIMULATION_FAILED
        data.stable = false
        data.built = true
    elseif sim_full.status == PSID.BUILD_FAILED
        data.stable = false
        data.built = false
    else
        @error "Cannot hanldle the PSID simulation status $(sim_full.status)"
    end
end

function _fill_terminal_data!(
    terminal_data_dict,
    results,
    device,
    orientation_details,
    data_collection_params,
    sys,
)
    @error "Function _fill_terminal_data! not defined for type $(typeof(device)) "
    return
end

function _fill_terminal_data!(
    terminal_data_dict,
    results,
    device::Union{PSY.Line, PSY.Transformer2W, PSY.DynamicBranch},
    orientation_details,
    data_collection_params,
    sys,
)
    @error "Line or transformer name was passed, running _fill_terminal_data for corresponding ARC"
    arc = PSY.get_arc(device)
    return _fill_terminal_data!(
        terminal_data_dict,
        results,
        arc::PSY.Arc,
        orientation_details,
        data_collection_params,
        sys,
    )
end

function _fill_terminal_data!(
    terminal_data_dict,
    sim_full,
    device::S,
    orientation_details,
    data_collection_params,
    sys,
) where {S <: PSY.StaticInjection}
    data_dict = Dict{Symbol, AbstractArray}()
    if sim_full.status == PSID.SIMULATION_FINALIZED
        results = PSID.read_results(sim_full)
        if length(data_collection_params.tsave) == 0
            save_indices = 1:length(unique(results.solution.t))
        else
            save_indices = indexin(data_collection_params.tsave, unique(results.solution.t))
        end
        device_name = PSY.get_name(device)
        dyn_injector = PSY.get_dynamic_injector(device)
        if orientation_details[:direction] == :in
            scale = -1.0
        elseif orientation_details[:direction] == :out
            scale = 1.0
        else
            @error "invalid value of direction: must be :in or :out"
            return
        end
        if dyn_injector === nothing
            @assert typeof(device) == PSY.Source
            data_dict[:ir] =
                scale .*
                PSID.get_source_real_current_series(results, device_name)[2][save_indices]
            data_dict[:ii] =
                scale .*
                PSID.get_source_imaginary_current_series(results, device_name)[2][save_indices]
        else
            data_dict[:ir] =
                scale .* PSID.get_real_current_series(results, device_name)[2][save_indices]
            data_dict[:ii] =
                scale .*
                PSID.get_imaginary_current_series(results, device_name)[2][save_indices]
        end
        data_dict[:p] =
            scale .* PSID.get_activepower_series(results, device_name)[2][save_indices]
        data_dict[:q] =
            scale .* PSID.get_reactivepower_series(results, device_name)[2][save_indices]
        bus_number_source = PSY.get_number(PSY.get_bus(device))
        V_surrogate =
            PSID.get_voltage_magnitude_series(results, bus_number_source)[2][save_indices]
        θ_surrogate =
            PSID.get_voltage_angle_series(results, bus_number_source)[2][save_indices]
        data_dict[:vr] = V_surrogate .* cos.(θ_surrogate)
        data_dict[:vi] = V_surrogate .* sin.(θ_surrogate)
    else
        @warn "System unstable, not recording any TerminalData for StaticInjection"
    end

    terminal_data_dict[device_name] = data_dict
end

function _fill_terminal_data!(
    terminal_data_dict,
    sim_full,
    device::PSY.Arc,
    orientation_details,
    data_collection_params,
    sys,
)
    data_dict = Dict{Symbol, AbstractArray}()
    if sim_full.status == PSID.SIMULATION_FINALIZED
        results = PSID.read_results(sim_full)
        if length(data_collection_params.tsave) == 0
            save_indices = 1:length(unique(results.solution.t))
        else
            save_indices = indexin(data_collection_params.tsave, unique(results.solution.t))
        end
        corresponding_branches =
            collect(PSY.get_components(x -> PSY.get_arc(x) == device, PSY.Branch, sys))
        device_name = PSY.get_name(device)
        Ir_from_to = zeros(length(save_indices))
        Ii_from_to = zeros(length(save_indices))
        for b in corresponding_branches
            branch_name = PSY.get_name(b)
            #TODO - get rid of the if/else logic below after this PSID issue is resolve: https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/issues/283
            if data_collection_params.all_lines_dynamic && typeof(b) == PSY.Line ||
               data_collection_params.all_branches_dynamic
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
        branch = corresponding_branches[1] #can get voltage info based on a single branch.
        if orientation_details[:direction] == :in
            scale = 1.0
        elseif orientation_details[:direction] == :out
            scale = -1.0
        else
            @error "invalid value of direction: must be :in or :out"
            @assert false
        end
        if orientation_details[:side] == :from
            data_dict[:ir] = Ir_from_to * scale
            data_dict[:ii] = Ii_from_to * scale
            bus_number_voltage_measurement =
                PSY.get_number(PSY.get_from(PSY.get_arc(branch)))
        elseif orientation_details[:side] == :to
            data_dict[:ir] = Ir_from_to * scale * -1.0
            data_dict[:ii] = Ii_from_to * scale * -1.0
            bus_number_voltage_measurement = PSY.get_number(PSY.get_to(PSY.get_arc(branch)))
        else
            @error "invalid value of side: must be :to or :from"
            @assert false
        end
        V =
            PSID.get_voltage_magnitude_series(results, bus_number_voltage_measurement)[2][save_indices]
        θ =
            PSID.get_voltage_angle_series(results, bus_number_voltage_measurement)[2][save_indices]
        data_dict[:vr] = V .* cos.(θ)
        data_dict[:vi] = V .* sin.(θ)
        data_dict[:p] = data_dict[:vr] .* data_dict[:ir] .- data_dict[:vi] .* data_dict[:ii]
        data_dict[:q] = data_dict[:vr] .* data_dict[:ii] .+ data_dict[:vi] .* data_dict[:ir]
    else
        ic = PSID.read_initial_conditions(sim_full)
        corresponding_branches =
            collect(PSY.get_components(x -> PSY.get_arc(x) == device, PSY.Branch, sys))
        device_name = PSY.get_name(device)
        Ir0_from_to = 0.0
        Ii0_from_to = 0.0
        for b in corresponding_branches
            branch_name = PSY.get_name(b)
            #TODO - get rid of the if/else logic below after this PSID issue is resolve: https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/issues/283
            if data_collection_params.all_lines_dynamic
                Ir0_from_to += ic[string("Line ", branch_name)][:Il_R]  # PSID.get_state_series(results, (branch_name, :Il_R))[2][save_indices]
                Ii0_from_to += ic[string("Line ", branch_name)][:Il_I] # PSID.get_state_series(results, (branch_name, :Il_I))[2][save_indices]
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
        branch = corresponding_branches[1] #can get voltage info based on a single branch.
        if orientation_details[:direction] == :in
            scale = 1.0
        elseif orientation_details[:direction] == :out
            scale = -1.0
        else
            @error "invalid value of direction: must be :in or :out"
            @assert false
        end
        if orientation_details[:side] == :from
            data_dict[:ir] = [Ir0_from_to * scale]
            data_dict[:ii] = [Ii0_from_to * scale]
            bus_number_voltage_measurement =
                PSY.get_number(PSY.get_from(PSY.get_arc(branch)))
        elseif orientation_details[:side] == :to
            data_dict[:ir] = [Ir0_from_to * scale * -1.0]
            data_dict[:ii] = [Ii0_from_to * scale * -1.0]
            bus_number_voltage_measurement = PSY.get_number(PSY.get_to(PSY.get_arc(branch)))
        else
            @error "invalid value of side: must be :to or :from"
            @assert false
        end
        Vr0 = ic["V_R"][bus_number_voltage_measurement]
        Vi0 = ic["V_I"][bus_number_voltage_measurement]
        data_dict[:vr] = [Vr0]
        data_dict[:vi] = [Vi0]
        data_dict[:p] = data_dict[:vr] .* data_dict[:ir] .- data_dict[:vi] .* data_dict[:ii]
        data_dict[:q] = data_dict[:vr] .* data_dict[:ii] .+ data_dict[:vi] .* data_dict[:ir]
    end
    terminal_data_dict[device_name] = data_dict
end

"""
Matches the operating point from the ground truth dataset when generating the dataset for a surrogate model. 
"""
function match_operating_point(sys, data_aux::TerminalData, surrogate_params)
    settings_unit_cache = deepcopy(sys.units_settings.unit_system)
    PSY.set_units_base_system!(sys, "SYSTEM_BASE")
    @assert data_aux.built == true
    for (_, v) in data_aux.device_terminal_data
        Vr0 = v[:vr][1]
        Vi0 = v[:vi][1]
        Ir0 = v[:ir][1]
        Ii0 = v[:ii][1]
        P0 = Vr0 * Ir0 + Vi0 * Ii0
        Q0 = Vi0 * Ir0 - Vr0 * Ii0
        Vm0 = sqrt(Vr0^2 + Vi0^2)
        θ0 = atan(Vi0, Vr0)
        _match_operating_point(sys, P0, Q0, Vm0, θ0, surrogate_params)
    end
    PSY.set_units_base_system!(sys, settings_unit_cache)
end

function _match_operating_point(sys, P0, Q0, Vm0, θ0, surrogate_params::SourceParams)
    for s in PSY.get_components(PSY.Source, sys)
        PSY.set_active_power!(s, P0)
        PSY.set_reactive_power!(s, Q0)
        PSY.set_internal_voltage!(s, Vm0)
        PSY.set_internal_angle!(s, θ0)
        b = PSY.get_bus(s)
        PSY.set_magnitude!(b, Vm0)
        PSY.set_angle!(b, θ0)
    end
end

function _match_operating_point(sys, P0, Q0, Vm0, θ0, surrogate_params::SourceLoadParams)
    for s in PSY.get_components(SourceLoad, sys)
        PSY.set_active_power!(s, P0)
        PSY.set_reactive_power!(s, Q0)
        #PSY.set_internal_voltage!(s, Vm0)
        #PSY.set_internal_angle!(s, θ0)
        b = PSY.get_bus(s)
        PSY.set_magnitude!(b, Vm0)
        PSY.set_angle!(b, θ0)
    end
end

function _match_operating_point(
    sys,
    P0,
    Q0,
    Vm0,
    θ0,
    surrogate_params::SteadyStateNODEParams,
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

#Note: Don't set reactive power in _match_operating_point for surrogates that can have multiple components
#The redistribution of reactive power is handled by the post-processing in PowerFlows.jl to match the surrogate implementation.
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
        PSY.set_active_power!(s, P0)
    end
end

#Note: Don't set reactive power in _match_operating_point for surrogates that can have multiple components
#The redistribution of reactive power is handled by the post-processing in PowerFlows.jl to match the surrogate implementation.
function _match_operating_point(sys, P0, Q0, Vm0, θ0, surrogate_params::ZIPParams)
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
        -1 * P0 * PSY.get_max_impedance_active_power(load) / total_P,
    )
    PSY.set_current_active_power!(
        load,
        -1 * P0 * PSY.get_max_current_active_power(load) / total_P,
    )
    PSY.set_constant_active_power!(
        load,
        -1 * P0 * PSY.get_max_constant_active_power(load) / total_P,
    )
end

#Note: Don't set reactive power in _match_operating_point for surrogates that can have multiple components
#The redistribution of reactive power is handled by the post-processing in PowerFlows.jl to match the surrogate implementation.
function _match_operating_point(sys, P0, Q0, Vm0, θ0, surrogate_params::MultiDeviceParams)
    P_device_available = []
    load_polarity_factor = -1.0
    for s in surrogate_params.static_devices
        if typeof(s) == ZIPParams
            device = PSY.get_component(PSY.StandardLoad, sys, s.name)
            active_power_available =
                PSY.get_max_impedance_active_power(device) +
                PSY.get_max_current_active_power(device) +
                PSY.get_max_constant_active_power(device)
            push!(P_device_available, load_polarity_factor * active_power_available)
        else
            device = PSY.get_component(PSY.Component, sys, s.name)
            push!(
                P_device_available,
                load_polarity_factor * PSY.get_max_active_power(device),
            )
        end
    end
    for s in surrogate_params.dynamic_devices
        device = PSY.get_component(PSY.StaticInjection, sys, s.name)
        push!(P_device_available, PSY.get_max_active_power(device))
    end
    P_total_available = sum(P_device_available)
    P_device = P_device_available ./ P_total_available .* P0
    ix = 1
    for s in surrogate_params.static_devices
        _match_operating_point(sys, P_device[ix], nothing, Vm0, θ0, s)
        ix += 1
    end
    for s in surrogate_params.dynamic_devices
        _match_operating_point(sys, P_device[ix], nothing, Vm0, θ0, s)
        ix += 1
    end
end

function instantiate_solver(inputs)
    if occursin("MethodOfSteps", inputs.solver)
        return solver_map(inputs.solver)
    else
        return solver_map(inputs.solver)
    end
end

function solver_map(key)
    d = Dict(
        "Rodas4" => OrdinaryDiffEq.Rodas4(),
        "Rodas5" => OrdinaryDiffEq.Rodas5(),
        "Rodas5P" => OrdinaryDiffEq.Rodas5P(),
        "Rodas4(autodiff=false)" => OrdinaryDiffEq.Rodas4(autodiff = false),
        "Rodas5(autodiff=false)" => OrdinaryDiffEq.Rodas5(autodiff = false),
        "Rodas5P(autodiff=false)" => OrdinaryDiffEq.Rodas5P(autodiff = false),
        "TRBDF2" => OrdinaryDiffEq.TRBDF2(),
        "Tsit5" => OrdinaryDiffEq.Tsit5(),
        "IDA" => Sundials.IDA(),
        "MethodOfSteps(Rodas5(autodiff=false))" =>
            DelayDiffEq.MethodOfSteps(OrdinaryDiffEq.Rodas5(autodiff = false)),
        "MethodOfSteps(Rodas5P(autodiff=false))" =>
            DelayDiffEq.MethodOfSteps(OrdinaryDiffEq.Rodas5P(autodiff = false)),
        "MethodOfSteps(Rodas5)" =>
            DelayDiffEq.MethodOfSteps(OrdinaryDiffEq.Rodas5()),
        "MethodOfSteps(Rodas5P)" =>
            DelayDiffEq.MethodOfSteps(OrdinaryDiffEq.Rodas5P()),
        #"IDA(linear_solver = :KLU)" => Sundials.IDA(linear_solver = :KLU),
    )
    return d[key]
end

function EmptyTrainDataSet(T::Type{TerminalData})
    return TerminalData()
end

function generate_empty_plot(T::Type{TerminalData})
    p = PlotlyJS.make_subplots(
        rows = 2,
        cols = 2,
        specs = [
            PlotlyJS.Spec() PlotlyJS.Spec()
            PlotlyJS.Spec() PlotlyJS.Spec()
        ],
        subplot_titles = ["vr" "vi" "ir" "ii"],
        vertical_spacing = 0.1,
    )
    return p
end

function add_data_trace!(p, data::TerminalData; mode = "lines", color = "", name = "")
    for (device_name, device_data_dict) in data.device_terminal_data
        Vr = device_data_dict[:vr]
        Vi = device_data_dict[:vi]
        Ir = device_data_dict[:ir]
        Ii = device_data_dict[:ii]
        #P = device_data_dict[:p]
        #Q = device_data_dict[:q]

        PlotlyJS.add_trace!(
            p,
            PlotlyJS.scatter(;
                x = data.tsteps,
                y = Vr,
                marker = PlotlyJS.attr(color = color),
                mode = mode,
                name = string(device_name, name),
            ),
            row = 1,
            col = 1,
        )
        PlotlyJS.add_trace!(
            p,
            PlotlyJS.scatter(;
                x = data.tsteps,
                y = Vi,
                marker = PlotlyJS.attr(color = color),
                mode = mode,
                #showlegend=false,
                name = string(device_name, name),
            ),
            row = 1,
            col = 2,
        )
        PlotlyJS.add_trace!(
            p,
            PlotlyJS.scatter(;
                x = data.tsteps,
                y = Ir,
                marker = PlotlyJS.attr(color = color),
                mode = mode,
                #showlegend=false,
                name = string(device_name, name),
            ),
            row = 2,
            col = 1,
        )
        PlotlyJS.add_trace!(
            p,
            PlotlyJS.scatter(;
                x = data.tsteps,
                y = Ii,
                marker = PlotlyJS.attr(color = color),
                mode = mode,
                #showlegend=false,
                name = string(device_name, name),
            ),
            row = 2,
            col = 2,
        )
    end
    return p
end
