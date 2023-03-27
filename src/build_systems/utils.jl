const SOURCE_X_TH = 1e-6
const SOURCE_R_TH = 1e-6   #If using a source to implement a load step for training, important to have non-zero resistance if lines are dynamic and there is capacitance at the bus. 
"""
    subsystem, location_of_data_collection = create_subsystem_from_buses(sys_full::PSY.System, subsystem_bus_numbers::Vector{Int64}) 

Takes a deepcopy of PSY system and removes all components that are not attached to the indicated buses. 
Leaves all branches that connect two buses within the surrogate AND branches that connect to the rest of the system. 
Adds source components at connecting buses --- these are used to perturb the surrogate when generating training data. 
# Inputs 
- `sys_full::PowerSystems.System`: The full starting system. 
- `surrogate_bus_numbers::Vector{Int64}`: Numbers of the buses to be included in the subsystem. These bus numbers should form a continuous (non-islanded) subsystem.
# Outputs
- `subsystem::PowerSystems.System`:
- `location_of_data_collection::Vector{Tuple{String, Symbol}}`: Tuple of branch name and either `:to` or `:from` for interpreting the polarity of data from those branches.
"""

function create_validation_system_from_buses(
    sys_full::PSY.System,
    surrogate_bus_numbers::Vector{Int64},
)
    nonsurrogate_bus_numbers =
        PSY.get_number.(
            PSY.get_components(
                x -> PSY.get_number(x) ∉ surrogate_bus_numbers,
                PSY.Bus,
                sys_full,
            ),
        )
    powerflow_df = PowerFlows.run_powerflow(sys_full)["flow_results"]
    sys = deepcopy(sys_full)
    (_check_connectivity(sys, nonsurrogate_bus_numbers) == false) && return false
    _remove_static_and_dynamic_injectors!(sys, surrogate_bus_numbers)
    connecting_branch_data =
        _get_connecting_branch_power(sys, surrogate_bus_numbers, powerflow_df) #inner terminal is the bus that is within the numbers provided. For the validation system that is the bus where the source is added. 
    _remove_internal_branches!(sys, surrogate_bus_numbers)

    #Keep the buses within the surrogate that are connecting to the rest of the system (this is where the surrogate component will attach)
    connecting_bus_numbers = Int64[]
    for b in connecting_branch_data
        if b.inner_terminal == :to
            branch = PSY.get_component(PSY.Branch, sys, b.connecting_branch_name)
            push!(connecting_bus_numbers, PSY.get_number(PSY.get_to(PSY.get_arc(branch))))
        elseif b.inner_terminal == :from
            branch = PSY.get_component(PSY.Branch, sys, b.connecting_branch_name)
            push!(connecting_bus_numbers, PSY.get_number(PSY.get_from(PSY.get_arc(branch))))
        end
    end
    #Remove buses that are both in the surrogate and not the connecting bus within the surrogate.
    _remove_buses!(sys, setdiff(surrogate_bus_numbers, connecting_bus_numbers))
    _add_sources!(sys, connecting_branch_data, :inner)
    _ensure_a_reference_bus!(sys)
    location_of_data_collection =
        [(x.connecting_branch_name, x.inner_terminal) for x in connecting_branch_data]
    return sys, location_of_data_collection
end

function create_train_system_from_buses(
    sys_full::PSY.System,
    surrogate_bus_numbers::Vector{Int64},
)
    nonsurrogate_bus_numbers =
        PSY.get_number.(
            PSY.get_components(
                x -> PSY.get_number(x) ∉ surrogate_bus_numbers,
                PSY.Bus,
                sys_full,
            ),
        )
    powerflow_df = PowerFlows.run_powerflow(sys_full)["flow_results"]
    sys = deepcopy(sys_full)
    (_check_connectivity(sys, surrogate_bus_numbers) == false) && return false
    _remove_static_and_dynamic_injectors!(sys, nonsurrogate_bus_numbers)
    connecting_branch_data =
        _get_connecting_branch_power(sys, nonsurrogate_bus_numbers, powerflow_df)       #inner terminal is the bus that is within the numbers provided. For the train system that is not the bus where the source is added. 
    _remove_internal_branches!(sys, nonsurrogate_bus_numbers)
    _remove_buses!(sys, nonsurrogate_bus_numbers)
    connecting_branches = [
        PSY.get_component(PSY.Branch, sys, b.connecting_branch_name) for
        b in connecting_branch_data
    ]
    connecting_arcs = unique([PSY.get_arc(branch) for branch in connecting_branches])
    for a in connecting_arcs
        PSY.remove_component!(sys, a)
    end
    for b in connecting_branches
        PSY.remove_component!(sys, b)
    end
    _add_sources!(sys, connecting_branch_data, :outer)
    _ensure_a_reference_bus!(sys)
    location_of_data_collection =
        [("source_$ix", :source) for (ix, x) in enumerate(connecting_branch_data)]
    return sys, location_of_data_collection
end

function _add_sources!(sys, connecting_branch_data, location)
    #Add a source to the connecting buses within the surrogate with appropriate active power set point. Ensure that the type of bus is either connecting_branch_data or PV
    for (ix, data) in enumerate(connecting_branch_data)
        if (data.inner_terminal == :from) && (location == :inner)   #Building the validation system
            bus = PSY.get_component(PSY.Bus, sys, data.from_bus_name)
            P_source = data.P_from_to
        elseif (data.inner_terminal == :to) && (location == :outer) #Building the train system
            bus = PSY.get_component(PSY.Bus, sys, data.from_bus_name)
            P_source = data.P_from_to * -1
        elseif (data.inner_terminal == :to) && (location == :inner)  #Building the validation system
            bus = PSY.get_component(PSY.Bus, sys, data.to_bus_name)
            P_source = data.P_to_from
        elseif (data.inner_terminal == :from) && (location == :outer) #Building the train system
            bus = PSY.get_component(PSY.Bus, sys, data.to_bus_name)
            P_source = data.P_to_from * -1
        else
            @error "Invalid value for inner_terminal: $(data.inner_terminal)"
        end
        @assert bus !== nothing
        PSY.get_bustype(bus) !== PSY.BusTypes.REF && PSY.set_bustype!(bus, PSY.BusTypes.PV)
        source = PSY.Source(
            name = string("source_", ix),
            active_power = P_source / 100,
            available = true,
            reactive_power = 0.0,
            bus = bus,
            R_th = SOURCE_R_TH,
            X_th = SOURCE_X_TH,
            internal_voltage = 0.0,
            internal_angle = 0.0,
        )
        PSY.add_component!(sys, source)
    end
end
function _ensure_a_reference_bus!(sys)
    reference_buses =
        PSY.get_components(x -> (PSY.get_bustype(x) == PSY.BusTypes.REF), PSY.Bus, sys)
    if length(reference_buses) == 0
        first_bus = collect(PSY.get_components(PSY.Bus, sys))[1]
        PSY.set_bustype!(first_bus, PSY.BusTypes.REF)
        @assert length(
            PSY.get_components(x -> (PSY.get_bustype(x) == PSY.BusTypes.REF), PSY.Bus, sys),
        ) == 1
    end
end

function _check_connectivity(sys, bus_numbers)
    for b in bus_numbers
        bus = collect(PSY.get_components(x -> PSY.get_number(x) == b, PSY.Bus, sys))[1]
        connected_branches = collect(
            PSY.get_components(
                x ->
                    typeof(x) <: PSY.Branch && ((
                        PSY.get_from(PSY.get_arc(x)) == bus ||
                        PSY.get_to(PSY.get_arc(x)) == bus
                    )),
                PSY.Component,
                sys,
            ),
        )
        x = filter(
            x ->
                PSY.get_number(PSY.get_from(PSY.get_arc(x))) ∈ bus_numbers &&
                    PSY.get_number(PSY.get_to(PSY.get_arc(x))) ∈ bus_numbers,
            connected_branches,
        )
        if length(bus_numbers) > 1 && length(x) == 0
            @error "The specified surrogate bus numbers do not give a validation system thas is connected"
            return false
        end
    end
    return true
end

function _remove_static_and_dynamic_injectors!(sys, bus_numbers)
    static_injectors = collect(
        PSY.get_components(x -> typeof(x) <: PSY.StaticInjection, PSY.Component, sys),
    )
    for static_injector in static_injectors
        if (PSY.get_number(PSY.get_bus(static_injector)) ∈ bus_numbers)
            dynamic_injector = PSY.get_dynamic_injector(static_injector)
            (dynamic_injector !== nothing) && PSY.remove_component!(sys, dynamic_injector)
            PSY.remove_component!(sys, static_injector)
        end
    end
end

function _remove_buses!(sys, bus_numbers)
    buses = collect(PSY.get_components(x -> PSY.get_number(x) ∈ bus_numbers, PSY.Bus, sys))
    for b in buses
        PSY.remove_component!(sys, b)
    end
end

#Removes only the branches within the bus_numbers (both sides)
function _remove_internal_branches!(sys, bus_numbers)
    internal_branches = collect(
        PSY.get_components(
            x ->
                (PSY.get_number(PSY.get_from(PSY.get_arc(x)))) ∈ bus_numbers &&
                    (PSY.get_number(PSY.get_to(PSY.get_arc(x)))) ∈ bus_numbers,
            PSY.Branch,
            sys,
        ),
    )
    for branch in internal_branches
        arc = PSY.get_arc(branch)
        PSY.remove_component!(sys, arc)
        PSY.remove_component!(sys, branch)
    end
end

function _get_connecting_branch_power(sys, bus_numbers, powerflow_df)
    #=     connecting_branch_data = NamedTuple{
            (:connecting_branch_name, :to_bus_name, :from_bus_name, :inner_terminal, :P_from_to, :P_to_from :Q_from_to, :Q_to_from),
            Tuple{String, String, String, Symbol, Float64, Float64, Float64, Float64},
        }[] =#
    connecting_branch_data = []
    branches = collect(PSY.get_components(x -> typeof(x) <: PSY.Branch, PSY.Component, sys))
    for branch in branches
        bus_number_from = PSY.get_number(PSY.get_from(PSY.get_arc(branch)))
        bus_number_to = PSY.get_number(PSY.get_to(PSY.get_arc(branch)))
        if (bus_number_from ∈ bus_numbers) && (bus_number_to ∉ bus_numbers)
            to_bus_name = PSY.get_name(PSY.get_to(PSY.get_arc(branch)))
            from_bus_name = PSY.get_name(PSY.get_from(PSY.get_arc(branch)))
            df_row = filter(row -> row.line_name == PSY.get_name(branch), powerflow_df)
            P_from_to = df_row.P_from_to[1]
            P_to_from = df_row.P_to_from[1]
            Q_from_to = df_row.Q_from_to[1]
            Q_to_from = df_row.Q_to_from[1]
            push!(
                connecting_branch_data,
                (
                    connecting_branch_name = PSY.get_name(branch),
                    to_bus_name = to_bus_name,
                    from_bus_name = from_bus_name,
                    inner_terminal = :from,
                    P_from_to = P_from_to,
                    P_to_from = P_to_from,
                    Q_from_to = Q_from_to,
                    Q_to_from = Q_to_from,
                ),
            )
        elseif (bus_number_from ∉ bus_numbers) && (bus_number_to ∈ bus_numbers)
            to_bus_name = PSY.get_name(PSY.get_to(PSY.get_arc(branch)))
            from_bus_name = PSY.get_name(PSY.get_from(PSY.get_arc(branch)))
            df_row = filter(row -> row.line_name == PSY.get_name(branch), powerflow_df)
            P_from_to = df_row.P_from_to[1]
            P_to_from = df_row.P_to_from[1]
            Q_from_to = df_row.Q_from_to[1]
            Q_to_from = df_row.Q_to_from[1]
            push!(
                connecting_branch_data,
                (
                    connecting_branch_name = PSY.get_name(branch),
                    to_bus_name = to_bus_name,
                    from_bus_name = from_bus_name,
                    inner_terminal = :to,
                    P_from_to = P_from_to,
                    P_to_from = P_to_from,
                    Q_from_to = Q_from_to,
                    Q_to_from = Q_to_from,
                ),
            )
        end
    end
    return connecting_branch_data
end
