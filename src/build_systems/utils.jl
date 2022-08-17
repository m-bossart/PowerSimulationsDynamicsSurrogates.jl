"""
create_subsystem_from_buses(sys_full::PSY.System, subsystem_bus_numbers::Vector{Int64}) 
-Takes a deepcopy of PSY system and removes all components that are not attached to the buses included in `subsystem_bus_numbers`. 
-Leaves all branches that connect two buses within the surrogate AND branches that connect to the rest of the system ("connecting buses")
-Adds source components at connecting buses --- these are used to perturb the surrogate when generating training data. 
"""
function create_subsystem_from_buses(
    sys_full::PSY.System,
    subsystem_bus_numbers::Vector{Int64},
)
    powerflow_df = PowerFlows.run_powerflow(sys_full)["flow_results"]
    sys = deepcopy(sys_full)

    for b in subsystem_bus_numbers
        bus = collect(PSY.get_components(PSY.Bus, sys, x -> PSY.get_number(x) == b))[1]
        connected_branches = collect(
            PSY.get_components(
                PSY.Component,
                sys,
                x ->
                    typeof(x) <: PSY.Branch && ((
                        PSY.get_from(PSY.get_arc(x)) == bus ||
                        PSY.get_to(PSY.get_arc(x)) == bus
                    )),
            ),
        )

        x = filter(
            x ->
                PSY.get_number(PSY.get_from(PSY.get_arc(x))) ∈ subsystem_bus_numbers &&
                    PSY.get_number(PSY.get_to(PSY.get_arc(x))) ∈ subsystem_bus_numbers,
            connected_branches,
        )
        if length(subsystem_bus_numbers) > 1 && length(x) == 0
            @error "The specified bus numbers do not give a continuous area"
            return false
        end
    end
    static_injectors = collect(
        PSY.get_components(PSY.Component, sys, x -> typeof(x) <: PSY.StaticInjection),
    )
    for static_injector in static_injectors
        if (PSY.get_number(PSY.get_bus(static_injector)) ∉ subsystem_bus_numbers)
            dynamic_injector = PSY.get_dynamic_injector(static_injector)
            (dynamic_injector !== nothing) && PSY.remove_component!(sys, dynamic_injector)
            PSY.remove_component!(sys, static_injector)
        end
    end
    connecting_bus_names = []
    connecting_bus_powers = []
    connecting_branch_names = Tuple{String, Symbol}[]
    branches = collect(PSY.get_components(PSY.Component, sys, x -> typeof(x) <: PSY.Branch))
    for branch in branches
        bus_number_from = PSY.get_number(PSY.get_from(PSY.get_arc(branch)))
        bus_number_to = PSY.get_number(PSY.get_to(PSY.get_arc(branch)))
        if (bus_number_from ∉ subsystem_bus_numbers) &&
           (bus_number_to ∉ subsystem_bus_numbers)
            arc = PSY.get_arc(branch)
            PSY.remove_component!(sys, arc)
            PSY.remove_component!(sys, branch)
        elseif (bus_number_from ∈ subsystem_bus_numbers) &&
               (bus_number_to ∉ subsystem_bus_numbers)
            connecting_bus_name = PSY.get_name(PSY.get_to(PSY.get_arc(branch)))
            df_row = filter(row -> row.line_name == PSY.get_name(branch), powerflow_df)
            P = df_row.P_to_from[1]
            push!(connecting_bus_names, connecting_bus_name)
            push!(connecting_bus_powers, P)
            push!(connecting_branch_names, (PSY.get_name(branch), :from))
        elseif (bus_number_from ∉ subsystem_bus_numbers) &&
               (bus_number_to ∈ subsystem_bus_numbers)
            connecting_bus_name = PSY.get_name(PSY.get_from(PSY.get_arc(branch)))
            df_row = filter(row -> row.line_name == PSY.get_name(branch), powerflow_df)
            P = df_row.P_from_to[1]
            push!(connecting_bus_names, connecting_bus_name)
            push!(connecting_bus_powers, P)
            push!(connecting_branch_names, (PSY.get_name(branch), :to))
        elseif (bus_number_from ∈ subsystem_bus_numbers) &&
               (bus_number_to ∈ subsystem_bus_numbers)
        else
            @error "Error in determining if branch should be removed"
        end
    end

    buses = collect(PSY.get_components(PSY.Component, sys, x -> typeof(x) <: PSY.Bus))
    for bus in buses
        if (PSY.get_number(bus) ∉ subsystem_bus_numbers) &&
           (PSY.get_name(bus) ∉ connecting_bus_names)
            PSY.remove_component!(sys, bus)
        end
    end

    reference_buses_in_surrogate = PSY.get_components(
        PSY.Bus,
        sys,
        x ->
            (PSY.get_bustype(x) == PSY.BusTypes.REF) &&
                (PSY.get_number(x) ∈ subsystem_bus_numbers),
    )
    for b in reference_buses_in_surrogate
        PSY.set_bustype!(b, PSY.BusTypes.PV)
    end

    for (i, connecting_bus_name) in enumerate(connecting_bus_names)
        bus = PSY.get_component(PSY.Bus, sys, connecting_bus_name)
        @assert bus !== nothing
        PSY.get_bustype(bus) !== PSY.BusTypes.REF && PSY.set_bustype!(bus, PSY.BusTypes.PV)
        source = PSY.Source(
            name = string("source_", i),
            active_power = connecting_bus_powers[i] / 100,
            available = true,
            reactive_power = 0.0,
            bus = bus,
            R_th = 1e-6, #0.0
            X_th = 1e-6, #5e-6
            internal_voltage = 0.0,
            internal_angle = 0.0,
        )
        PSY.add_component!(sys, source)
    end

    reference_buses_outside_surrogate = PSY.get_components(
        PSY.Bus,
        sys,
        x ->
            (PSY.get_bustype(x) == PSY.BusTypes.REF) &&
                (PSY.get_number(x) ∉ subsystem_bus_numbers),
    )
    if length(reference_buses_outside_surrogate) == 0
        first_bus = collect(
            PSY.get_components(
                PSY.Bus,
                sys,
                x -> (PSY.get_number(x) ∉ subsystem_bus_numbers),
            ),
        )[1]
        PSY.set_bustype!(first_bus, PSY.BusTypes.REF)
        @assert length(
            PSY.get_components(
                PSY.Bus,
                sys,
                x ->
                    (PSY.get_bustype(x) == PSY.BusTypes.REF) &&
                        (PSY.get_number(x) ∉ subsystem_bus_numbers),
            ),
        ) == 1
    end
    return sys, connecting_branch_names
end
