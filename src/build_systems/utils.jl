"""
    create_surrogate_training_system(sys_full::PSY.System, surrogate_bus_numbers::Vector{Int64}) 
-Takes a deepcopy of PSY system and removes all components that are not attached to the buses included in `surrogate_bus_numbers`. 
-Leaves all branches that connect two buses within the surrogate AND branches that connect to the rest of the system ("connecting buses")
-Adds source components at connecting buses --- these are used to perturb the surrogate when generating training data. 
"""
#TODO - don't bother with the surrogate area name part - just list the bus numbers. 
#Later, when doing testing you will want to do the opposite -- pass the opposite list of bus numbers? 
#New name - "create_subsystem_from_buses(sys_full, bus_numbers)"
#TODO - should we be setting a bus as the reference? Is there a general rule for this? 
function create_subsystem_from_buses(sys_full::PSY.System, surrogate_bus_numbers::Vector{Int64})
    powerflow_df = PowerFlows.run_powerflow(sys_full)["flow_results"]
    sys = deepcopy(sys_full)
    static_injectors = collect(
        PSY.get_components(PSY.Component, sys, x -> typeof(x) <: PSY.StaticInjection),
    )
    for static_injector in static_injectors
        if (PSY.get_number(PSY.get_bus(static_injector)) ∉ surrogate_bus_numbers) 
            dynamic_injector = PSY.get_dynamic_injector(static_injector)
            (dynamic_injector !== nothing) && PSY.remove_component!(sys, dynamic_injector)
            PSY.remove_component!(sys, static_injector)
        end
    end
    connecting_bus_names = []
    connecting_bus_powers= [] 
    branches = collect(PSY.get_components(PSY.Component, sys, x -> typeof(x) <: PSY.Branch))
    for branch in branches
        bus_number_from = PSY.get_number(PSY.get_from(PSY.get_arc(branch)))
        bus_number_to = PSY.get_number(PSY.get_to(PSY.get_arc(branch)))
        if (bus_number_from ∉  surrogate_bus_numbers) && (bus_number_to ∉ surrogate_bus_numbers )  
            arc = PSY.get_arc(branch)
            PSY.remove_component!(sys, arc)
            PSY.remove_component!(sys, branch)
        elseif (bus_number_from ∈ surrogate_bus_numbers ) && (bus_number_to ∉ surrogate_bus_numbers )
            connecting_bus_name = PSY.get_name(PSY.get_to(PSY.get_arc(branch)))
            push!(connecting_bus_names, connecting_bus_name)
            #TODO - set Powers (same as below )
        elseif  (bus_number_from ∉ surrogate_bus_numbers ) && (bus_number_to ∈ surrogate_bus_numbers )
            connecting_bus_name = PSY.get_name(PSY.get_from(PSY.get_arc(branch)))
            display(powerflow_df)
            display(branch)
            df_row = filter(row -> row.line_name == PSY.get_name(branch), powerflow_df)
            P = df_row.P_to_from[1] #TODO - opposite of expectation... 
            #P = powerflow_df[(powerflow_df.line_name == connecting_bus_name), :] # [connecting_bus_name, !] df[(df.A .> 500) .& (300 .< df.C .< 400), :]
            display(P)
            push!(connecting_bus_names, connecting_bus_name)
            push!(connecting_bus_powers, P )
        elseif  (bus_number_from ∈ surrogate_bus_numbers ) && (bus_number_to ∈ surrogate_bus_numbers )
            @info "Branch is within surrogate, leaving"
        else 
            @error "Error in determining if branch should be removed"
        end 
    end

    buses = collect(PSY.get_components(PSY.Component, sys, x -> typeof(x) <: PSY.Bus))
    for bus in buses
        if (PSY.get_number(bus) ∉ surrogate_bus_numbers) &&
           (PSY.get_name(bus) ∉ connecting_bus_names)
            PSY.remove_component!(sys, bus)
        end
    end

    #TODO - if there is exactly one reference bus outside the surrogate don't do anything. 
    #If the reference bus is in the surrogate, error
    #If there are no references bus, set one of the connecting buses to be reference.

    #PSY.set_bustype!(bus, PSY.BusTypes.REF)    #TODO - determine the appropriate action here
    #PSY.set_bustype!(bus, PSY.BusTypes.PV)     


    @warn connecting_bus_names
    for (i,connecting_bus_name) in enumerate(connecting_bus_names)
        bus = PSY.get_component(PSY.Bus, sys, connecting_bus_name)
        #for branch in PSY.get_components(PSY.Branch)
        @assert bus !== nothing 

        source = PSY.Source(
            name = string("source_", i),
            active_power = connecting_bus_powers[i],                    #Only need to set this appropriately (if it is a PV bus...) 
            available = true,
            reactive_power = 0.0,
            bus = bus,
            R_th = 1e-6, #0.0
            X_th = 1e-6, #5e-6
            internal_voltage = 0.0,     #Voltage and angle will come from the Bus when system is initialized. 
            internal_angle = 0.0,
        )
        PSY.add_component!(sys, source)
    end 

    for b in PSY.get_components(PSY.Bus, sys)
        #PSY.set_bustype!(b, PSY.BusTypes.PV) 
        @warn PSY.get_bustype(b)
    end 

    for s in PSY.get_components(PSY.Source, sys)
     
        @warn PSY.get_active_power(s)
    end 
#=     connecting_branches = find_connecting_branches(sys, surrogate_area_name)
    for (fault_id, PVSDatas) in pvs_data
        for PVSData in PVSDatas
            branch = PSY.get_components_by_name(PSY.ACBranch, sys, PVSData.branch_name)
            @assert length(branch) == 1
            branch = branch[1]
            branch_name = PSY.get_name(branch)
            if PVSData.from_or_to == "from"
                bus_to_add_pvs = PSY.get_from(PSY.get_arc(branch))
            elseif PVSData.from_or_to == "to"
                bus_to_add_pvs = PSY.get_to(PSY.get_arc(branch))
            else
                @error "invalid value of from_or_to"
            end
            @warn bus_to_add_pvs
            pvs_name = string(fault_id, "_", branch_name)

            PSY.set_bustype!(bus_to_add_pvs, PSY.BusTypes.REF)
            source = PSY.Source(
                name = pvs_name,
                active_power = PVSData.P0,
                available = false,
                reactive_power = PVSData.Q0,
                bus = bus_to_add_pvs,
                R_th = 1e-6, #0.0
                X_th = 1e-6, #5e-6
                internal_voltage = PVSData.V0,
                internal_angle = PVSData.θ0,
            )
            @warn "Voltage Bias when building PVS", PVSData.internal_voltage_bias
            pvs = PSY.PeriodicVariableSource(
                name = PSY.get_name(source),
                R_th = PSY.get_R_th(source),
                X_th = PSY.get_X_th(source),
                internal_voltage_bias = PVSData.internal_voltage_bias,
                internal_voltage_frequencies = PVSData.internal_voltage_frequencies,
                internal_voltage_coefficients = PVSData.internal_voltage_coefficients,
                internal_angle_bias = PVSData.internal_angle_bias,
                internal_angle_frequencies = PVSData.internal_angle_frequencies,
                internal_angle_coefficients = PVSData.internal_angle_coefficients,
            )
            PSY.add_component!(sys, source)
            PSY.add_component!(sys, pvs, source)
        end
    end =#
    return sys
end

#= 
function label_area!(sys::PSY.System, bus_numbers, area_name::String)
    buses = collect(PSY.get_components(PSY.Bus, sys))
    areas = collect(PSY.get_components(PSY.Area, sys))
    for area in areas
        if PSY.get_name(area) == area_name
            @error "area already exists"
            return 0
        end
    end
    surrogate_area = PSY.Area(; name = area_name)
    PSY.add_component!(sys, surrogate_area)
    for bus in buses
        if PSY.get_number(bus) in bus_numbers
            PSY.set_area!(bus, surrogate_area)
        end
    end
end =#