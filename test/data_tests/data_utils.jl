function add_source_to_ref(sys::PSY.System)
    for g in PSY.get_components(StaticInjection, sys)
        isa(g, ElectricLoad) && continue
        g.bus.bustype == ACBusTypes.REF &&
            error("A device is already attached to the REF bus")
    end

    slack_bus = [b for b in PSY.get_components(Bus, sys) if b.bustype == ACBusTypes.REF][1]
    inf_source = Source(
        name = "InfBus", #name
        available = true, #availability
        active_power = 0.0,
        reactive_power = 0.0,
        bus = slack_bus, #bus
        R_th = 0.0,
        X_th = 5e-6, #Xth
    )
    PSY.add_component!(sys, inf_source)
    return
end

function add_sources_to_buses(sys::PSY.System)
    for g in PSY.get_components(StaticInjection, sys)
        isa(g, ElectricLoad) && continue
        g.bus.bustype == ACBusTypes.REF &&
            error("A device is already attached to the REF bus")
    end

    buses = PSY.get_components(PSY.ACBus, sys)
    display(buses)
    for (i, b) in enumerate(buses)
        inf_source = Source(
            name = string("InfBus", i), #name
            available = true, #availability
            active_power = 0.0,
            reactive_power = 0.0,
            bus = b, #bus
            R_th = 0.0,
            X_th = 5e-6, #Xth
        )
        PSY.add_component!(sys, inf_source)
    end
    return
end
