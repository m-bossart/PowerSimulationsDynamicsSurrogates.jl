@testset "one bus system" begin
    sys = System("test/data_tests/OneBus.raw")
    slack_bus = [b for b in PSY.get_components(Bus, sys) if b.bustype == ACBusTypes.REF][1]
    inf_source = Source(
        name = "InfBus", #name
        available = true, #availability
        active_power = -0.5,
        reactive_power = 0.0,
        bus = slack_bus, #bus
        R_th = 0.0,
        X_th = 5e-6, #Xth
    )
    PSY.add_component!(sys, inf_source)
    for s in get_components(ThermalStandard, sys)
        dyngen = dyn_gen_second_order(s)
        add_component!(sys, dyngen, s)
    end

    d = collect(get_components(DynamicInjection, sys))[1]
    s = collect(get_components(Source, sys))[1]
    #p = ControlReferenceChange(0.5, d, :P_ref, 1.01) 
    p = SourceBusVoltageChange(0.5, s, :V_ref, 1.01)
    sim = Simulation!(
        MassMatrixModel,
        sys,
        pwd(),
        (0.0, 1.0),
        [p];
        initialize_simulation = true,
    )
    show_states_initial_value(sim)
    using OrdinaryDiffEq
    execute!(sim, Rodas5(), saveat = 0:0.01:1.0)
    res = read_results(sim)
    show_states_initial_value(res)
    v1 = PSID.get_voltage_magnitude_series(res, 101)
    solve_powerflow(PowerFlows.ACPowerFlow(), sys)["bus_results"]
    @test v1[2][1] == 1.05
end
