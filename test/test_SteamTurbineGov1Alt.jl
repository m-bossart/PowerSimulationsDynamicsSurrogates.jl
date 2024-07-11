@testset "Build and Execute Simulation with SteamTurbineGov1Alt" begin
    sys = build_system(PSIDSystems, "psid_11bus_andes")
    gen = get_component(ThermalStandard, sys, "generator-2-1")
    dyn_gen = get_component(DynamicGenerator, sys, "generator-2-1")
    pert_ref = PSID.ControlReferenceChange(1.0, dyn_gen, :P_ref, 0.02)
    sim = Simulation!(MassMatrixModel, sys, pwd(), (0.0, 10.0), pert_ref)
    @test sim.status == PSID.BUILT
    @test execute!(sim, Rodas5(); saveat = 0.01) == PSID.SIMULATION_FINALIZED
    results = read_results(sim)
    t, δ3 = get_state_series(results, ("generator-2-1", :δ))
    pm = get_prime_mover(dyn_gen)
    dyn_gen_new = DynamicGenerator(
        name = get_name(dyn_gen),
        ω_ref = get_ω_ref(dyn_gen),
        machine = get_machine(dyn_gen),
        shaft = get_shaft(dyn_gen),
        avr = get_avr(dyn_gen),
        prime_mover = SteamTurbineGov1Alt(
            R = get_R(pm),
            T1 = get_T1(pm),
            valve_position_limits = get_valve_position_limits(pm),
            T2 = get_T2(pm),
            T3 = get_T3(pm),
            D_T = get_D_T(pm),
            DB_h = get_DB_h(pm),
            DB_l = get_DB_l(pm),
            T_rate = get_T_rate(pm),
            P_ref = get_P_ref(pm),
        ),
        pss = get_pss(dyn_gen),
        base_power = get_base_power(dyn_gen),
    )
    remove_component!(sys, dyn_gen)
    add_component!(sys, dyn_gen_new, gen)
    pert_state = PSID.PerturbState(1.0, 63, (0.02 - 0.03897134141384735))
    sim = Simulation!(MassMatrixModel, sys, pwd(), (0.0, 10.0), pert_state)
    @test sim.status == PSID.BUILT
    @test execute!(sim, Rodas5(); saveat = 0.01) == PSID.SIMULATION_FINALIZED
    results = read_results(sim)
    t_alt, δ3_alt = get_state_series(results, ("generator-2-1", :δ))
    @test sum(abs.(δ3 .- δ3_alt)) < 0.02
    #using PlotlyJS
    #plot([scatter(x=t,y=δ3),scatter(x=t_alt,y=δ3_alt) ])
end
