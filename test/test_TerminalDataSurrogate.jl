
#This probably belongs in a different testset. It is not used now but might be useful in the future
#For defininig a solution prediction surrogate that corresponds 1:1 with an existing dynamic component. 
#= @testset "Test Initialization Function within the Surrogate" begin
    dyn_gen = DynamicInverter(
        name = "gen_perturb",
        ω_ref = 1.0,
        converter = converter_high_power(),
        outer_control = outer_control_droop(),
        inner_control = inner_control(),
        dc_source = dc_source_lv(),
        freq_estimator = pll(),
        filter = filt(),
    )
    _, refs = PSIDS.init_underlying_device(dyn_gen, 0.8, 0.2, 1.2, 0.11)
    @test refs["P_ref"] == 0.8047222222221517
    @test refs["Q_ref"] == 0.2944444444443753
    @test refs["V_ref"] == 1.300616306901241
end =#

@testset "Add TerminalDataSurrogate to system" begin
    sys = System(100)
    bus = ACBus(nothing)
    set_bustype!(bus, ACBusTypes.SLACK)
    add_component!(sys, bus)
    source = Source(nothing)
    set_bus!(source, bus)
    add_component!(sys, source)
    surrogate = TerminalDataSurrogate(nothing)
    add_component!(sys, surrogate, source)
    @test get_components(TerminalDataSurrogate, sys) !== nothing
end

@testset "Build and Execute Simulation with TerminalDataSurrogate" begin
    tspan = (0.0, 1.0)
    tstep = 0.01
    tsteps = tspan[1]:tstep:tspan[2]
    Random.seed!(1234)
    function gain(x)
        return x .* 10.0
    end
    v0_path = Lux.Chain(Lux.Dense(2, 2))
    i0_path = Lux.Chain(Lux.Dense(2, 2))
    v_path = Lux.Chain(
        Lux.FlattenLayer(),
        Lux.WrappedFunction(gain),
        Lux.Dense(10, 100),
        Lux.Dense(100, 100),
        Lux.Dense(100, 2),
    )
    i_path = Lux.Chain(Lux.FlattenLayer(), Lux.Dense(10, 100), Lux.Dense(100, 2))
    model = Lux.Chain(Lux.Parallel(+, v0_path, i0_path, v_path, i_path))

    rng = Random.default_rng()
    Random.seed!(rng, 0)
    ps, st = Lux.setup(rng, model)

    sys_full = build_system(PSIDSystems, "psid_11bus_andes")
    sys, _ = PSIDS.create_validation_system_from_buses(sys_full, [1, 2, 5, 6, 7, 8])
    show_components(sys, ACBus, [:number])
    for source in get_components(Source, sys)
        s = TerminalDataSurrogate(
            name = get_name(source),
            τ = 0.3,
            window_size = 5,
            fc = 0.0,
            steadystate_offset_correction = true,
            ext = Dict{String, Any}("model" => model, "ps" => ps, "st" => st),
        )
        display(s)
        add_component!(sys, s, source)
    end

    tspan = (0.0, 10.0)
    tfault = 5.0
    pert = PSID.BranchImpedanceChange(tfault, Line, "BUS 10-BUS 11-i_1", 1.5)

    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, pert)
    show_components(sys, Line)
    @test sim.status == PSID.BUILT
    @test execute!(
        sim,
        MethodOfSteps(Rodas5P(autodiff = false)),
        saveat = 0.1,
        reltol = 1e-8,
        abstol = 1e-11,
    ) == PSID.SIMULATION_FINALIZED
    results = read_results(sim)

    vbus = get_voltage_magnitude_series(results, 8)
    vr_surr = get_state_series(results, ("source_1", :vr))
    ir_surr = get_state_series(results, ("source_1", :ir))

    set_fc!(get_component(TerminalDataSurrogate, sys, "source_1"), 0.1)     #add filtering
    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, pert)
    show_components(sys, Line)
    @test sim.status == PSID.BUILT
    @test execute!(
        sim,
        MethodOfSteps(Rodas5P(autodiff = false)),
        saveat = 0.1,
        reltol = 1e-8,
        abstol = 1e-11,
    ) == PSID.SIMULATION_FINALIZED
    results = read_results(sim)

    vbus_filt = get_voltage_magnitude_series(results, 8)
    vr_surr_filt = get_state_series(results, ("source_1", :vr))
    ir_surr_filt = get_state_series(results, ("source_1", :ir))

    @test sum(vr_surr[2]) == 95.37524456676604
    @test sum(vr_surr_filt[2]) == 95.34139542343429

    #using PlotlyJS  #add PlotlyJS to test environment
    #display(
    #    PlotlyJS.plot([
    #        PlotlyJS.scatter(x = vbus[1], y = vbus[2], name = "no filtering"),
    #        PlotlyJS.scatter(x = vbus_filt[1], y = vbus_filt[2], name = "filter (0.1s)"),
    #    ]),
    #)
    #display(
    #    PlotlyJS.plot([
    #        PlotlyJS.scatter(x = vr_surr[1], y = vr_surr[2], name = "no filtering"),
    #        PlotlyJS.scatter(
    #            x = vr_surr_filt[1],
    #            y = vr_surr_filt[2],
    #            name = "filter (0.1s)",
    #        ),
    #    ]),
    #)
    #display(
    #    PlotlyJS.plot([
    #        PlotlyJS.scatter(x = ir_surr[1], y = ir_surr[2], name = "no filtering"),
    #        PlotlyJS.scatter(
    #            x = ir_surr_filt[1],
    #            y = ir_surr_filt[2],
    #            name = "filter (0.1s)",
    #        ),
    #    ]),
    #)
end

@testset "Test Basepower scaling" begin
    tspan = (0.0, 1.0)
    tstep = 0.01
    tsteps = tspan[1]:tstep:tspan[2]
    Random.seed!(1234)
    function gain(x)
        return x .* 10.0
    end
    v0_path = Lux.Chain(Lux.Dense(2, 2))
    i0_path = Lux.Chain(Lux.Dense(2, 2))
    v_path = Lux.Chain(
        Lux.FlattenLayer(),
        Lux.WrappedFunction(gain),
        Lux.Dense(10, 100),
        Lux.Dense(100, 100),
        Lux.Dense(100, 2),
    )
    i_path = Lux.Chain(Lux.FlattenLayer(), Lux.Dense(10, 100), Lux.Dense(100, 2))
    model = Lux.Chain(Lux.Parallel(+, v0_path, i0_path, v_path, i_path))

    rng = Random.default_rng()
    Random.seed!(rng, 0)
    ps, st = Lux.setup(rng, model)

    sys_full = build_system(PSIDSystems, "psid_11bus_andes")
    solve_powerflow(ACPowerFlow(), sys_full)
    sys, _ = PSIDS.create_validation_system_from_buses(sys_full, [1, 2, 5, 6, 7, 8])
    show_components(sys, ACBus, [:number])
    s = collect(get_components(Source, sys))[1]
    solve_ac_powerflow!(sys)
    #set_base_power!(s, 100.0/2)

    set_active_power!(s, get_active_power(s) / 2)
    set_reactive_power!(s, get_reactive_power(s) / 2)

    #@assert false 
    s_new = Source(
        name = "newname",
        available = get_available(s),
        bus = get_bus(s),
        active_power = get_active_power(s),
        reactive_power = get_reactive_power(s),
        R_th = get_R_th(s),
        X_th = get_X_th(s),
    )
    add_component!(sys, s_new)
    display(sys)
    sources = get_components(Source, sys)
    for source in sources
        s = TerminalDataSurrogate(
            name = get_name(source),
            τ = 0.3,
            window_size = 5,
            fc = 0.0,
            steadystate_offset_correction = true,
            ext = Dict{String, Any}("model" => model, "ps" => ps, "st" => st),
        )
        display(s)
        add_component!(sys, s, source)
        s_in_system = get_component(TerminalDataSurrogate, sys, get_name(source))

        set_base_power!(s_in_system, 100.0 / length(sources))
    end

    tspan = (0.0, 10.0)
    tfault = 5.0
    pert = PSID.BranchImpedanceChange(tfault, Line, "BUS 10-BUS 11-i_1", 1.5)

    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, pert)
    show_components(sys, Line)
    @test sim.status == PSID.BUILT
    @test execute!(
        sim,
        MethodOfSteps(Rodas5P(autodiff = false)),
        saveat = 0.1,
        reltol = 1e-8,
        abstol = 1e-11,
    ) == PSID.SIMULATION_FINALIZED
    results = read_results(sim)

    vbus = get_voltage_magnitude_series(results, 8)
    vr_surr = get_state_series(results, ("source_1", :vr))
    ir_surr = get_state_series(results, ("source_1", :ir))

    set_fc!(get_component(TerminalDataSurrogate, sys, "source_1"), 0.1)     #add filtering
    set_fc!(get_component(TerminalDataSurrogate, sys, "newname"), 0.1)     #add filtering

    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, pert)
    show_components(sys, Line)
    @test sim.status == PSID.BUILT
    @test execute!(
        sim,
        MethodOfSteps(Rodas5P(autodiff = false)),
        saveat = 0.1,
        reltol = 1e-8,
        abstol = 1e-11,
    ) == PSID.SIMULATION_FINALIZED
    results = read_results(sim)

    vbus_filt = get_voltage_magnitude_series(results, 8)
    vr_surr_filt = get_state_series(results, ("source_1", :vr))
    ir_surr_filt = get_state_series(results, ("source_1", :ir))

    @test sum(vr_surr[2]) == 95.37487158388667
    @test sum(vr_surr_filt[2]) == 95.34139540068229

    #using PlotlyJS  #add PlotlyJS to test environment
    #display(
    #    PlotlyJS.plot([
    #        PlotlyJS.scatter(x = vbus[1], y = vbus[2], name = "no filtering"),
    #        PlotlyJS.scatter(x = vbus_filt[1], y = vbus_filt[2], name = "filter (0.1s)"),
    #    ]),
    #)
    #display(
    #    PlotlyJS.plot([
    #        PlotlyJS.scatter(x = vr_surr[1], y = vr_surr[2], name = "no filtering"),
    #        PlotlyJS.scatter(
    #            x = vr_surr_filt[1],
    #            y = vr_surr_filt[2],
    #            name = "filter (0.1s)",
    #        ),
    #    ]),
    #)
    #display(
    #    PlotlyJS.plot([
    #        PlotlyJS.scatter(x = ir_surr[1], y = ir_surr[2], name = "no filtering"),
    #        PlotlyJS.scatter(
    #            x = ir_surr_filt[1],
    #            y = ir_surr_filt[2],
    #            name = "filter (0.1s)",
    #        ),
    #    ]),
    #)
end

@testset "Sensitivity - TerminalDataSurrogate" begin
    Random.seed!(1234)
    function gain(x)
        return x .* 10.0
    end
    v0_path = Lux.Chain(Lux.Dense(2, 2))
    i0_path = Lux.Chain(Lux.Dense(2, 2))
    v_path = Lux.Chain(Lux.FlattenLayer(), Lux.WrappedFunction(gain), Lux.Dense(10, 2))
    i_path = Lux.Chain(Lux.FlattenLayer(), Lux.Dense(10, 2))
    model = Lux.Chain(Lux.Parallel(+, v0_path, i0_path, v_path, i_path))
    Lux.f64(model)  #TODO - how to deal with Float32 vs Float64
    rng = Random.default_rng()
    Random.seed!(rng, 0)
    ps, st = Lux.setup(rng, model)
    sys_full = build_system(PSIDSystems, "psid_11bus_andes")
    sys, _ = PSIDS.create_validation_system_from_buses(sys_full, [1, 2, 5, 6, 7, 8])
    for source in get_components(Source, sys)
        s = TerminalDataSurrogate(
            name = get_name(source),
            τ = 0.3,
            window_size = 5,
            fc = 0.0,
            steadystate_offset_correction = true,
            ext = Dict{String, Any}("model" => model, "ps" => ps, "st" => st),
        )
        add_component!(sys, s, source)
    end
    tspan = (0.0, 1.0)
    tfault = 0.5
    statgen = get_component(ThermalStandard, sys, "generator-4-1")

    dyngen = get_component(DynamicGenerator, sys, "generator-4-1")
    #pert = PSID.BranchImpedanceChange(tfault, Line, "BUS 10-BUS 11-i_1", 1.5)  #not yet compatible with PSID sensitivity 
    pert = PSID.ControlReferenceChange(tfault, dyngen, :P_ref, 0.04)
    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, pert)

    #GET GROUND TRUTH DATA 
    execute!(
        sim,
        MethodOfSteps(Rodas5(autodiff = false));
        abstol = 1e-6,
        reltol = 1e-6,
        saveat = 0.5,
    )
    res = read_results(sim)
    t, δ_gt = get_state_series(res, ("generator-4-1", :δ))

    function plot_traces(δ, δ_gt)
        display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
    end
    EnzymeRules.inactive(::typeof(plot_traces), args...) = nothing
    function f_loss(p, states, δ_gt, aux)
        #plot_traces(states[1], δ_gt)
        return sum(abs.(states[1] - δ_gt))
    end
    sum(sim.inputs_init.ybus_rectangular)
    sum(sim.inputs.ybus_rectangular)

    sim.inputs.ybus_rectangular - sim.inputs_init.ybus_rectangular
    f_forward, f_grad, f_forward_zygote = get_sensitivity_functions(
        sim,
        [("source_1", :θ)],
        [("generator-4-1", :δ)],
        MethodOfSteps(Rodas5(autodiff = false)),
        f_loss;
        sensealg = ForwardDiffSensitivity(),
        abstol = 1e-2,
        reltol = 1e-2,
        saveat = 0.5,
    )
    θ = get_parameter_values(sim, [("source_1", :θ)])
    @test θ[1] == -0.03485802561044693

    @test isapprox(f_forward(θ, [pert], δ_gt, []), 0.0, atol = 1e-4)
    @test f_forward(θ * 1.1, [pert], δ_gt, []) ==
          f_forward_zygote(θ * 1.1, [pert], δ_gt, []) ==
          0.00014677407628382877
    @test f_grad(θ, [pert], δ_gt, [])[1] ==
          Zygote.gradient(p -> f_forward_zygote(p, [pert], δ_gt, []), θ)[1][1] ==
          9.26538894066859
end

@testset "Sensitivity - TerminalDataSurrogate - ReverseDiffAdjoint" begin
    Random.seed!(1234)
    function gain(x)
        return x .* 10.0
    end
    v0_path = Lux.Chain(Lux.Dense(2, 2))
    i0_path = Lux.Chain(Lux.Dense(2, 2))
    v_path = Lux.Chain(Lux.FlattenLayer(), Lux.WrappedFunction(gain), Lux.Dense(10, 2))
    i_path = Lux.Chain(Lux.FlattenLayer(), Lux.Dense(10, 2))
    model = Lux.Chain(Lux.Parallel(+, v0_path, i0_path, v_path, i_path))
    Lux.f64(model)  #TODO - how to deal with Float32 vs Float64
    #model = Lux.Experimental.@debug_mode model_1

    rng = Random.default_rng()
    Random.seed!(rng, 0)
    ps, st = Lux.setup(rng, model)
    sys_full = build_system(PSIDSystems, "psid_11bus_andes")
    sys, _ = PSIDS.create_validation_system_from_buses(sys_full, [1, 2, 5, 6, 7, 8])
    for source in get_components(Source, sys)
        s = TerminalDataSurrogate(
            name = get_name(source),
            τ = 0.3,
            window_size = 5,
            fc = 0.0,
            steadystate_offset_correction = true,
            ext = Dict{String, Any}("model" => model, "ps" => ps, "st" => st),
        )
        add_component!(sys, s, source)
        s_in_system = get_component(TerminalDataSurrogate, sys, get_name(source))
        set_base_power!(s_in_system, 100.0)
    end
    tspan = (0.0, 1.0)
    tfault = 0.6
    gen = get_component(ThermalStandard, sys, "generator-4-1")
    dyn_gen = get_component(DynamicGenerator, sys, "generator-4-1")

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
    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan)
    P_ref_prev = get_P_ref(dyn_gen_new)
    P_rev_new = 0.04
    state_index = PSID.make_global_state_map(sim.inputs)["generator-4-1"][:P_ref]
    pert = PSID.PerturbState(tfault, state_index, (P_rev_new - P_ref_prev))

    #pert = PSID.BranchImpedanceChange(tfault, Line, "BUS 10-BUS 11-i_1", 1.5)  #not yet compatible with PSID sensitivity 
    #pert = PSID.ControlReferenceChange(tfault, dyngen, :P_ref, 0.04)
    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, pert)

    #GET GROUND TRUTH DATA 
    execute!(
        sim,
        MethodOfSteps(Rodas5(autodiff = false));
        abstol = 1e-6,
        reltol = 1e-6,
        saveat = 0.1,
    )
    res = read_results(sim)
    t, δ_gt = get_state_series(res, ("generator-4-1", :δ))

    function plot_traces(δ, δ_gt)
        display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
    end
    EnzymeRules.inactive(::typeof(plot_traces), args...) = nothing
    function f_loss(p, states, δ_gt, aux)
        #plot_traces(states[1], δ_gt)
        return sum(abs.(states[1] - δ_gt))
    end
    sum(sim.inputs_init.ybus_rectangular)
    sum(sim.inputs.ybus_rectangular)

    sim.inputs.ybus_rectangular - sim.inputs_init.ybus_rectangular
    f_forward, f_grad, f_forward_zygote = get_sensitivity_functions(
        sim,
        [("source_1", :θ)],
        [("generator-4-1", :δ)],
        MethodOfSteps(Rodas5(autodiff = false)),
        f_loss;
        sensealg = ForwardDiffSensitivity(),
        abstol = 1e-2,
        reltol = 1e-2,
        saveat = 0.1,
    )
    θ = get_parameter_values(sim, [("source_1", :θ)])
    @test θ[1] == -0.03485802561044693

    @test isapprox(f_forward(θ, [pert], δ_gt, []), 0.0, atol = 1e-4)
    @test f_forward(θ * 1.1, [pert], δ_gt, []) ==
          f_forward_zygote(θ * 1.1, [pert], δ_gt, []) ==
          0.00010788599295485923
    grads_forward = Zygote.gradient(p -> f_forward_zygote(p, [pert], δ_gt, []), θ)[1]
    @test f_grad(θ, [pert], δ_gt, [])[1] == grads_forward[1] == -170.10717294865816

    f_forward, f_grad, f_forward_zygote = get_sensitivity_functions(
        sim,
        [("source_1", :θ)],
        [("generator-4-1", :δ)],
        MethodOfSteps(Rodas5(autodiff = false)),
        f_loss;
        sensealg = ReverseDiffAdjoint(),
        abstol = 1e-2,
        reltol = 1e-2,
        saveat = 0.1,
    )

    @test Zygote.gradient(p -> f_forward_zygote(p, [pert], δ_gt, []), θ * 1.001)[1][30] ==
          2.4040948119363748e-5
    @test Zygote.gradient(p -> f_forward_zygote(p, [pert], δ_gt, []), θ * 2.0)[1][30] ==
          0.0004445025697350502
end
