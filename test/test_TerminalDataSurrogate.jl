
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

    #@show PowerFlows.solve_powerflow(PowerFlows.ACPowerFlow(), sys)["bus_results"]
    #@show PowerFlows.solve_powerflow(PowerFlows.ACPowerFlow(), sys)["flow_results"]

    tspan = (0.0, 10.0)
    tfault = 5.0
    pert = PSID.BranchImpedanceChange(tfault, Line, "BUS 10-BUS 11-i_1", 1.5)

    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, pert)
    show_components(sys, Line)
    @test sim.status == PSID.BUILT
    @test execute!(sim, MethodOfSteps(Rodas5(autodiff = false))) ==
          PSID.SIMULATION_FINALIZED
    results = read_results(sim)

    vbus = get_voltage_magnitude_series(results, 8)
    vr_surr = get_state_series(results, ("source_1", :vr))
    ir_surr = get_state_series(results, ("source_1", :ir))

    set_fc!(get_component(TerminalDataSurrogate, sys, "source_1"), 0.1)     #add filtering
    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, pert)
    show_components(sys, Line)
    @test sim.status == PSID.BUILT
    @test execute!(sim, MethodOfSteps(Rodas5(autodiff = false))) ==
          PSID.SIMULATION_FINALIZED
    results = read_results(sim)

    vbus_filt = get_voltage_magnitude_series(results, 8)
    vr_surr_filt = get_state_series(results, ("source_1", :vr))
    ir_surr_filt = get_state_series(results, ("source_1", :ir))

    using PlotlyJS  #add PlotlyJS to test environment
    display(
        PlotlyJS.plot([
            PlotlyJS.scatter(x = vbus[1], y = vbus[2], name = "no filtering"),
            PlotlyJS.scatter(x = vbus_filt[1], y = vbus_filt[2], name = "filter (0.1s)"),
        ]),
    )
    display(
        PlotlyJS.plot([
            PlotlyJS.scatter(x = vr_surr[1], y = vr_surr[2], name = "no filtering"),
            PlotlyJS.scatter(
                x = vr_surr_filt[1],
                y = vr_surr_filt[2],
                name = "filter (0.1s)",
            ),
        ]),
    )
    display(
        PlotlyJS.plot([
            PlotlyJS.scatter(x = ir_surr[1], y = ir_surr[2], name = "no filtering"),
            PlotlyJS.scatter(
                x = ir_surr_filt[1],
                y = ir_surr_filt[2],
                name = "filter (0.1s)",
            ),
        ]),
    )
end

@testset "Sensitivity - TerminalDataSurrogate" begin
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
    tspan = (0.0, 10.0)
    tfault = 5.0
    statgen = get_component(ThermalStandard, sys, "generator-4-1")

    dyngen = get_component(DynamicGenerator, sys, "generator-4-1")
    pert = PSID.BranchImpedanceChange(tfault, Line, "BUS 10-BUS 11-i_1", 1.5)
    pert = PSID.ControlReferenceChange(tfault, dyngen, :P_ref, 0.04)
    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, pert)

    #GET GROUND TRUTH DATA 
    execute!(
        sim,
        MethodOfSteps(Rodas5(autodiff = false));
        abstol = 1e-9,
        reltol = 1e-9,
        dtmax = 0.005,
        saveat = 0.005,
    )
    res = read_results(sim)
    t, δ_gt = get_state_series(res, ("generator-4-1", :δ))

    function f_loss(δ, δ_gt)
        #display(plot([scatter(;x = t, y = δ_gt), scatter(;x= t, y = δ)]))
        #display(plot([scatter(;x = t, y = δ_gt .- δ)]))
        return sum(abs.(δ - δ_gt))
    end
    sum(sim.inputs_init.ybus_rectangular)
    sum(sim.inputs.ybus_rectangular)

    sim.inputs.ybus_rectangular - sim.inputs_init.ybus_rectangular
    f_forward, f_grad = get_sensitivity_functions(
        sim,
        [("source_1", :θ)],
        [("generator-4-1", :δ)],
        MethodOfSteps(Rodas5(autodiff = false)),
        f_loss;
        sensealg = ForwardDiffSensitivity(),
        abstol = 1e-9,
        reltol = 1e-9,
        dtmax = 0.005,
        saveat = 0.005,
    )
    using ComponentArrays
    p = ComponentArray(ps)
    loss_zero = f_forward(p, δ_gt)
    @test isapprox(loss_zero, 0.0, atol = 5e-3)

    #TODO - perturbation is not being cleared properly for branch impedance change?
    #TODO - get non-zero losses when changing the parameters
    #TODO - try the gradient once Enzyme is made compatible with DDE problem: https://github.com/SciML/DiffEqBase.jl/issues/1061
end
