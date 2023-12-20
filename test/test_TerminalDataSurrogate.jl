
#This probably belongs in a different testset. It is not used now but might be useful in the future
#For defininig a solution prediction surrogate that corresponds 1:1 with an existing dynamic component. 
@testset "Test Initialization Function within the Surrogate" begin
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
end

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
    sys = System("test/data_tests/ThreeBus.raw")

    add_source_to_ref(sys)
    gen = collect(get_components(ThermalStandard, sys))[1]
    dyn_gen = dyn_gen_second_order(gen)
    add_component!(sys, dyn_gen, gen)
    b = collect(get_components(x -> get_number(x) == 102, Bus, sys))[1]
    source = Source(;
        name = "test",
        available = true,
        bus = b,
        active_power = 0.2,
        reactive_power = 0.1,
        R_th = 0.0,
        X_th = 0.0,
        internal_voltage = 1.0,
        internal_angle = 0.0,
    )
    add_component!(sys, source)
    v0_path = Lux.Chain(Lux.Dense(2, 2))
    i0_path = Lux.Chain(Lux.Dense(2, 2))
    v_path = Lux.Chain(Lux.FlattenLayer(), Lux.Dense(10, 2))
    i_path = Lux.Chain(Lux.FlattenLayer(), Lux.Dense(10, 2))
    model = Lux.Chain(Lux.Parallel(+, v0_path, i0_path, v_path, i_path))
    rng = Random.default_rng()
    Random.seed!(rng, 0)
    ps, st = Lux.setup(rng, model)
    s = TerminalDataSurrogate(
        name = "test",
        τ = 0.1,
        window_size = 5,
        ext = Dict{String, Any}("model" => model, "ps" => ps, "st" => st),
    )
    add_component!(sys, s, source)

    #@show PowerFlows.solve_powerflow(PowerFlows.ACPowerFlow(), sys)["bus_results"]
    #@show PowerFlows.solve_powerflow(PowerFlows.ACPowerFlow(), sys)["flow_results"]

    tspan = (0.0, 5.0)
    pert = PSID.BranchImpedanceChange(0.1, Line, "BUS 1-BUS 2-i_1", 1.5)

    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, pert)
    @test sim.status == PSID.BUILT
    @test execute!(sim, MethodOfSteps(Rodas5(autodiff = false))) ==
          PSID.SIMULATION_FINALIZED
    results = read_results(sim)

    v = get_voltage_magnitude_series(results, 101)
end
