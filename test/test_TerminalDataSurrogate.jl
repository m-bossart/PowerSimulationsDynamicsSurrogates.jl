@testset "Add TerminalDataSurrogate to system" begin
    sys = System(100)
    bus = Bus(nothing)
    set_bustype!(bus, BusTypes.SLACK)
    add_component!(sys, bus)
    surrogate = TerminalDataSurrogate(nothing)
    set_bus!(surrogate, bus)
    add_component!(sys, surrogate)
    @test get_components(TerminalDataSurrogate, sys) !== nothing
end

@testset "Build and Execute Simulation with TerminalDataSurrogate @ zero output" begin
    tspan = (0.0, 1.0)
    tstep = 0.01
    tsteps = tspan[1]:tstep:tspan[2]
    Random.seed!(1234)
    sys = System("test/data_tests/ThreeBus.raw")
    display(sys)
    add_source_to_ref(sys)
    gen = collect(get_components(ThermalStandard, sys))[1]
    dyn_gen = dyn_gen_second_order(gen)
    add_component!(sys, dyn_gen, gen)

    b = collect(get_components(x -> get_number(x) == 102, Bus, sys))[1]
    m = FullyConnected(
        nn_structure = [(18, 2, false, "tanh")],
        nn_parameters = zeros(36),
        input_min = ones(18) .* -1,
        input_max = ones(18),
        input_lims = (-1.0, 1.0),
        target_min = [-1.0, -1.0],
        target_max = [1.0, 1.0],
        target_lims = (-1.0, 1.0),
    )
    s = TerminalDataSurrogate(
        name = "test",
        available = true,
        bus = b,
        active_power = 0.0,
        reactive_power = 0.0,
        active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = 0.0, max = 1.0),
        internal_voltage = 1.0,
        internal_angle = 0.0,
        model_architecture = m,
        underlying_dynamic_model = GenericDER(nothing),
        θ_ref_frame = 0.0,
        n_past_timesteps = 3,   # size of input: 5*n_past_timesteps + 3
    )
    PowerFlows.solve_powerflow(PowerFlows.ACPowerFlow(), sys)["bus_results"]
    add_component!(sys, s)
    tspan = (0.0, 1.0)
    sim = Simulation!(ResidualModel, sys, pwd(), tspan)
    @test sim.status == PSID.BUILT
    @test execute!(sim, IDA()) == PSID.SIMULATION_FINALIZED
    results = read_results(sim)
end

@testset "Add callbacks for caching values" begin
    tspan = (0.0, 1.0)
    tstep = 0.01
    tsteps = tspan[1]:tstep:tspan[2]
    Random.seed!(1234)
    sys = System("test/data_tests/ThreeBus.raw")
    display(sys)
    add_source_to_ref(sys)
    gen = collect(get_components(ThermalStandard, sys))[1]
    dyn_gen = dyn_gen_second_order(gen)
    add_component!(sys, dyn_gen, gen)

    b = collect(get_components(x -> get_number(x) == 102, Bus, sys))[1]
    m = FullyConnected(
        nn_structure = [(18, 2, false, "tanh")],
        nn_parameters = zeros(36),
        input_min = ones(18) .* -1,
        input_max = ones(18),
        input_lims = (-1.0, 1.0),
        target_min = [-1.0, -1.0],
        target_max = [1.0, 1.0],
        target_lims = (-1.0, 1.0),
    )
    s = TerminalDataSurrogate(
        name = "test",
        available = true,
        bus = b,
        active_power = 0.0,
        reactive_power = 0.0,
        active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = 0.0, max = 1.0),
        internal_voltage = 1.0,
        internal_angle = 0.0,
        model_architecture = m,
        underlying_dynamic_model = GenericDER(nothing),
        θ_ref_frame = 0.0,
        n_past_timesteps = 3,   # size of input: 5*n_past_timesteps + 3
    )
    PowerFlows.solve_powerflow(PowerFlows.ACPowerFlow(), sys)["bus_results"]
    add_component!(sys, s)
    cbs = PSID.Perturbation[]
    for s in PSY.get_components(TerminalDataSurrogate, sys)
        push!(cbs, SolutionSurrogateCacheValues(s))
    end
    tspan = (0.0, 1.0)
    sim = Simulation!(ResidualModel, sys, pwd(), tspan, cbs)
    @test sim.status == PSID.BUILT
    @test execute!(sim, IDA()) == PSID.SIMULATION_FINALIZED
    results = read_results(sim)
end

@testset "Time performance hit" begin
    tspan = (0.0, 1.0)
    tstep = 0.01
    tsteps = tspan[1]:tstep:tspan[2]
    Random.seed!(1234)
    sys = System("test/data_tests/144Bus.json")
    gentrip = GeneratorTrip(0.1, get_component(DynamicInjection, sys, "GFM_Battery_2"))
    tspan = (0.0, 1.0)
    sim = Simulation!(ResidualModel, sys, pwd(), tspan, gentrip)
    @test sim.status == PSID.BUILT
    @test execute!(sim, IDA()) == PSID.SIMULATION_FINALIZED
    results_no_surrogates = read_results(sim)
    solve_time_no_surrogates = results_no_surrogates.time_log[:timed_solve_time]

    for (ix, b) in
        enumerate(get_components(x -> PSY.get_bustype(x) == BusTypes.PV, PSY.Bus, sys))
        m = FullyConnected(
            nn_structure = [(18, 2, false, "tanh")],
            nn_parameters = zeros(36),
            input_min = ones(18) .* -1,
            input_max = ones(18),
            input_lims = (-1.0, 1.0),
            target_min = [-1.0, -1.0],
            target_max = [1.0, 1.0],
            target_lims = (-1.0, 1.0),
        )
        s = TerminalDataSurrogate(
            name = string("test", ix),
            available = true,
            bus = b,
            active_power = 0.0,
            reactive_power = 0.0,
            active_power_limits = (min = 0.0, max = 1.0),
            reactive_power_limits = (min = 0.0, max = 1.0),
            internal_voltage = 1.0,
            internal_angle = 0.0,
            model_architecture = m,
            underlying_dynamic_model = GenericDER(nothing),
            θ_ref_frame = 0.0,
            n_past_timesteps = 3,   # size of input: 5*n_past_timesteps + 3
        )
        add_component!(sys, s)
    end
    display(sys)
    cbs = PSID.Perturbation[]
    for s in PSY.get_components(TerminalDataSurrogate, sys)
        push!(cbs, SolutionSurrogateCacheValues(s))
    end
    push!(cbs, gentrip)
    sim = Simulation!(ResidualModel, sys, pwd(), tspan, cbs)
    @test sim.status == PSID.BUILT
    @test execute!(sim, IDA()) == PSID.SIMULATION_FINALIZED
    results_surrogates = read_results(sim)
    solve_time_surrogates = results_surrogates.time_log[:timed_solve_time]
    println("solve time without surrogates: ", solve_time_no_surrogates)
    println("solve time wtih surrogates: ", solve_time_surrogates)
    @test solve_time_surrogates < solve_time_no_surrogates * 2  #TODO - make tighter
    #p = plot(get_voltage_magnitude_series(results_no_surrogates, 2), label ="no surrogates")
    #plot!(p, get_voltage_magnitude_series(results_surrogates, 2), label ="surrogates")
    #display(p)
end
