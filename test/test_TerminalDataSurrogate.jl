
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
    scaler = MinMaxScaler(
        scale_input = true,
        input_min = ones(8) .* -1,
        input_max = ones(8),
        input_lims = (-1.0, 1.0),
        scale_target = true,
        target_min = [-1.0, -1.0],
        target_max = [1.0, 1.0],
        target_lims = (-1.0, 1.0),
    )
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
    s = TerminalDataSurrogate(
        name = "test",
        available = true,
        bus = b,
        active_power = 0.2,
        reactive_power = 0.1,
        active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = 0.0, max = 1.0),
        internal_voltage = 1.0,
        internal_angle = 0.0,
        model_architecture = [FFNN(18, 2, false, "tanh")],
        model_parameters = zeros(36),
        underlying_dynamic_model = dyn_gen,
        data_scaler = scaler,
        θ_ref_frame = 0.0,
        n_past_timesteps = 3,
    )
    add_component!(sys, s)
    display(sys)
    @show PowerFlows.solve_powerflow(PowerFlows.ACPowerFlow(), sys)["bus_results"]

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
    scaler = MinMaxScaler(
        scale_input = true,
        input_min = ones(8) .* -1,
        input_max = ones(8),
        input_lims = (-1.0, 1.0),
        scale_target = true,
        target_min = [-1.0, -1.0],
        target_max = [1.0, 1.0],
        target_lims = (-1.0, 1.0),
    )
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
        model_architecture = [FFNN(18, 2, false, "tanh")],
        model_parameters = zeros(36),
        underlying_dynamic_model = dyn_gen,
        data_scaler = scaler,
        θ_ref_frame = 0.0,
        n_past_timesteps = 3,   # size of input: 5*n_past_timesteps + 3
    )
    PowerFlows.solve_powerflow(PowerFlows.ACPowerFlow(), sys)["bus_results"]
    add_component!(sys, s)
    cbs = PSID.Perturbation[]
    for s in PSY.get_components(TerminalDataSurrogate, sys)
        push!(cbs, TerminalDataSurrogateCacheValues(s))
    end
    tspan = (0.0, 1.0)
    sim = Simulation!(ResidualModel, sys, pwd(), tspan, cbs)
    @test sim.status == PSID.BUILT
    @test execute!(sim, IDA()) == PSID.SIMULATION_FINALIZED
    results = read_results(sim)
end

#= @testset "Time performance hit" begin
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
        scaler = MinMaxScaler(
            scale_input = true,
            input_min = ones(8) .* -1,
            input_max = ones(8),
            input_lims = (-1.0, 1.0),
            scale_target = true,
            target_min = [-1.0, -1.0],
            target_max = [1.0, 1.0],
            target_lims = (-1.0, 1.0),
        )
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
            model_architecture = [FFNN(18, 2, false, "tanh")],
            model_parameters = zeros(36),
            underlying_dynamic_model = dyn_gen,
            data_scaler = scaler,
            θ_ref_frame = 0.0,
            n_past_timesteps = 3,
        )
        add_component!(sys, s)
    end
    display(sys)
    cbs = PSID.Perturbation[]
    for s in PSY.get_components(TerminalDataSurrogate, sys)
        push!(cbs, TerminalDataSurrogateCacheValues(s))
    end
    push!(cbs, gentrip)
    sim = Simulation!(ResidualModel, sys, pwd(), tspan, cbs)
    @test sim.status == PSID.BUILT
    @test execute!(sim, IDA()) == PSID.SIMULATION_FINALIZED
    results_surrogates = read_results(sim)
    solve_time_surrogates = results_surrogates.time_log[:timed_solve_time]
    println("solve time without surrogates: ", solve_time_no_surrogates)
    println("solve time wtih surrogates: ", solve_time_surrogates)
    @test solve_time_surrogates < solve_time_no_surrogates * 2.5  #TODO - make tighter
    #p = plot(get_voltage_magnitude_series(results_no_surrogates, 2), label ="no surrogates")
    #plot!(p, get_voltage_magnitude_series(results_surrogates, 2), label ="surrogates")
    #display(p)
end
 =#
#TODO - add testset for different types of model architectures. 
#TODO - test serialization/deserialization 
