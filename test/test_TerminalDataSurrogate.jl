
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

    add_source_to_ref(sys)
    gen = collect(get_components(ThermalStandard, sys))[1]
    dyn_gen = dyn_gen_second_order(gen)
    add_component!(sys, dyn_gen, gen)
    b = collect(get_components(x -> get_number(x) == 102, Bus, sys))[1]
    v0_path = Lux.Chain(Lux.Dense(2, 2))
    i0_path = Lux.Chain(Lux.Dense(2, 2))
    v_path = Lux.Chain(Lux.FlattenLayer(), Lux.Dense(10,2))
    i_path = Lux.Chain(Lux.FlattenLayer(), Lux.Dense(10,2))
    model = Lux.Chain(Lux.Parallel(+, v0_path, i0_path, v_path, i_path))
    rng = Random.default_rng()
    Random.seed!(rng, 0)
    ps, st = Lux.setup(rng, model) 

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
        τ = 0.1, 
        window_size = 5, 
        ext = Dict{String, Any}("model"=>model, "ps"=>ps, "st"=> st),
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
#TODO - test serialization/deserialization; in a different test set.... 