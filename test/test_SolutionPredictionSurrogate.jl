@testset "Add SolutionPredictionSurrogate to system" begin
    sys = System(100)
    bus = Bus(nothing)
    set_bustype!(bus, BusTypes.SLACK)
    add_component!(sys, bus)
    surrogate = SolutionPredictionSurrogate(nothing)
    set_bus!(surrogate, bus)
    add_component!(sys, surrogate)
    @test get_components(SolutionPredictionSurrogate, sys) !== nothing
end

@testset "Build and Execute Simulation with SolutionPredictionSurrogate @ zero output" begin
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
    s = SolutionPredictionSurrogate(
        name = "test",
        available = true,
        bus = b,
        active_power = 0.0,
        reactive_power = 0.0,
        active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = 0.0, max = 1.0),
        internal_voltage = 1.0,
        internal_angle = 0.0,
        nn_structure = [(8, 2, false, "tanh")],
        nn_parameters = zeros(16),
        nn_type = :dense,
        length_cache = 1,
        nn_features = :direct,
        input_min = -1 .* ones(8),
        input_max = ones(8),
        input_lims = (-1.0, 1.0),
        target_min = -1 .* ones(2),
        target_max = ones(2),
        target_lims = (-1.0, 1.0),
    )
    PowerFlows.run_powerflow(sys)["bus_results"]
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
    s = SolutionPredictionSurrogate(
        name = "test",
        available = true,
        bus = b,
        active_power = 0.0,
        reactive_power = 0.0,
        active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = 0.0, max = 1.0),
        internal_voltage = 1.0,
        internal_angle = 0.0,
        nn_structure = [(13, 2, false, "tanh")],
        nn_parameters = zeros(26),
        nn_type = :dense,
        length_cache = 2,
        nn_features = :direct,
        input_min = -1 .* ones(13),
        input_max = ones(13),
        input_lims = (-1.0, 1.0),
        target_min = -1 .* ones(2),
        target_max = ones(2),
        target_lims = (-1.0, 1.0),
    )
    PowerFlows.run_powerflow(sys)["bus_results"]
    add_component!(sys, s)
    cbs = PSID.Perturbation[]
    for s in PSY.get_components(SolutionPredictionSurrogate, sys)
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

    for b in get_components(x -> PSY.get_bustype(x) == BusTypes.PV, PSY.Bus, sys)
        s = SolutionPredictionSurrogate(
            name = string("test", PSY.get_name(b)),
            available = true,
            bus = b,
            active_power = 0.0,
            reactive_power = 0.0,
            active_power_limits = (min = 0.0, max = 1.0),
            reactive_power_limits = (min = 0.0, max = 1.0),
            internal_voltage = 1.0,
            internal_angle = 0.0,
            nn_structure = [(13, 2, false, "tanh")],   #Make this a more realistic size NN 
            nn_parameters = zeros(26),
            nn_type = :dense,
            length_cache = 2,                       #Make the length of cache 
            nn_features = :direct,
            input_min = -1 .* ones(13),
            input_max = ones(13),
            input_lims = (-1.0, 1.0),
            target_min = -1 .* ones(2),
            target_max = ones(2),
            target_lims = (-1.0, 1.0),
        )
        add_component!(sys, s)
    end
    display(sys)
    cbs = PSID.Perturbation[]
    for s in PSY.get_components(SolutionPredictionSurrogate, sys)
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
    @test solve_time_surrogates < solve_time_no_surrogates * 1.5  #TODO - make tighter
    #p = plot(get_voltage_magnitude_series(results_no_surrogates, 2), label ="no surrogates")
    #plot!(p, get_voltage_magnitude_series(results_surrogates, 2), label ="surrogates")
    #display(p)
end

#Goals:Friday
#2) Write a performance test. Time a medium size system without any surrogates. Add a bunch of surrogates all over the place which output zero.
#Compare tstops for the two simulations (should be the same)
#Compare the run time (hopefully not to big of a performance hit)

# With random parameterization, see if we can build and execute a simulation.
# Somehow we want to check that the caching of prior values is correct... 

#= @testset "Build Simulation with SolutionPredictionSurrogate" begin
    Random.seed!(1234)
    sys = System("test/data_tests/ThreeBus.raw")
    add_source_to_ref(sys)
    #display(run_powerflow(sys)["flow_results"])
    tspan = (0.0, 1.0)
    tstep = 0.01
    tsteps = tspan[1]:tstep:tspan[2]
    input_min = [0.0, 0.0]
    input_max = [1.0, 1.0]
    input_lims = (-1.0, 1.0)
    target_min = [0.0, 0.0]
    target_max = [1.0, 1.0]
    target_lims = (-1.0, 1.0)

    #THESE LAYERS DO THE SCALING WITHIN THE FLUX NNs
    layer_initializer_input = Parallel(
        vcat,
        (x) -> PowerSimulationsDynamicsSurrogates.min_max_normalization(
            x,
            input_min[2],
            input_max[2],
            input_lims[2],
            input_lims[1],
        ),
        (x) -> PowerSimulationsDynamicsSurrogates.min_max_normalization(
            x,
            target_min,
            target_max,
            target_lims[2],
            target_lims[1],
        ),
    )  #Second output only 
    layer_node_input = Parallel(
        vcat,
        (x) -> x,
        (x) -> PowerSimulationsDynamicsSurrogates.min_max_normalization(
            x,
            input_min,
            input_max,
            input_lims[2],
            input_lims[1],
        ),
        (x) -> x,
    )

    initializer = Chain(layer_initializer_input, Dense(3, 5, tanh; bias = true))
    node = Chain(layer_node_input, Dense(7, 3, tanh; bias = true))

    function SteadyStateNODE_simple(source)
        return SteadyStateNODE(
            name = get_name(source),
            initializer_structure = [(3, 5, true, "tanh")],
            initializer_parameters = Flux.destructure(initializer)[1],
            node_structure = [(7, 3, true, "tanh")],
            node_parameters = Flux.destructure(node)[1],
            input_min = input_min,
            input_max = input_max,
            input_lims = input_lims,
            target_min = target_min,
            target_max = target_max,
            target_lims = target_lims,
        )
    end

    source_surrogate = [s for s in get_components(Source, sys)][1]
    ssnode = SteadyStateNODE_simple(source_surrogate)
    add_component!(sys, ssnode, source_surrogate)
    for g in PSY.get_components(Generator, sys)
        dyn_g = dyn_gen_second_order(g)
        add_component!(sys, dyn_g, g)
    end

    sim = Simulation!(
        ResidualModel,
        sys,
        pwd(),
        tspan,
        PSID.BranchTrip(0.5, PSY.Line, "BUS 1-BUS 2-i_1"),
    )
    surrogate_wrapper = filter(
        x -> typeof(x) == PSID.DynamicWrapper{SteadyStateNODE},
        sim.inputs.dynamic_injectors,
    )[1]

    target_flux = rand(2)
    target_psid = copy(target_flux)
    input_flux = rand(2)
    input_psid = copy(input_flux)
    r_flux = rand(3)
    r_psid = copy(r_flux)
    ref_flux = rand(2)
    ref_psid = copy(ref_flux)

    target_scaled =
        PowerSimulationsDynamicsSurrogates._target_scale(surrogate_wrapper, target_psid)
    input_scaled =
        PowerSimulationsDynamicsSurrogates._input_scale(surrogate_wrapper, input_psid)
    @test isapprox(
        initializer((input_flux[2], target_flux)),
        PowerSimulationsDynamicsSurrogates._forward_pass_initializer(
            surrogate_wrapper,
            input_scaled[2],
            target_scaled,
        );
        atol = 1e-14,
    )

    @test isapprox(
        node((r_flux, input_flux, ref_flux)),
        PowerSimulationsDynamicsSurrogates._forward_pass_node(
            surrogate_wrapper,
            r_psid,
            input_scaled,
            ref_psid,
        );
        atol = 1e-14,
    )

    @test isapprox(
        surrogate_wrapper.ext["initializer_error"],
        [
            0.8620012947569395,
            -1.5949281018523056,
            0.4659432191773033,
            0.3382842845487015,
            -1.0813953248504833,
        ],
        atol = 1e-10,
    )
    @test execute!(sim, IDA(), saveat = tsteps) == PSID.SIMULATION_FINALIZED
    results = read_results(sim)

    #Plot results - for debug only 
    p = plot()
    for b in get_components(Bus, sys)
        plot!(
            p,
            PSID.get_voltage_magnitude_series(results, get_number(b)),
            label = string(get_number(b)),
        )
    end
    r1 = PSID.get_state_series(results, ("InfBus", :r1))
    p2 = plot()
    plot!(p2, r1, label = "r1")
    display(plot(p, p2))
end
 =#
