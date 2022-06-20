@testset "Add SteadyStateNODE to system" begin
    sys = System(100)
    bus = Bus(nothing)
    set_bustype!(bus, BusTypes.SLACK)
    add_component!(sys, bus)
    source = Source(nothing)
    set_bus!(source, bus)
    add_component!(sys, source)
    ssnode = SteadyStateNODE(nothing)
    add_component!(sys, ssnode, source)
    @test get_components(SteadyStateNODE, sys) !== nothing
end

@testset "Compare to flux + execute simulation - Untrained Params" begin
    Random.seed!(1234)
    sys = System("test/data_tests/ThreeBus.raw")
    add_source_to_ref(sys)
    #display(run_powerflow(sys)["flow_results"])
    tspan = (0.0, 1.0)
    tstep = 0.01
    tsteps = tspan[1]:tstep:tspan[2]
    x_scale = [1.1, 1.0, 1.0, 1.0]
    x_bias = [0.0, 0.0, -0.1, 0.0]
    exogenous_scale = [1.1, 1.0]
    exogenous_bias = [0.0, 0.1]

    initializer = Chain((x) -> x .* x_scale .+ x_bias, Dense(4, 3, tanh; bias = true))

    node = Chain(
        Parallel(
            +,
            Chain(
                (x) -> x .* exogenous_scale .+ exogenous_bias,
                Dense(2, 3, tanh; bias = true),
            ),
            Dense(3, 3, tanh; bias = true),
        ),
        Dense(3, 3, tanh; bias = true),
    )
    observer = Dense(3, 2, tanh; bias = true)

    function SteadyStateNODE_simple(source)
        return SteadyStateNODE(
            name = get_name(source),
            initializer_structure = [(4, 3, true, "tanh")],
            initializer_parameters = Flux.destructure(initializer)[1],
            node_structure_exogenous = [(2, 3, true, "tanh")],
            node_structure_states = [(3, 3, true, "tanh")],
            node_structure_common = [(3, 3, true, "tanh")],
            node_parameters = Flux.destructure(node)[1],
            observer_structure = [(3, 2, true, "tanh")],
            observer_parameters = Flux.destructure(observer)[1],
            x_scale = x_scale,
            x_bias = x_bias,
            exogenous_scale = exogenous_scale,
            exogenous_bias = exogenous_bias,
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
    x_flux = rand(4)
    x_psid = copy(x_flux)
    ex_flux = rand(2)
    ex_psid = copy(ex_flux)
    r_flux = rand(3)
    r_psid = copy(r_flux)

    x_scaled = PowerSimulationsDynamicsSurrogates._x_scale(surrogate_wrapper, x_psid)
    ex_scaled =
        PowerSimulationsDynamicsSurrogates._exogenous_scale(surrogate_wrapper, ex_psid)
    @test isapprox(
        initializer(x_flux),
        PowerSimulationsDynamicsSurrogates._forward_pass_initializer(
            surrogate_wrapper,
            x_scaled,
        );
        atol = 1e-14,
    )
    @test isapprox(
        node((ex_flux, r_flux)),
        PowerSimulationsDynamicsSurrogates._forward_pass_node(
            surrogate_wrapper,
            r_psid,
            ex_scaled,
        );
        atol = 1e-14,
    )
    @test isapprox(
        observer(r_flux),
        PowerSimulationsDynamicsSurrogates._forward_pass_observer(
            surrogate_wrapper,
            r_psid,
        );
        atol = 1e-14,
    )
    @test surrogate_wrapper.ext["epsilon"] == [-0.03433842411754684, 1.2399628875880429]
    @test surrogate_wrapper.ext["initializer_error"] ==
          [-1.6297954841407067, 4.169328812631354, 0.19831182702606404]

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
