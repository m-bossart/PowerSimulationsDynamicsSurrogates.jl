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

@testset "Compare to Flux" begin
    x_scale = [1.1, 1.0]
    x_bias = [0.0, 1.0]
    initializer = Chain((x) -> x .* x_scale .+ x_bias, Dense(2, 10, tanh; bias = true))
    node = Chain(
        Parallel(+, Dense(2, 2; bias = true), Dense(2, 2; bias = true)),
        Dense(2, 2; bias = true),
    )
    observer = Dense(2, 2; bias = true)

    ssnode = SteadyStateNODE(
        name = "test",
        initializer_structure = [(2, 10, true, "tanh")],
        initializer_parameters = Flux.destructure(initializer)[1],
        node_structure_exogenous = [(2, 2, true, "tanh")],
        node_structure_states = [(2, 2, true, "tanh")],
        node_structure_common = [(2, 2, true, "tanh")],
        node_parameters = Flux.destructure(node)[1],
        observer_structure = [(2, 2, true, "tanh")],
        observer_parameters = Flux.destructure(observer)[1],
        x_scale = x_scale,
        x_bias = x_bias,
        exogenous_scale = [1.0, 1.0],
        exogenous_bias = [0.0, 0.0],
    )
    x = rand(2)
    ex = rand(2)
    r = rand(2)
    #wrapper = get_wrapper_from_sim()       #TODO - get the wrapper from simulation! 

    #@test isapprox(initializer(x), _forward_pass_initializer(wrapper, x); atol = 1e-14)     
    #@test isapprox(node((r,ex)), _forward_pass_node(wrapper, r, ex); atol = 1e-14)    
    #@test isapprox(observer(x), _forward_pass_observer(wrapper, r); atol = 1e-14)     

end

@testset "Execute Simulations" begin
    sys = System("test/data_tests/TwoBus.raw")
    add_sources_to_buses(sys)
    tspan = (0.0, 1.0)
    tstep = 0.01
    tsteps = tspan[1]:tstep:tspan[2]
    x_scale = [1.0, 1.0, 1.0, 1.0]
    x_bias = [0.0, 0.0, 0.0, 0.0]
    initializer = Chain((x) -> x .* x_scale .+ x_bias, Dense(4, 10, tanh; bias = true))
    node = Chain(
        Parallel(+, Dense(2, 2; bias = true), Dense(2, 2; bias = true)),
        Dense(2, 2; bias = true),
    )
    observer = Dense(2, 2; bias = true)

    function SteadyStateNODE_simple(source)
        return SteadyStateNODE(
            name = get_name(source),
            initializer_structure = [(4, 2, true, "tanh")],
            initializer_parameters = zeros(length(Flux.destructure(initializer)[1])),   #TODO- debug initialization with non-zero parameters
            node_structure_exogenous = [(2, 2, true, "tanh")],
            node_structure_states = [(2, 2, true, "tanh")],
            node_structure_common = [(2, 2, true, "tanh")],
            node_parameters = zeros(length(Flux.destructure(node)[1])),
            observer_structure = [(2, 2, true, "tanh")],
            observer_parameters = zeros(length(Flux.destructure(observer)[1])),
            x_scale = x_scale,
            x_bias = x_bias,
            exogenous_scale = [1.0, 1.0],
            exogenous_bias = [0.0, 0.0],
        )
    end

    source_surrogate = [s for s in get_components(Source, sys)][1]
    source_IB = [s for s in get_components(Source, sys)][2]
    ssnode = SteadyStateNODE_simple(source_surrogate)
    add_component!(sys, ssnode, source_surrogate)

    sim = Simulation!(
        ResidualModel,
        sys,
        pwd(),
        tspan,
        SourceBusVoltageChange(0.5, source_IB, :V_ref, 1.04),
    )
    @test execute!(sim, IDA(), saveat = tsteps) == PSID.SIMULATION_FINALIZED
    results = read_results(sim)
    p = plot()
    for b in get_components(Bus, sys)
        plot!(p, PSID.get_voltage_magnitude_series(results, get_number(b)))
    end
    display(p)
end
