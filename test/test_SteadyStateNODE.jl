@testset "Add SteadyStateNODE to system" begin
    sys = System(100)
    bus = ACBus(nothing)
    set_bustype!(bus, ACBusTypes.SLACK)
    add_component!(sys, bus)
    source = Source(nothing)
    set_bus!(source, bus)
    add_component!(sys, source)
    ssnode = SteadyStateNODE(nothing)
    add_component!(sys, ssnode, source)
    @test get_components(SteadyStateNODE, sys) !== nothing
end

@testset "Execute simulation - Untrained Params" begin
    sys = System("test/data_tests/ThreeBus.raw")
    add_source_to_ref(sys)
    #display(run_powerflow(sys)["flow_results"])
    tspan = (0.0, 1.0)
    tstep = 0.01
    tsteps = tspan[1]:tstep:tspan[2]
    input_min = [0.0, 0.0]
    input_max = [1.0, 1.0]
    input_lims = (0.0, 1.0)
    target_min = [0.0, 0.0]
    target_max = [1.0, 1.0]
    target_lims = (0.0, 1.0)

    layer_initializer_input = Lux.Parallel(
        vcat,
        WrappedFunction(
            (x) -> PowerSimulationsDynamicsSurrogates.min_max_normalization(
                x,
                input_min[2],
                input_max[2],
                input_lims[2],
                input_lims[1],
            ),
        ),
        WrappedFunction(
            (x) -> PowerSimulationsDynamicsSurrogates.min_max_normalization(
                x,
                target_min,
                target_max,
                target_lims[2],
                target_lims[1],
            ),
        ),
    )  #Second output only 
    layer_node_input = Lux.Parallel(
        vcat,
        WrappedFunction((x) -> x),
        WrappedFunction(
            (x) -> PowerSimulationsDynamicsSurrogates.min_max_normalization(
                x,
                input_min,
                input_max,
                input_lims[2],
                input_lims[1],
            ),
        ),
        WrappedFunction((x) -> x),
    )

    initializer = Lux.f64(Lux.Chain(layer_initializer_input, Lux.Dense(3, 5, tanh)))
    node = Lux.f64(Lux.Chain(layer_node_input, Lux.Dense(7, 3, tanh)))
    rng = Random.default_rng()
    Random.seed!(rng, 1234)
    ps_initializer, st_initializer = Lux.setup(rng, initializer)
    Random.seed!(rng, 3)
    ps_node, st_node = Lux.setup(rng, node)
    @test isapprox(
        [0.39178863167762756, -0.3014172315597534],
        ComponentArray(ps_initializer)[1:2],
        atol = 1e-7,
    )
    @test isapprox(
        [0.2890770435333252, 0.37868717312812805],
        ComponentArray(ps_node)[1:2],
        atol = 1e-7,
    )

    source_surrogate = [s for s in get_components(Source, sys)][1]
    ssnode = SteadyStateNODE(
        name = get_name(source_surrogate),
        n_states = 3,
        ext = Dict{String, Any}(
            "model_node" => node,
            "model_init" => initializer,
            "ps_node" => ps_node,
            "ps_init" => ps_initializer,
            "st_node" => st_node,
            "st_init" => st_initializer,
        ),
    )
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

    @test isapprox(
        surrogate_wrapper.ext["initializer_error"],
        [
            -0.2624221835960246,
            -1.375540778144861,
            -0.2441003958223837,
            1.297850317193571,
            -0.011542250287076644,
        ],
        atol = 1e-7,
    )

    @test execute!(sim, IDA(), saveat = tsteps) == PSID.SIMULATION_FINALIZED
    results = read_results(sim)
    t, δ = get_state_series(results, ("generator-102-1", :δ))
    @test isapprox(δ[1], 0.7050620746521667, atol = 1e-7)
    @test isapprox(δ[end], 0.5345286521782711, atol = 2e-7)
    #Plot results - for debug only 
    #=     p = plot()
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
        display(plot(p, p2)) =#
end
EnzymeRules.inactive(::typeof(PowerSystems._get_multiplier), args...) = nothing #TODO - don't rely on inactivating get multiplier! 
@testset "Sensitivity function- Untrained Params" begin
    sys = System("test/data_tests/ThreeBus.raw")
    add_source_to_ref(sys)
    tspan = (0.0, 2.0)
    tstep = 0.01
    tsteps = tspan[1]:tstep:tspan[2]
    input_min = [0.0, 0.0]
    input_max = [1.0, 1.0]
    input_lims = (0.0, 1.0)
    target_min = [0.0, 0.0]
    target_max = [1.0, 1.0]
    target_lims = (0.0, 1.0)

    layer_initializer_input =
        Lux.Parallel(vcat, WrappedFunction((x) -> x), WrappedFunction((x) -> x))  #Second output only 
    layer_node_input = Lux.Parallel(
        vcat,
        WrappedFunction((x) -> x),
        WrappedFunction((x) -> x),
        WrappedFunction((x) -> x),
    )

    initializer = Lux.f64(Lux.Chain(layer_initializer_input, Lux.Dense(3, 5, tanh)))
    node = Lux.f64(Lux.Chain(layer_node_input, Lux.Dense(7, 3, tanh)))
    rng = Random.default_rng()
    Random.seed!(rng, 1234)
    ps_initializer, st_initializer = Lux.setup(rng, initializer)
    Random.seed!(rng, 3)
    ps_node, st_node = Lux.setup(rng, node)

    source_surrogate = [s for s in get_components(Source, sys)][1]
    ssnode = SteadyStateNODE(
        name = get_name(source_surrogate),
        n_states = 3,
        ext = Dict{String, Any}(
            "model_node" => node,
            "model_init" => initializer,
            "ps_node" => ps_node,
            "ps_init" => ps_initializer,
            "st_node" => st_node,
            "st_init" => st_initializer,
        ),
    )
    add_component!(sys, ssnode, source_surrogate)
    for g in PSY.get_components(Generator, sys)
        dyn_g = dyn_gen_second_order(g)
        add_component!(sys, dyn_g, g)
    end
    dyn_g = collect(PSY.get_components(DynamicGenerator, sys))[1]
    pert = ControlReferenceChange(0.5, dyn_g, :P_ref, 0.799999)

    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, pert)

    @test execute!(sim, Rodas5(), saveat = 0.005) == PSID.SIMULATION_FINALIZED
    results = read_results(sim)
    t, δ = get_state_series(results, ("generator-102-1", :δ))

    function plot_traces(δ, δ_gt)
        display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
    end
    EnzymeRules.inactive(::typeof(plot_traces), args...) = nothing
    function f_loss(p, δ, δ_gt, aux)
        #plot_traces(δ[1], δ_gt)
        return sum(abs.(δ[1] - δ_gt))
    end

    f_forward, f_grad, f_forward_zygote = get_sensitivity_functions(
        sim,
        [("InfBus", :θ, "node")],
        [("generator-102-1", :δ)],
        Rodas5(),
        f_loss;
        sensealg = ForwardDiffSensitivity(),
        abstol = 1e-6,
        reltol = 1e-6,
        #dtmax = 0.005,
        saveat = 0.005,
    )
    θ = get_parameter_values(sim, [("InfBus", :θ, "node")])
    @test θ[1] == 0.2890770435333252
    @test isapprox(f_forward(θ, [pert], δ, []), 0.0, atol = 1e-4)
    rand_scaler = [rand() / 100.0 + 1.0 for x in 1:length(θ)]
    @test isapprox(f_forward(θ .* rand_scaler, [pert], δ, []), 0.0, atol = 1e-4)
    @test f_forward(θ .* rand_scaler, [pert], δ, []) ==
          f_forward_zygote(θ .* rand_scaler, [pert], δ, [])
    @test f_grad(θ * 1.001, [pert], δ, [])[1] ==
          Zygote.gradient(p -> f_forward_zygote(p, [pert], δ, []), θ * 1.001)[1][1] ==
          -0.0044582109185284935
end

@testset "Sensitivity function- ReverseDiffAdjoint" begin
    sys = System("test/data_tests/ThreeBus.raw")
    add_source_to_ref(sys)
    tspan = (0.0, 2.0)
    tstep = 0.01
    tsteps = tspan[1]:tstep:tspan[2]
    input_min = [0.0, 0.0]
    input_max = [1.0, 1.0]
    input_lims = (0.0, 1.0)
    target_min = [0.0, 0.0]
    target_max = [1.0, 1.0]
    target_lims = (0.0, 1.0)

    layer_initializer_input =
        Lux.Parallel(vcat, WrappedFunction((x) -> x), WrappedFunction((x) -> x))  #Second output only 
    layer_node_input = Lux.Parallel(
        vcat,
        WrappedFunction((x) -> x),
        WrappedFunction((x) -> x),
        WrappedFunction((x) -> x),
    )

    initializer = Lux.f64(Lux.Chain(layer_initializer_input, Lux.Dense(3, 5, tanh)))
    node = Lux.f64(Lux.Chain(layer_node_input, Lux.Dense(7, 3, tanh)))
    rng = Random.default_rng()
    Random.seed!(rng, 1234)
    ps_initializer, st_initializer = Lux.setup(rng, initializer)
    Random.seed!(rng, 3)
    ps_node, st_node = Lux.setup(rng, node)

    source_surrogate = [s for s in get_components(Source, sys)][1]
    ssnode = SteadyStateNODE(
        name = get_name(source_surrogate),
        n_states = 3,
        ext = Dict{String, Any}(
            "model_node" => node,
            "model_init" => initializer,
            "ps_node" => ps_node,
            "ps_init" => ps_initializer,
            "st_node" => st_node,
            "st_init" => st_initializer,
        ),
    )
    add_component!(sys, ssnode, source_surrogate)
    for g in PSY.get_components(Generator, sys)
        dyn_g = dyn_gen_second_order(g)
        add_component!(sys, dyn_g, g)
    end
    gen = get_component(ThermalStandard, sys, "generator-102-1")

    dyn_gen = collect(PSY.get_components(DynamicGenerator, sys))[1]
    pm = get_prime_mover(dyn_gen)
    dyn_gen_new = DynamicGenerator(
        name = get_name(dyn_gen),
        ω_ref = get_ω_ref(dyn_gen),
        machine = get_machine(dyn_gen),
        shaft = get_shaft(dyn_gen),
        avr = get_avr(dyn_gen),
        prime_mover = TGFixedAlt(efficiency = get_efficiency(pm), P_ref = get_P_ref(pm)),
        pss = get_pss(dyn_gen),
        base_power = get_base_power(dyn_gen),
    )
    remove_component!(sys, dyn_gen)
    add_component!(sys, dyn_gen_new, gen)
    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan)
    P_ref_prev = get_P_ref(dyn_gen_new)
    P_rev_new = 0.799999
    state_index = PSID.make_global_state_map(sim.inputs)["generator-102-1"][:P_ref]
    pert_state = PSID.PerturbState(1.0, state_index, (P_rev_new - P_ref_prev))

    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, pert_state)
    @test execute!(sim, Rodas5(), saveat = 0.005) == PSID.SIMULATION_FINALIZED

    results = read_results(sim)
    t, δ = get_state_series(results, ("generator-102-1", :δ))

    function plot_traces(δ, δ_gt)
        display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
    end
    EnzymeRules.inactive(::typeof(plot_traces), args...) = nothing
    function f_loss(p, δ, δ_gt, aux)
        #plot_traces(δ[1], δ_gt)
        return sum(abs.(δ[1] - δ_gt))
    end

    f_forward, f_grad, f_forward_zygote = get_sensitivity_functions(
        sim,
        [("InfBus", :θ, "node")],
        [("generator-102-1", :δ)],
        Rodas5(autodiff = false),
        f_loss;
        sensealg = ForwardDiffSensitivity(),  #TODO - fails with ReverseDiffAdjoint() and autodiff=true
        abstol = 1e-6,
        reltol = 1e-6,
        #dtmax = 0.005,
        saveat = 0.005,
    )
    θ = get_parameter_values(sim, [("InfBus", :θ, "node")])
    @test θ[1] == 0.2890770435333252
    @test isapprox(f_forward(θ, [pert_state], δ, []), 0.0, atol = 1e-4)
    rand_scaler = [rand() / 100.0 + 1.0 for x in 1:length(θ)]
    @test isapprox(f_forward(θ .* rand_scaler, [pert_state], δ, []), 0.0, atol = 1e-4)
    @test f_forward(θ .* rand_scaler, [pert_state], δ, []) ==
          f_forward_zygote(θ .* rand_scaler, [pert_state], δ, [])
    grads_forward = Zygote.gradient(p -> f_forward_zygote(p, [pert_state], δ, []), θ * 1.001)[1]
    @test f_grad(θ * 1.001, [pert_state], δ, [])[1] == grads_forward[1] == 0.0005295818045851775
    _, _, f_forward_zygote = get_sensitivity_functions(
        sim,
        [("InfBus", :θ, "node")],
        [("generator-102-1", :δ)],
        Rodas5(autodiff = false),         #TODO - fails with ReverseDiffAdjoint() and autodiff=true
        f_loss;
        sensealg = ReverseDiffAdjoint(),
        abstol = 1e-6,
        reltol = 1e-6,
        #dtmax = 0.005,
        saveat = 0.005,
    )
    grads_reverse = Zygote.gradient(p -> f_forward_zygote(p, [pert_state], δ, []), θ * 1.001)[1]
    pert_state_new = PSID.PerturbState(0.5, state_index, (P_rev_new - P_ref_prev))
    grads_reverse_2 =
        Zygote.gradient(p -> f_forward_zygote(p, [pert_state_new], δ, []), θ * 1.001)[1]
    grads_reverse_3 =
        Zygote.gradient(p -> f_forward_zygote(p, [pert_state], δ, []), θ * 1.001)[1]
    @test grads_reverse[1] .- grads_reverse_3[1] == 0.0  #repeatable for same computation
    @test grads_reverse[1] .- grads_reverse_2[1] != 0.0  #time of perturbation impacts gradient
    #plt = plot()
    #add_trace!(plt, scatter(y=grads_forward))
    #add_trace!(plt, scatter(y=grads_reverse))
    #display(plt)   #Compare gradient from ForwardDiffSensitivity() and ReverseDiffAdjoint()
    @test LinearAlgebra.norm(abs.(grads_reverse .- grads_forward)) < 0.007
end
