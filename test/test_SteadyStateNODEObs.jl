@testset "Add SteadyStateNODEObs to system" begin
    sys = System(100)
    bus = ACBus(nothing)
    set_bustype!(bus, ACBusTypes.SLACK)
    add_component!(sys, bus)
    source = Source(nothing)
    set_bus!(source, bus)
    add_component!(sys, source)
    ssnode = SteadyStateNODEObs(nothing)
    add_component!(sys, ssnode, source)
    @test get_components(SteadyStateNODEObs, sys) !== nothing
end

@testset "Compare to flux + execute simulation - Untrained Params" begin
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
    layer_initializer_input = Flux.Parallel(
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
    layer_node_input = Flux.Parallel(
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
    layer_observer_output =
        (x) -> PowerSimulationsDynamicsSurrogates.min_max_normalization_inverse(
            x,
            target_min,
            target_max,
            target_lims[2],
            target_lims[1],
        )

    initializer =
        Flux.f64(Flux.Chain(layer_initializer_input, Flux.Dense(3, 5, tanh; bias = true)))
    node = Flux.f64(Flux.Chain(layer_node_input, Flux.Dense(7, 3, tanh; bias = true)))
    observer =
        Flux.f64(Flux.Chain(Flux.Dense(3, 2, tanh; bias = true), layer_observer_output))

    function SteadyStateNODEObs_simple(source)
        return SteadyStateNODEObs(
            name = get_name(source),
            initializer_structure = [(3, 5, true, "tanh")],
            initializer_parameters = Flux.destructure(initializer)[1],
            node_structure = [(7, 3, true, "tanh")],
            node_parameters = Flux.destructure(node)[1],
            observer_structure = [(3, 2, true, "tanh")],
            observer_parameters = Flux.destructure(observer)[1],
            input_min = input_min,
            input_max = input_max,
            input_lims = input_lims,
            target_min = target_min,
            target_max = target_max,
            target_lims = target_lims,
        )
    end

    source_surrogate = [s for s in get_components(Source, sys)][1]
    ssnode = SteadyStateNODEObs_simple(source_surrogate)
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
        x -> typeof(x) == PSID.DynamicWrapper{SteadyStateNODEObs},
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
        observer(r_flux),
        PowerSimulationsDynamicsSurrogates._target_scale_inverse(
            surrogate_wrapper,
            PowerSimulationsDynamicsSurrogates._forward_pass_observer(
                surrogate_wrapper,
                r_psid,
            ),
        );
        atol = 1e-14,
    )

    @test isapprox(
        surrogate_wrapper.ext["initializer_error"],
        [
            3.7711905762469193,
            -3.1373439536687755,
            1.3036265514829255,
            -0.969121363441726,
            -5.0936732894054035,
        ],
        atol = 1e-10,
    )
    @test execute!(sim, IDA(), saveat = tsteps) == PSID.SIMULATION_FINALIZED
    results = read_results(sim)

    #=     #Plot results - for debug only 
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
        display(plot(p, p2)) =#
end
