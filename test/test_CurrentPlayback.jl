include(joinpath(dirname(@__FILE__), "data_tests/dynamic_test_data.jl"))

function add_source_to_ref(sys::PSY.System)
    for g in PSY.get_components(StaticInjection, sys)
        isa(g, ElectricLoad) && continue
        g.bus.bustype == ACBusTypes.REF &&
            error("A device is already attached to the REF bus")
    end

    slack_bus = [b for b in PSY.get_components(Bus, sys) if b.bustype == ACBusTypes.REF][1]
    inf_source = Source(
        name = "InfBus", #name
        available = true, #availability
        active_power = 0.0,
        reactive_power = 0.0,
        bus = slack_bus, #bus
        R_th = 0.0,
        X_th = 1e-6, #Xth
    )
    PSY.add_component!(sys, inf_source)
    return
end

@testset "CurrentPlayback from Source" begin
    tspan = (0.0, 2.0)
    step = 1e-2
    tsteps = tspan[1]:step:tspan[2]
    sys_dir = joinpath(dirname(@__FILE__), "data_tests/OMIB.raw")
    sys = System(sys_dir, runchecks = false)

    path = (joinpath(pwd(), "test_CurrentPlayback"))
    !isdir(path) && mkdir(path)

    try
        add_source_to_ref(sys)

        #Attach dynamic generator
        gen = [g for g in get_components(Generator, sys)][1]
        case_gen = dyn_gen_second_order(gen)
        add_component!(sys, case_gen, gen)
        show_components(sys, Line, [:b])

        gen = get_component(DynamicInjection, sys, "generator-102-1") #show_components(sys, ThermalStandard)
        pert = ControlReferenceChange(1.0, gen, :V_ref, 0.5)

        #Define Simulation Problem
        sim = Simulation(
            MassMatrixModel,
            sys, # system
            path,
            tspan,
            pert,
        )

        #Solve problem
        @test execute!(sim, Rodas5(), abstol = 1e-9, reltol = 1e-9, saveat = tsteps) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)
        ir_gt = PSID.get_source_real_current_series(results, "InfBus")
        v1_gt = PSID.get_voltage_magnitude_series(results, 101)
        v2_gt = PSID.get_voltage_magnitude_series(results, 102)
        using PlotlyJS
        display(plot([PlotlyJS.scatter(x = ir_gt[1], y = ir_gt[2], name = "isource-gt")]))
        display(
            plot([
                PlotlyJS.scatter(x = v1_gt[1], y = v1_gt[2], name = "v1-gt"),
                PlotlyJS.scatter(x = v2_gt[1], y = v2_gt[2], name = "v2-gt"),
            ]),
        )

        sys_2 = deepcopy(sys)
        source = get_component(Source, sys_2, "InfBus")
        b1 = get_component(Bus, sys_2, "BUS 1")
        P = get_active_power(source)
        CP = CurrentPlayback(
            name = "cp-1",
            available = true,
            bus = b1,
            active_power = get_active_power(source),
            reactive_power = get_reactive_power(source),
            base_power = get_base_power(source),
            playback_name = "InfBus",  #"BUS 1-BUS 2-i_1",#"InfBus", 
            playback_type = :Source, #:Branch, #:Source,
            playback_result = results,
            reverse_current_polarity = false,
        )
        remove_component!(sys_2, source)
        add_component!(sys_2, CP)

        sim2 = Simulation!(MassMatrixModel, sys_2, path, tspan, pert)
        #TODO - Simulation with CurrentPlayback hangs indefinitely
        #@test execute!(sim2, Rodas5(), abstol = 1e-2, reltol = 1e-2, saveat = tsteps) == PSID.SIMULATION_FINALIZED
        #results2 = read_results(sim2)
        #v1_pb = PSID.get_voltage_magnitude_series(results2, 101)
        #v2_pb = PSID.get_voltage_magnitude_series(results2, 102)
        #using PlotlyJS
        #display(plot([PlotlyJS.scatter(x= v1_pb[1], y= v1_pb[2], name ="v1-playback"), PlotlyJS.scatter(x= v2_pb[1], y= v2_pb[2], name="v2-playback") ]))
        #display(plot([PlotlyJS.scatter(x= v2_gt[1], y= v2_gt[2], name ="v2-gt"), PlotlyJS.scatter(x= v2_pb[1], y= v2_pb[2], name="v2-playback") ]))
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
