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
        X_th = 5e-6, #Xth
    )
    PSY.add_component!(sys, inf_source)
    return
end

@testset "Test SourceLoad" begin
    tspan = (0.0, 6.0)
    step = 1e-2
    tsteps = tspan[1]:step:tspan[2]
    sys_dir = joinpath(dirname(@__FILE__), "data_tests/OMIB.raw")
    sys = System(sys_dir, runchecks = false)
    path = (joinpath(pwd(), "test_SourceLoad"))
    !isdir(path) && mkdir(path)

    try
        #Build a base system with a Generator at reference bus and source at PV bus.
        add_source_to_ref(sys)
        gen = [g for g in get_components(Generator, sys)][1]
        case_gen = dyn_gen_second_order(gen)
        add_component!(sys, case_gen, gen)
        b1 = get_component(Bus, sys, "BUS 1")
        b2 = get_component(Bus, sys, "BUS 2")
        set_bustype!(b1, ACBusTypes.PV)
        set_bustype!(b2, ACBusTypes.REF)
        solve_ac_powerflow!(sys)
        sys_load = deepcopy(sys)
        sys_sourceload = deepcopy(sys)
        sys_source = deepcopy(sys)
        pert = PSID.BranchImpedanceChange(0.1, Line, "BUS 1-BUS 2-i_1", 1.5)

        #SOURCE SYSTEM
        sim_source = Simulation!(ResidualModel, sys_source, path, (0.0, 0.15), pert)
        @test execute!(sim_source, IDA(), saveat = tsteps) == PSID.SIMULATION_FINALIZED
        results_source = read_results(sim_source)
        Vt_source = get_voltage_magnitude_series(results_source, 102)

        #SOURCELOAD SYSTEM
        source = get_component(Source, sys_sourceload, "InfBus")
        b1 = get_component(Bus, sys_sourceload, "BUS 1")
        P = get_active_power(source)
        sourceload = SourceLoad(
            name = "load-1",
            available = true,
            bus = b1,
            base_power = get_base_power(source),
            constant_active_power = 0.0,
            constant_reactive_power = 0.0,
            impedance_active_power = P,
            impedance_reactive_power = 0.0,
            current_active_power = 0.0,
            current_reactive_power = 0.0,
            max_constant_active_power = 1.0,
            max_constant_reactive_power = 1.0,
            max_impedance_active_power = 1.0,
            max_impedance_reactive_power = 1.0,
            max_current_active_power = 1.0,
            max_current_reactive_power = 1.0,
        )
        remove_component!(sys_sourceload, source)
        add_component!(sys_sourceload, sourceload)

        sim_sourceload = Simulation!(ResidualModel, sys_sourceload, path, (0.0, 0.15), pert)
        @test execute!(sim_sourceload, IDA(), saveat = tsteps) == PSID.SIMULATION_FINALIZED
        results_sourceload = read_results(sim_sourceload)
        Vt_sourceload = get_voltage_magnitude_series(results_sourceload, 102)

        #LOAD SYSTEM 
        b1 = get_component(Bus, sys_load, "BUS 1")
        set_bustype!(b1, ACBusTypes.PQ)
        source = get_component(Source, sys_load, "InfBus")
        P = get_active_power(source)
        Q = get_reactive_power(source)
        load = StandardLoad(
            name = "load-1",
            available = true,
            bus = b1,
            base_power = get_base_power(source),
            constant_active_power = 0.0,
            constant_reactive_power = 0.0,
            impedance_active_power = -P,
            impedance_reactive_power = -Q,
            current_active_power = 0.0,
            current_reactive_power = 0.0,
            max_constant_active_power = 1.0,
            max_constant_reactive_power = 1.0,
            max_impedance_active_power = 1.0,
            max_impedance_reactive_power = 1.0,
            max_current_active_power = 1.0,
            max_current_reactive_power = 1.0,
        )
        remove_component!(sys_load, source)
        add_component!(sys_load, load)
        sim_load = Simulation!(ResidualModel, sys_load, path, (0.0, 0.15), pert)
        @test execute!(sim_load, IDA(), saveat = tsteps) == PSID.SIMULATION_FINALIZED
        results_load = read_results(sim_load)
        Vt_load = get_voltage_magnitude_series(results_load, 102)

        #Test that Load system response matches SourceLoad system response
        @test LinearAlgebra.norm(Vt_load[2] .- Vt_sourceload[2], Inf) <= 1e-10
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
