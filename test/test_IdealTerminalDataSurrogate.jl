@testset "Add IdealTerminalDataSurrogate to system" begin
    sys = System(100)
    bus = ACBus(nothing)
    set_bustype!(bus, ACBusTypes.SLACK)
    add_component!(sys, bus)
    source = Source(nothing)
    set_bus!(source, bus)
    add_component!(sys, source)
    surrogate = IdealTerminalDataSurrogate(nothing)
    add_component!(sys, surrogate, source)
    @test get_components(IdealTerminalDataSurrogate, sys) !== nothing
end

@testset "Build and Execute IdealTerminalDataSurrogate Simulation with Pref Perturbation " begin
    tspan = (0.0, 10.0)
    Random.seed!(1234)
    tfault = 5.0

    #Capacitance at connecting bus must be zero so that the current state of the line matches exactly the current injected by the surrogate
    sys_full = build_system(PSIDSystems, "psid_11bus_andes")
    l1 = get_component(Line, sys_full, "BUS 7-BUS 8-i_1")
    set_b!(l1, (from = get_b(l1).from, to = 0.0))
    l2 = get_component(Line, sys_full, "BUS 7-BUS 8-i_2")
    set_b!(l2, (from = get_b(l2).from, to = 0.0))

    dyngen = get_component(DynamicGenerator, sys_full, "generator-4-1")
    pert = PSID.ControlReferenceChange(tfault, dyngen, :P_ref, 0.04)
    sim =
        Simulation!(MassMatrixModel, sys_full, pwd(), tspan, pert; all_lines_dynamic = true)   #Lines must be dynamic to have current as dynamic state

    execute!(sim, Rodas5P(); abstol = 1e-11, saveat = 0.01)
    res_ode = read_results(sim)
    t_gt, v4_gt = get_state_series(res_ode, ("BUS 8-BUS 9-i_1", :Il_R))
    reference_solution = res_ode.solution
    gsm = PSID.make_global_state_map(sim.inputs)
    ir_ixs = [gsm["BUS 7-BUS 8-i_1"][:Il_R], gsm["BUS 7-BUS 8-i_2"][:Il_R]]
    ii_ixs = [gsm["BUS 7-BUS 8-i_1"][:Il_I], gsm["BUS 7-BUS 8-i_2"][:Il_I]]

    sys, _ = PSIDS.create_validation_system_from_buses(sys_full, [1, 2, 5, 6, 7, 8])
    for source in get_components(Source, sys)
        s = IdealTerminalDataSurrogate(
            name = get_name(source),
            τ = 0.5,
            reference_solution = reference_solution,
            reference_solution_ir_indices = ir_ixs,
            reference_solution_ii_indices = ii_ixs,
        )
        add_component!(sys, s, source)
    end
    dyngen = get_component(DynamicGenerator, sys, "generator-4-1")
    pert = PSID.ControlReferenceChange(5.0, dyngen, :P_ref, 0.04)
    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, pert; all_lines_dynamic = true)
    show_states_initial_value(sim)
    execute!(sim, MethodOfSteps(Rodas5P()); abstol = 1e-11, saveat = 0.01)
    res_dde = read_results(sim)
    t_surr, v4_surr = get_state_series(res_dde, ("BUS 8-BUS 9-i_1", :Il_R))
    println(
        "ODE solution in $(length(reference_solution)) steps, DDE solution in $(length(res_dde.solution)) steps",
    )
    @test norm(v4_surr .- v4_gt) < 9e-5

    #plt = PlotlyJS.plot()
    #add_trace!(plt, PlotlyJS.scatter(x=t_gt, y = v4_gt, name="GT- ODE"))
    #add_trace!(plt, PlotlyJS.scatter(x=t_surr, y = v4_surr, name="Ideal Surrogate"))        
    #display(plt)
end

@testset "Build and Execute IdealTerminalDataSurrogate Simulation with Pref Perturbation " begin
    tspan = (0.0, 10.0)
    Random.seed!(1234)
    tfault = 5.0
    fault_branch = "BUS 8-BUS 9-i_2"

    #Capacitance at connecting bus must be zero so that the current state of the line matches exactly the current injected by the surrogate
    sys_full = build_system(PSIDSystems, "psid_11bus_andes")
    l1 = get_component(Line, sys_full, "BUS 7-BUS 8-i_1")
    set_b!(l1, (from = get_b(l1).from, to = 0.0))
    l2 = get_component(Line, sys_full, "BUS 7-BUS 8-i_2")
    set_b!(l2, (from = get_b(l2).from, to = 0.0))
    for l in get_components(Line, sys_full)  #Lines must be dynamic to have current as dynamic state BUT tripped line can't be dynamic
        if get_name(l) !== fault_branch
            dyn_branch = DynamicBranch(l)
            add_component!(sys_full, dyn_branch)
        end
    end
    pert = PSID.BranchTrip(tfault, Line, fault_branch)
    sim = Simulation!(MassMatrixModel, sys_full, pwd(), tspan, pert)

    execute!(sim, Rodas5P(); abstol = 1e-11, saveat = 0.01)
    res_ode = read_results(sim)
    t_gt, v4_gt = get_state_series(res_ode, ("BUS 8-BUS 9-i_1", :Il_R))
    reference_solution = res_ode.solution
    gsm = PSID.make_global_state_map(sim.inputs)
    ir_ixs = [gsm["BUS 7-BUS 8-i_1"][:Il_R], gsm["BUS 7-BUS 8-i_2"][:Il_R]]
    ii_ixs = [gsm["BUS 7-BUS 8-i_1"][:Il_I], gsm["BUS 7-BUS 8-i_2"][:Il_I]]

    sys, _ = PSIDS.create_validation_system_from_buses(sys_full, [1, 2, 5, 6, 7, 8])
    for source in get_components(Source, sys)
        s = IdealTerminalDataSurrogate(
            name = get_name(source),
            τ = 0.5,
            reference_solution = reference_solution,
            reference_solution_ir_indices = ir_ixs,
            reference_solution_ii_indices = ii_ixs,
        )
        add_component!(sys, s, source)
    end
    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, pert)
    show_states_initial_value(sim)
    execute!(sim, MethodOfSteps(Rodas5P()); abstol = 1e-11, saveat = 0.01)
    res_dde = read_results(sim)
    t_surr, v4_surr = get_state_series(res_dde, ("BUS 8-BUS 9-i_1", :Il_R))
    println(
        "ODE solution in $(length(reference_solution)) steps, DDE solution in $(length(res_dde.solution)) steps",
    )
    @test norm(v4_surr .- v4_gt) < 0.3

    #plt = PlotlyJS.plot()
    #add_trace!(plt, PlotlyJS.scatter(x=t_gt, y = v4_gt, name="GT- ODE"))
    #add_trace!(plt, PlotlyJS.scatter(x=t_surr, y = v4_surr, name="Ideal Surrogate"))       
    #display(plt)
end

@testset "Test with model evaluation" begin
    sim_times = []
    for n in [10, 100000]
        tspan = (0.0, 10.0)
        Random.seed!(1234)
        tfault = 5.0

        #Capacitance at connecting bus must be zero so that the current state of the line matches exactly the current injected by the surrogate
        sys_full = build_system(PSIDSystems, "psid_11bus_andes")
        l1 = get_component(Line, sys_full, "BUS 7-BUS 8-i_1")
        set_b!(l1, (from = get_b(l1).from, to = 0.0))
        l2 = get_component(Line, sys_full, "BUS 7-BUS 8-i_2")
        set_b!(l2, (from = get_b(l2).from, to = 0.0))

        dyngen = get_component(DynamicGenerator, sys_full, "generator-4-1")
        pert = PSID.ControlReferenceChange(tfault, dyngen, :P_ref, 0.04)
        sim = Simulation!(
            MassMatrixModel,
            sys_full,
            pwd(),
            tspan,
            pert;
            all_lines_dynamic = true,
        )   #Lines must be dynamic to have current as dynamic state

        execute!(sim, Rodas5P(); abstol = 1e-11, saveat = 0.01)
        res_ode = read_results(sim)
        t_gt, v4_gt = get_state_series(res_ode, ("BUS 8-BUS 9-i_1", :Il_R))
        reference_solution = res_ode.solution
        gsm = PSID.make_global_state_map(sim.inputs)
        ir_ixs = [gsm["BUS 7-BUS 8-i_1"][:Il_R], gsm["BUS 7-BUS 8-i_2"][:Il_R]]
        ii_ixs = [gsm["BUS 7-BUS 8-i_1"][:Il_I], gsm["BUS 7-BUS 8-i_2"][:Il_I]]
        model = Lux.Chain(Lux.Dense(2, n), Lux.Dense(n, 2))
        rng = Random.default_rng()
        Random.seed!(rng, 0)
        ps, st = Lux.setup(rng, model)
        sys, _ = PSIDS.create_validation_system_from_buses(sys_full, [1, 2, 5, 6, 7, 8])
        for source in get_components(Source, sys)
            s = IdealTerminalDataSurrogate(
                name = get_name(source),
                τ = 0.5,
                reference_solution = reference_solution,
                reference_solution_ir_indices = ir_ixs,
                reference_solution_ii_indices = ii_ixs,
                ext = Dict{String, Any}("model" => model, "ps" => ps, "st" => st),
            )
            add_component!(sys, s, source)
        end
        dyngen = get_component(DynamicGenerator, sys, "generator-4-1")
        pert = PSID.ControlReferenceChange(5.0, dyngen, :P_ref, 0.04)
        sim =
            Simulation!(MassMatrixModel, sys, pwd(), tspan, pert; all_lines_dynamic = true)
        show_states_initial_value(sim)
        execute!(sim, MethodOfSteps(Rodas5P()); abstol = 1e-11, saveat = 0.01)
        res_dde = read_results(sim)
        println("Solution time: $(res_dde.time_log[:timed_solve_time])")
        t_surr, v4_surr = get_state_series(res_dde, ("BUS 8-BUS 9-i_1", :Il_R))
        println(
            "ODE solution in $(length(reference_solution)) steps, DDE solution in $(length(res_dde.solution)) steps",
        )
        @test norm(v4_surr .- v4_gt) < 9e-5
        push!(sim_times, res_dde.time_log[:timed_solve_time])
    end
    @test sim_times[2] / sim_times[1] > 5.0
    #plt = PlotlyJS.plot()
    #add_trace!(plt, PlotlyJS.scatter(x=t_gt, y = v4_gt, name="GT- ODE"))
    #add_trace!(plt, PlotlyJS.scatter(x=t_surr, y = v4_surr, name="Ideal Surrogate"))        
    #display(plt)
end
