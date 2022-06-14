include(joinpath(dirname(@__FILE__), "data/dynamic_test_data.jl"))

@testset "PSY Chirp Source 1" begin
    sys = PSY.System(100)
    bus = PSY.Bus(nothing)
    PSY.set_bustype!(bus, BusTypes.SLACK)
    PSY.add_component!(sys, bus)
    source = PSY.Source(nothing)
    PSY.set_bus!(source, bus)
    PSY.add_component!(sys, source)
    cvs = ChirpVariableSource(nothing)
    PSY.add_component!(sys, cvs, source)
    @test get_components(ChirpVariableSource, sys).length !== 0
end

@testset "PSY Chirp Source 2" begin
    sys = PSY.System(100)
    bus = PSY.Bus(nothing)
    PSY.set_bustype!(bus, BusTypes.REF)
    PSY.add_component!(sys, bus)
    source = PSY.Source(nothing)
    PSY.set_bus!(source, bus)
    PSY.add_component!(sys, source)
    cvs = ChirpVariableSource(nothing)
    PSY.add_component!(sys, cvs, source)
    @test get_components(ChirpVariableSource, sys).length !== 0
    # sys2, result = PSY.validate_serialization(sys)
    # @test result
end

function cvs_simple(source)
    return ChirpVariableSource(
        name = PSY.get_name(source),
        R_th = PSY.get_R_th(source),
        X_th = PSY.get_X_th(source),
        ω1=2.0*pi,
        ω2=2.0*pi*3,
        tstart=1.0,
        N=4.0,
        V_amp=0.1,
        θ_amp=0.01,
    )
end

function add_source_to_ref(sys::PSY.System)
    for g in PSY.get_components(StaticInjection, sys)
        isa(g, ElectricLoad) && continue
        g.bus.bustype == BusTypes.REF &&
            error("A device is already attached to the REF bus")
    end

    slack_bus = [b for b in PSY.get_components(Bus, sys) if b.bustype == BusTypes.REF][1]
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

@testset "Chirp PSID ResidualModel" begin
    tspan = (0.0, 1.0);
    step = 1e-1
    tsteps = tspan[1]:step:tspan[2]

    sys_dir = joinpath(dirname(@__FILE__), "data/OMIB.raw")
    sys = System(sys_dir, runchecks = false)

    
    path = (joinpath(pwd(), "test_ChirpSource"))
    !isdir(path) && mkdir(path)

    try
        add_source_to_ref(sys)

        #Attach dynamic generator
        gen = [g for g in get_components(Generator, sys)][1]
        case_gen = dyn_gen_second_order(gen)
        add_component!(sys, case_gen, gen)

        #Attach periodic variable source
        source = [s for s in get_components(Source, sys)][1]
        cvs = cvs_simple(source)
        add_component!(sys, cvs, source)

        #Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            sys, # system
            path,
            tspan,
        )

        # Test Initial Condition
        # diff_val = [0.0]
        # res = PSID.get_init_values_for_comparison(sim)
        # for (k, v) in test09_x0_init
        #     diff_val[1] += LinearAlgebra.norm(res[k] - v)
        # end
        # @test (diff_val[1] < 1e-3)

        #Solve problem
        @test execute!(sim, IDA(), saveat = tsteps) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for angles
        Vt_sim = get_state_series(results, ("InfBus", :Vt))
        θt_sim = get_state_series(results, ("InfBus", :θt))

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

# @testset "Chirp PSID  MassMatrixModel" begin
#     tspan = (0.0, 1.0);
#     step = 1e-1
#     tsteps = tspan[1]:step:tspan[2]

#     sys_dir = joinpath(dirname(@__FILE__), "data/OMIB.raw")
#     sys = System(sys_dir, runchecks = false)

    
#     path = (joinpath(pwd(), "test_ChirpSource"))
#     !isdir(path) && mkdir(path)
#     try

#         add_source_to_ref(sys)

#         #Attach dynamic generator
#         gen = [g for g in get_components(Generator, sys)][1]
#         case_gen = dyn_gen_second_order(gen)
#         add_component!(sys, case_gen, gen)

#         #Attach periodic variable source
#         source = [s for s in get_components(Source, sys)][1]
#         cvs = cvs_simple(source)
#         add_component!(sys, cvs, source)

#         #Define Simulation Problem
#         sim = Simulation!(
#             MassMatrixModel,
#             sys, # system
#             path,
#             tspan,
#         )

#         #Solve problem
#         @test execute!(sim, Rodas4(), saveat = tsteps) == PSID.SIMULATION_FINALIZED
#         results = read_results(sim)

#         # Obtain data for angles
#         Vt_sim = get_state_series(results, ("InfBus", :Vt))
#         θt_sim = get_state_series(results, ("InfBus", :θt))

#     finally
#         @info("removing test files")
#         rm(path, force = true, recursive = true)
#     end
# end