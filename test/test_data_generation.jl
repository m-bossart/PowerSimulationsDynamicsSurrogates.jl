@testset "2 bus system- generate terminal data from line" begin
    sys = System("test/data_tests/OMIB.raw")
    add_source_to_ref(sys)
    for g in PSY.get_components(Generator, sys)
        dyn_g = inv_case78(g)
        add_component!(sys, dyn_g, g)
    end
    for l in get_components(Line, sys)
        PSY.set_b!(l, (from = 0.0, to = PSY.get_b(l)[2]))
    end
    perturbations =
        Vector{Union{PowerSimulationsDynamics.Perturbation, SurrogatePerturbation}}[
            #TODO - the post processing functions for the static source cannot capture the current when a change in source voltage is implemented
            #Once this issue is solve in PSID we can confirm that the currents are the same for the VStep perturbation     
            #https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/issues/269
            #=         [
                        VStep(source_name = "InfBus", t_step = 0.5, V_step = 1.0),
                        VStep(source_name = "InfBus", t_step = 0.7, V_step = 0.95),
                    ], =#
            [
                Chirp(
                    source_name = "InfBus",
                    ω1 = 2 * pi * 3,
                    ω2 = 2 * pi * 3,
                    tstart = 0.1,
                    N = 0.5,
                    V_amp = 0.01,
                    ω_amp = 0.01,
                ),
            ],
            [
                PVS(
                    source_name = "InfBus",
                    internal_voltage_frequencies = [2 * pi * 3],
                    internal_voltage_coefficients = [(0.001, 0.01)],
                    internal_angle_frequencies = [2 * pi * 3],
                    internal_angle_coefficients = [(0.0, 0.01)],
                ),
            ],
        ]

    #Define the ways to change the operating point of the system        
    operating_points = [
        GenerationLoadScale(generation_scale = 1.0, load_scale = 1.0),
        GenerationLoadScale(generation_scale = 1.1, load_scale = 1.1),
        #  RandomOperatingPointXiao(),
    ]

    Random.seed!(2)
    dataset = generate_surrogate_data(
        TerminalData,
        sys,
        sys,
        perturbations,
        operating_points,
        Dict{String, Dict{Symbol, Symbol}}(
            "BUS 1 -> BUS 2" => Dict{Symbol, Symbol}(:direction => :in, :side => :from),
        ),
        GenerateDataParams(
            all_lines_dynamic = true,
            tstops = 0:0.001:1.0,
            tsave = 0:0.001:1.0,
            solver = "Rodas4",
            formulation = "MassMatrix",
        ),
    )
    p = generate_empty_plot(TerminalData)
    for d in dataset
        add_data_trace!(p, d)
    end
    display(p)
end

@testset "1 bus system- generate terminal data from source" begin
    sys = System("test/data_tests/OneBus.raw")
    slack_bus = [b for b in PSY.get_components(Bus, sys) if b.bustype == BusTypes.REF][1]
    inf_source = Source(
        name = "InfBus", #name
        available = true, #availability
        active_power = -0.5,
        reactive_power = -0.1,
        bus = slack_bus, #bus
        R_th = 0.0,
        X_th = 5e-6, #Xth
    )
    PSY.add_component!(sys, inf_source)
    for s in get_components(ThermalStandard, sys)
        dyngen = inv_case78(s)
        add_component!(sys, dyngen, s)
    end

    display(sys)
    #Showcase two different ways to perturb the system (through a new component, through a PSID perturbation)
    perturbations =
        Vector{Union{PowerSimulationsDynamics.Perturbation, SurrogatePerturbation}}[
        #TODO - the post processing functions for the static source cannot capture the current when a change in source voltage is implemented
        #Once this issue is solve in PSID we can confirm that the currents are the same for the VStep perturbation     
        #https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/issues/269
        #=         [
                    VStep(source_name = "InfBus", t_step = 0.5, V_step = 1.0),
                    VStep(source_name = "InfBus", t_step = 0.7, V_step = 0.95),
                ],  =#
        [
            Chirp(
                source_name = "InfBus",
                ω1 = 2 * pi * 3,
                ω2 = 2 * pi * 3,
                tstart = 0.1,
                N = 0.5,
                V_amp = 0.01,
                ω_amp = 0.01,
            ),
        ],]

    #Define the ways to change the operating point of the system        
    operating_points = [
        ScaleSource(
            source_name = "InfBus",
            V_scale = 1.0,
            θ_scale = 1.0,
            P_scale = 1.0,
            Q_scale = 1.0,
        ),
        ScaleSource(
            source_name = "InfBus",
            V_scale = 1.1,
            θ_scale = 1.0,
            P_scale = 0.9,
            Q_scale = 1.2,
        ),
    ]

    Random.seed!(2)
    dataset = generate_surrogate_data(
        TerminalData,
        sys,
        sys,
        perturbations,
        operating_points,
        Dict{String, Dict{Symbol, Symbol}}(
            "InfBus" => Dict{Symbol, Symbol}(:direction => :in),
        ),
        GenerateDataParams(
            all_lines_dynamic = true,
            tstops = 0:0.001:1.0,
            tsave = 0:0.001:1.0,
            solver = "Rodas4",
            formulation = "MassMatrix",
        ),
    )
    p = generate_empty_plot(TerminalData)
    for d in dataset
        add_data_trace!(p, d)
    end
    display(p)
end

@testset "1 bus system- generate FullSolutionData" begin
    sys = System("test/data_tests/OneBus.raw")
    slack_bus = [b for b in PSY.get_components(Bus, sys) if b.bustype == BusTypes.REF][1]
    inf_source = Source(
        name = "InfBus", #name
        available = true, #availability
        active_power = -0.5,
        reactive_power = -0.1,
        bus = slack_bus, #bus
        R_th = 0.0,
        X_th = 5e-6, #Xth
    )
    PSY.add_component!(sys, inf_source)
    for s in get_components(ThermalStandard, sys)
        dyngen = inv_case78(s)
        add_component!(sys, dyngen, s)
    end

    #Showcase two different ways to perturb the system (through a new component, through a PSID perturbation)
    perturbations =
        Vector{Union{PowerSimulationsDynamics.Perturbation, SurrogatePerturbation}}[
        #TODO - the post processing functions for the static source cannot capture the current when a change in source voltage is implemented
        #Once this issue is solve in PSID we can confirm that the currents are the same for the VStep perturbation     
        #https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/issues/269
        #=         [
                    VStep(source_name = "InfBus", t_step = 0.5, V_step = 1.0),
                    VStep(source_name = "InfBus", t_step = 0.7, V_step = 0.95),
                ],  =#
        [
            Chirp(
                source_name = "InfBus",
                ω1 = 2 * pi * 3,
                ω2 = 2 * pi * 3,
                tstart = 0.1,
                N = 0.5,
                V_amp = 0.01,
                ω_amp = 0.01,
            ),
        ],]

    #Define the ways to change the operating point of the system        
    operating_points = [
        ScaleSource(
            source_name = "InfBus",
            V_scale = 1.0,
            θ_scale = 1.0,
            P_scale = 1.0,
            Q_scale = 1.0,
        ),
        ScaleSource(
            source_name = "InfBus",
            V_scale = 1.1,
            θ_scale = 1.0,
            P_scale = 0.9,
            Q_scale = 1.2,
        ),
    ]

    Random.seed!(2)
    dataset = generate_surrogate_data(
        FullSolutionData,
        sys,
        sys,
        perturbations,
        operating_points,
        Dict{String, Dict{Symbol, Symbol}}(),
        GenerateDataParams(
            all_lines_dynamic = true,
            tstops = 0:0.001:1.0,
            tsave = 0:0.001:1.0,
            solver = "Rodas4",
            formulation = "MassMatrix",
        ),
    )
    p = generate_empty_plot(FullSolutionData)
    for d in dataset
        add_data_trace!(p, d)
    end
    display(p)
end

@testset "2 bus system- generate AllStatesData from generator" begin
    sys = System("test/data_tests/OMIB.raw")
    add_source_to_ref(sys)
    for g in PSY.get_components(Generator, sys)
        dyn_g = inv_case78(g)
        add_component!(sys, dyn_g, g)
    end
    for l in get_components(Line, sys)
        PSY.set_b!(l, (from = 0.0, to = PSY.get_b(l)[2]))
    end
    perturbations =
        Vector{Union{PowerSimulationsDynamics.Perturbation, SurrogatePerturbation}}[
        #TODO - the post processing functions for the static source cannot capture the current when a change in source voltage is implemented
        #Once this issue is solve in PSID we can confirm that the currents are the same for the VStep perturbation     
        #https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/issues/269
        #=         [
                    VStep(source_name = "InfBus", t_step = 0.5, V_step = 1.0),
                    VStep(source_name = "InfBus", t_step = 0.7, V_step = 0.95),
                ], =#
        [
            Chirp(
                source_name = "InfBus",
                ω1 = 2 * pi * 3,
                ω2 = 2 * pi * 3,
                tstart = 0.1,
                N = 0.5,
                V_amp = 0.01,
                ω_amp = 0.01,
            ),
        ],]

    #Define the ways to change the operating point of the system        
    operating_points = [
        GenerationLoadScale(generation_scale = 1.0, load_scale = 1.0),
        #GenerationLoadScale(generation_scale = 1.1, load_scale = 1.1),
        #  RandomOperatingPointXiao(),
    ]

    Random.seed!(2)
    dataset = generate_surrogate_data(
        AllStatesData,
        sys,
        sys,
        perturbations,
        operating_points,
        Dict{String, Dict{Symbol, Symbol}}("generator-102-1" => Dict{Symbol, Symbol}()),
        GenerateDataParams(
            all_lines_dynamic = true,
            tstops = 0:0.001:1.0,
            tsave = 0:0.001:1.0,
            solver = "Rodas4",
            formulation = "MassMatrix",
        ),
    )
    p = generate_empty_plot(AllStatesData)
    for d in dataset
        add_data_trace!(p, d)
    end
    display(p)
end

@testset "1 bus system- generate FullSolutionData from initial conditions" begin
    sys = System("test/data_tests/OneBus.raw")
    slack_bus = [b for b in PSY.get_components(Bus, sys) if b.bustype == BusTypes.REF][1]
    inf_source = Source(
        name = "InfBus", #name
        available = true, #availability
        active_power = -0.5,
        reactive_power = -0.1,
        bus = slack_bus, #bus
        R_th = 0.0,
        X_th = 5e-6, #Xth
    )
    PSY.add_component!(sys, inf_source)
    for s in get_components(ThermalStandard, sys)
        dyngen = inv_case78(s)
        add_component!(sys, dyngen, s)
    end

    sim = Simulation!(MassMatrixModel, sys, pwd(), (0.0, 1.0))
    ics = [PSID.get_initial_conditions(sim)]

    operating_points = [GenerationLoadScale(generation_scale = 1.0, load_scale = 1.0)]

    Random.seed!(2)
    dataset = generate_surrogate_data(
        FullSolutionData,
        sys,
        sys,
        ics,
        operating_points,
        Dict{String, Dict{Symbol, Symbol}}(),
        GenerateDataParams(
            all_lines_dynamic = true,
            tstops = 0:0.001:1.0,
            tsave = 0:0.001:1.0,
            solver = "Rodas4",
            formulation = "MassMatrix",
        ),
    )
    p = generate_empty_plot(FullSolutionData)
    for d in dataset
        add_data_trace!(p, d)
    end
    display(p)
end
