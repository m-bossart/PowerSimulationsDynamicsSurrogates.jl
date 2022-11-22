@testset "2 bus system- generate data from line" begin
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
        sys,
        sys,
        perturbations,
        operating_points,
        SteadyStateNODEDataParams(
            location_of_data_collection = [("BUS 1-BUS 2-i_1", :from)],
        ),
        GenerateDataParams(
            all_lines_dynamic = true,
            tstops = 0:0.001:1.0,
            tsave = 0:0.001:1.0,
            solver = "Rodas4",
            formulation = "MassMatrix",
        ),
    )

    p1 = plot()
    p2 = plot()
    p3 = plot()
    p4 = plot()
    for d in dataset
        display(length(d.tsteps))
        plot!(p1, d.tsteps, d.surrogate_real_voltage[1, :], ylabel = "real voltage")
        plot!(p2, d.tsteps, d.surrogate_imag_voltage[1, :], ylabel = "imag voltage")
        plot!(p3, d.tsteps, d.real_current[1, :], ylabel = "real current")
        plot!(p4, d.tsteps, d.imag_current[1, :], ylabel = "imag current")
    end
    display(plot(p1, p2, p3, p4, title = "data from branch"))
end

@testset "2 bus system- generate data from source" begin
    sys = System("test/data_tests/OMIB.raw")
    add_source_to_ref(sys)
    for g in PSY.get_components(Generator, sys)
        dyn_g = inv_case78(g)
        add_component!(sys, dyn_g, g)
    end
    for l in get_components(Line, sys)
        PSY.set_b!(l, (from = 0.0, to = PSY.get_b(l)[2]))
    end
    for b in get_components(Line, sys)
        @warn b
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
        sys,
        sys,
        perturbations,
        operating_points,
        SteadyStateNODEDataParams(location_of_data_collection = [("InfBus", :source)]),
        GenerateDataParams(
            all_lines_dynamic = true,
            tstops = 0:0.001:1.0,
            tsave = 0:0.001:1.0,
            solver = "Rodas4",
            formulation = "MassMatrix",
        ),
    )

    p1 = plot()
    p2 = plot()
    p3 = plot()
    p4 = plot()
    for d in dataset
        plot!(p1, d.tsteps, d.surrogate_real_voltage[1, :], ylabel = "real voltage")
        plot!(p2, d.tsteps, d.surrogate_imag_voltage[1, :], ylabel = "imag voltage")
        plot!(p3, d.tsteps, d.real_current[1, :], ylabel = "real current")
        plot!(p4, d.tsteps, d.imag_current[1, :], ylabel = "imag current")
    end
    display(plot(p1, p2, p3, p4, title = "data from source"))
end
