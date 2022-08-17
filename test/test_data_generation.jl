@testset "2 bus system" begin
    sys = System("test/data_tests/OMIB.raw")
    add_source_to_ref(sys)
    for g in PSY.get_components(Generator, sys)
        dyn_g = dyn_gen_second_order(g)
        add_component!(sys, dyn_g, g)
    end

    #Showcase two different ways to perturb the system (through a new component, through a PSID perturbation)
    perturbations =
        Vector{Union{PowerSimulationsDynamics.Perturbation, SurrogatePerturbation}}[
            [
                PVS(
                    source_name = "InfBus",
                    internal_voltage_frequencies = [2 * pi * 3],
                    internal_voltage_coefficients = [(0.001, 0.01)],
                    internal_angle_frequencies = [2 * pi * 3],
                    internal_angle_coefficients = [(0.0, 0.01)],
                ),
            ],
            [
                VStep(source_name = "InfBus", t_step = 0.5, ΔV = 0.05),
                VStep(source_name = "InfBus", t_step = 0.7, ΔV = -0.05),
            ],
        ]

    #Define the ways to change the operating point of the system        
    operating_points = [
        GenerationLoadScale(generation_scale = 1.0, load_scale = 1.0),
        GenerationLoadScale(generation_scale = 1.1, load_scale = 1.1),
        RandomOperatingPointXiao(),
    ]

    Random.seed!(2)
    dataset = generate_surrogate_data(
        sys,
        sys,
        perturbations,
        operating_points,
        SteadyStateNODEDataParams(connecting_branch_names = [("BUS 1-BUS 2-i_1", :from)]),
        GenerateDataParams(),
    )

    p = plot()
    for d in dataset
        plot!(p, d.tsteps, d.groundtruth_current[1, :])
    end
    display(p)
end

@info "9 bus system test is unstable; ignore simulation error message"
@testset "9 bus system" begin
    sys = System("test/data_tests/9BusSystem.json")

    perturbations =
        Vector{Union{PowerSimulationsDynamics.Perturbation, SurrogatePerturbation}}[[
            RandomLoadChange(time = 0.2, load_multiplier_range = (0.8, 1.0)),
        ]]
    operating_points = [
        GenerationLoadScale(generation_scale = 1.0, load_scale = 1.0),
        RandomOperatingPointXiao(),
        RandomOperatingPointXiao(),
    ]
    Random.seed!(100)
    dataset = generate_surrogate_data(
        sys,
        sys,
        perturbations,
        operating_points,
        SteadyStateNODEDataParams(connecting_branch_names = [("Bus 2-Bus 7-i_8", :from)]),
        GenerateDataParams(solver_tols = (1e-2, 1e-2)),
    )
    p = plot()
    for d in dataset
        if d.stable == true
            plot!(p, d.tsteps, d.groundtruth_current[1, :])
        end
    end
    @test [d.stable for d in dataset] == [0, 0, 0]
    display(p)
end
