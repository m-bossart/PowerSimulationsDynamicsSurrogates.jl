
sys = System("test/data_tests/OMIB.raw")
add_source_to_ref(sys)
for g in PSY.get_components(Generator, sys)
    dyn_g = dyn_gen_second_order(g)
    add_component!(sys, dyn_g, g)
end

#Showcase two different ways to perturb the system (through a new component, through a PSID perturbation)
perturbations = [
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
]

dataset = generate_surrogate_data(
    sys,
    perturbations,
    operating_points,
    "SteadyStateNODEData",
    GenerateDataParams(),
)

p = plot()
for d in dataset
    plot!(p, d.tsteps, d.groundtruth_current[1, :])
end
display(p)
