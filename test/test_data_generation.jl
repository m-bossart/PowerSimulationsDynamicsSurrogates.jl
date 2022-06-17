
sys = System("test/data_tests/OMIB.raw")
add_source_to_ref(sys)
for g in PSY.get_components(Generator, sys)
    dyn_g = dyn_gen_second_order(g)
    add_component!(sys, dyn_g, g)
end

#Showcase two different ways to perturb the system (through a new component, through a PSID perturbation)
perturbations = [
    [("InfBus", PVS([2 * pi * 3], [(0.001, 0.01)], [2 * pi * 3], [(0.0, 0.01)]))],
    [("InfBus", VStep(0.5, 0.05)), ("InfBus", VStep(0.7, -0.05))],
]

#Define the ways to change the operating point of the system        
operating_points = [GenerationLoadScale(1.0, 1.0), GenerationLoadScale(1.1, 1.1)]

dataset = generate_train_data(
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
