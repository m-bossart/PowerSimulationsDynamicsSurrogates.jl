#@testset "2 bus system" begin
sys = System("test/data_tests/OMIB.raw")
add_source_to_ref(sys)
for g in PSY.get_components(Generator, sys)
    dyn_g = inv_case78(g)
    add_component!(sys, dyn_g, g)
end
for l in get_components(Line, sys)
    PSY.set_b!(l, (from = 0.0, to = PSY.get_b(l)[2]))
end

#Showcase two different ways to perturb the system (through a new component, through a PSID perturbation)
perturbations = Vector{Union{PowerSimulationsDynamics.Perturbation, SurrogatePerturbation}}[
    [
        VStep(source_name = "InfBus", t_step = 0.5, V_step = 1.0),
        VStep(source_name = "InfBus", t_step = 0.7, V_step = 0.95),
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
    SteadyStateNODEDataParams(connecting_branch_names = [("BUS 1-BUS 2-i_1", :from)]),
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
    plot!(p1, d.tsteps, d.surrogate_real_voltage[1, :])
    plot!(p2, d.tsteps, d.opposite_real_voltage[1, :])
    plot!(p3, d.tsteps, d.branch_real_current[1, :])
    plot!(p4, d.tsteps, d.branch_imag_current[1, :])
end
display(plot(p1, p2, p3, p4))
#end

#TODO - add a test for a larger (stable) system
