# To generate train data using a reduced system requires a source at the surrogate bus to perturb. 
# More work needed for this case, skip for now. 



#= sys = System("test/data_tests/TwoBusone_gen.raw")
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
for g in PSY.get_components(Generator, sys)
    dyn_g = dyn_gen_second_order(g)       #
    add_component!(sys, dyn_g, g)
end
solve_powerflow(sys)["bus_results"]
##
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

 =#
