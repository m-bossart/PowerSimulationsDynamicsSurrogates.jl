include(joinpath(dirname(@__FILE__), "data_tests/dynamic_test_data.jl"))

@testset "PSY Frequency Source 1" begin
    sys = PSY.System(100)
    bus = PSY.Bus(nothing)
    PSY.set_bustype!(bus, BusTypes.SLACK)
    PSY.add_component!(sys, bus)
    source = PSY.Source(nothing)
    PSY.set_bus!(source, bus)
    PSY.add_component!(sys, source)
    cvs = FrequencySource(nothing)
    PSY.add_component!(sys, cvs, source)
    @test get_components(FrequencySource, sys).length !== 0
end

@testset "PSY Frequency Source 2" begin
    sys = PSY.System(100)
    bus = PSY.Bus(nothing)
    PSY.set_bustype!(bus, BusTypes.REF)
    PSY.add_component!(sys, bus)
    source = PSY.Source(nothing)
    PSY.set_bus!(source, bus)
    PSY.add_component!(sys, source)
    cvs = FrequencySource(nothing)
    PSY.add_component!(sys, cvs, source)
    @test get_components(FrequencySource, sys).length !== 0
    # sys2, result = PSY.validate_serialization(sys)
    # @test result
end

function cvs_simple(source)
    return FrequencySource(
        name = PSY.get_name(source),
        R_th = PSY.get_R_th(source),
        X_th = PSY.get_X_th(source),
        V_ref = 1.0,
        ω_ref = 1.0,
        Tv = 1.0,
        Tω = 0.5,
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

#@testset "Frequency Source PSID ResidualModel" begin
tspan = (0.0, 15.0)
step = 1e-2
tsteps = tspan[1]:step:tspan[2]

sys_dir = joinpath(dirname(@__FILE__), "data_tests/OMIB.raw")
sys = System(sys_dir, runchecks = false)

path = (joinpath(pwd(), "test_FrequencySource"))
!isdir(path) && mkdir(path)

#try
add_source_to_ref(sys)

#Attach dynamic generator
gen = [g for g in get_components(Generator, sys)][1]
case_gen = dyn_gen_second_order(gen)
add_component!(sys, case_gen, gen)

#Attach periodic variable source
source = [s for s in get_components(Source, sys)][1]
cvs = cvs_simple(source)
add_component!(sys, cvs, source)

#s_device = get_component(Source, sys, "InfBus")
#s_change = SourceBusVoltageChange(1.0, s_device, :V_ref, 1.02)

fs_device = get_component(FrequencySource, sys, "InfBus")
fs_change = ControlReferenceChange(2.0, fs_device, :ω_ref, 1.01)
#Define Simulation Problem
sim = Simulation!(
    ResidualModel,
    sys, # system
    path,
    tspan,
    fs_change,
)
x0_init = read_initial_conditions(sim)
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

# Obtain data for source
Vt_source = get_state_series(results, ("InfBus", :Vt))
θt_source = get_state_series(results, ("InfBus", :θt))
ω_source = get_state_series(results, ("InfBus", :ω))

# Obtain data for get
Vt_gen = get_voltage_magnitude_series(results, 102)
θt_gen = get_voltage_angle_series(results, 102)
ω_gen = get_state_series(results, ("generator-102-1", :ω))
δ_gen = get_state_series(results, ("generator-102-1", :δ))
P_gen = get_activepower_series(results, "generator-102-1")
p4 = plot(P_gen, label = "P-gen")
p1 = plot(Vt_gen, label = "Vt-gen")
p2 = plot(θt_gen, label = "θt-gen")
p3 = plot(ω_gen, label = "ω-gen")

plot!(p1, Vt_source, label = "Vt-source")
plot!(p2, θt_source, label = "θt-source")
plot!(p3, ω_source, label = "ω-source")
plot!(p2, δ_gen, label = "δ-gen")

display(plot(p1, p2, p3, p4, layout = (2, 2)))
#finally
@info("removing test files")
rm(path, force = true, recursive = true)
#end
#end
