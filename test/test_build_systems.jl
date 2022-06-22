#TODO - test a system with two connecting buses, make sure the power flow matches exactly the original power flow 

sys = System("test/data_tests/OMIB.raw")
add_source_to_ref(sys)
display(sys)
#run_powerflow(sys)["bus_results"]
sys_train = create_subsystem_from_buses(sys, [102])
display(sys_train)
run_powerflow(sys_train)["flow_results"]
#= #simulation!(sys)

sim = Simulation!(
    ResidualModel,
    sys,
    pwd(),
    (0.0,1.0),
    #SourceBusVoltageChange(0.5, source_IB, :V_ref, 1.04),
) =#