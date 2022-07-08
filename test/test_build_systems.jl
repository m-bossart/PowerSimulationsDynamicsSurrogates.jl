@testset "9 bus system" begin
    sys_full = System("test/data_tests/9BusSystem.json")
    surrogate_buses = [1, 4, 5]#[1,4,5]

    non_surrogate_buses =
        get_number.(get_components(Bus, sys_full, x -> get_number(x) ∉ surrogate_buses))
    sys_full_flow = run_powerflow(sys_full)["flow_results"]

    sys_train, connected_branches_names_1 = create_subsystem_from_buses(sys_full, surrogate_buses)
    sys_train_flow = run_powerflow(sys_train)["flow_results"]
    display(sys_train_flow)
    display(run_powerflow(sys_train)["bus_results"])
    display(sys_train)
    sys_test, connected_branches_names_2 = create_subsystem_from_buses(sys_full, non_surrogate_buses)
    sys_test_flow = run_powerflow(sys_test)["flow_results"]
    display(sys_test_flow)
    for (ix,b) in enumerate(connected_branches_names_1)
        @test b[1] == connected_branches_names_2[ix][1]
    end 
   
    for sys_full_row in eachrow(sys_full_flow)
        for sys_train_row in eachrow(sys_train_flow)
            if sys_full_row.line_name == sys_train_row.line_name
                for (i, c) in enumerate(sys_full_row)
                    if i > 1
                        @test isapprox(c, sys_train_row[i], atol = 1e-4)
                    end
                end
            end
        end
        for sys_test_row in eachrow(sys_test_flow)
            if sys_full_row.line_name == sys_test_row.line_name
                for (i, c) in enumerate(sys_full_row)
                    if i > 1
                        @test isapprox(c, sys_test_row[i], atol = 1e-4)
                    end
                end
            end
        end
    end
    @test create_subsystem_from_buses(sys_full, [1, 6]) == false
    @test create_subsystem_from_buses(sys_full, [1, 9]) == false
end
