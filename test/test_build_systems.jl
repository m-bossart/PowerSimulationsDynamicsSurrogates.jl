@testset "9 bus system" begin
    sys_full = System("test/data_tests/9BusSystem.json")
    surrogate_buses = [2]

    sys_full_flow = run_powerflow(sys_full)["flow_results"]
    sys_train, connected_branches_names_1 =
        create_train_system_from_buses(sys_full, surrogate_buses)
    sys_train_flow = run_powerflow(sys_train)["flow_results"]
    sys_test, connected_branches_names_2 =
        create_validation_system_from_buses(sys_full, surrogate_buses)
    sys_test_flow = run_powerflow(sys_test)["flow_results"]

    display(sys_full_flow)
    display(sys_train_flow)
    display(sys_test_flow)

    for (ix, b) in enumerate(connected_branches_names_1)
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
    @warn "Error below is expected"
    @test create_validation_system_from_buses(sys_full, [2, 3, 4, 5, 7, 8, 9]) == false
    @warn "Error below is expected"
    @test create_validation_system_from_buses(sys_full, [2, 3, 4, 5, 6, 7, 8]) == false
end

@testset "3 bus line system" begin
    sys_full = System("test/data_tests/3busline.json")
    surrogate_buses = [102, 103]

    sys_full_flow = run_powerflow(sys_full)["flow_results"]
    sys_train, connected_branches_names_1 =
        create_train_system_from_buses(sys_full, surrogate_buses)
    sys_train_flow = run_powerflow(sys_train)["flow_results"]
    sys_test, connected_branches_names_2 =
        create_validation_system_from_buses(sys_full, surrogate_buses)
    sys_test_flow = run_powerflow(sys_test)["flow_results"]

    display(sys_full_flow)
    display(sys_train_flow)
    display(sys_test_flow)

    for (ix, b) in enumerate(connected_branches_names_1)
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
    @warn "Error below is expected"
    @test create_validation_system_from_buses(sys_full, [102]) == false
end

@testset "3 bus line system" begin
    sys_full = System("test/data_tests/3busline.json")
    surrogate_buses = [103]

    sys_full_flow = run_powerflow(sys_full)["flow_results"]
    sys_train, connected_branches_names_1 =
        create_train_system_from_buses(sys_full, surrogate_buses)
    sys_train_flow = run_powerflow(sys_train)["flow_results"]
    sys_test, connected_branches_names_2 =
        create_validation_system_from_buses(sys_full, surrogate_buses)
    sys_test_flow = run_powerflow(sys_test)["flow_results"]

    display(sys_full_flow)
    display(sys_train_flow)
    display(sys_test_flow)

    for (ix, b) in enumerate(connected_branches_names_1)
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
    @warn "Error below is expected"
    @test create_validation_system_from_buses(sys_full, [102]) == false
end

#TODO - Can we compare the df with bus_results as well? 
#TODO - Add a test for a system where the train system has a generator and source at the same bus (requires PR in PowerFlows.jl to distribute P, Q between Source and generator)
#TODO - Add a test for a system where the train system is a single bus --> try manual initialization? 
