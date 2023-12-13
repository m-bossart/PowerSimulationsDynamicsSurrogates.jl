@testset "9 bus system - simulation stable" begin
    sys = System("test/data_tests/nine_bus_inv_gen.json")
    sys_flow_results = solve_powerflow(PowerFlows.ACPowerFlow(), sys)["flow_results"]
    perturbations =
        Vector{Union{PowerSimulationsDynamics.Perturbation, SurrogatePerturbation}}[[
            PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0)),
        ]]
    operating_points = [
        RandomOperatingPointXiao(
            generator_voltage_range = (0.96, 1.04),
            generator_power_range = (0.0, 1.0),
            load_multiplier_range = (0.5, 1.5),
        ),
    ]
    Random.seed!(2)
    generate_data_params = GenerateDataParams(
        all_lines_dynamic = true,
        tspan = (0.0, 2.0),
        tstops = 0:0.001:2.0,
        tsave = 0:0.001:2.0,
        solver = "Rodas4",
        formulation = "MassMatrix",
    )
    dataset_aux = generate_surrogate_data(
        TerminalData,
        sys,
        sys,
        perturbations,
        operating_points,
        Dict{String, Dict{Symbol, Symbol}}(
            "Bus_7 -> Bus_8" => Dict{Symbol, Symbol}(:direction => :in, :side => :from),
        ),
        generate_data_params,
    )
    Random.seed!(2)
    dataset = generate_surrogate_data(
        TerminalData,
        sys,
        sys,
        perturbations,
        operating_points,
        Dict{String, Dict{Symbol, Symbol}}(
            "Bus_7 -> Bus_8" => Dict{Symbol, Symbol}(:direction => :in, :side => :from),
        ),
        generate_data_params,
        dataset_aux = dataset_aux,
    )
    @test dataset[1].stable == true
    @test dataset_aux[1].stable == true
end

@testset "9 bus system - simulation not built" begin
    sys = System("test/data_tests/nine_bus_inv_gen.json")
    sys_flow_results = solve_powerflow(PowerFlows.ACPowerFlow(), sys)["flow_results"]
    perturbations =
        Vector{Union{PowerSimulationsDynamics.Perturbation, SurrogatePerturbation}}[[
            PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0)),
        ]]
    operating_points = [
        RandomOperatingPointXiao(
            generator_voltage_range = (0.5, 2.0),
            generator_power_range = (0.0, 1.0),
            load_multiplier_range = (0.5, 1.5),
        ),
    ]
    Random.seed!(2)
    generate_data_params = GenerateDataParams(
        all_lines_dynamic = true,
        tspan = (0.0, 2.0),
        tstops = 0:0.001:2.0,
        tsave = 0:0.001:2.0,
        solver = "Rodas4",
        formulation = "MassMatrix",
    )
    dataset_aux = generate_surrogate_data(
        TerminalData,
        sys,
        sys,
        perturbations,
        operating_points,
        Dict{String, Dict{Symbol, Symbol}}(
            "Bus_7 -> Bus_8" => Dict{Symbol, Symbol}(:direction => :in, :side => :from),
        ),
        generate_data_params,
    )
    Random.seed!(2)
    dataset = generate_surrogate_data(
        TerminalData,
        sys,
        sys,
        perturbations,
        operating_points,
        Dict{String, Dict{Symbol, Symbol}}(
            "Bus_7 -> Bus_8" => Dict{Symbol, Symbol}(:direction => :in, :side => :from),
        ),
        generate_data_params,
        dataset_aux = dataset_aux,
    )
    @test dataset[1].built == false
    @test dataset_aux[1].built == false
end

#TODO - find a better way to make this system unstable 
#By changing the PLL parameters, takes a long time to run. 
#=  @testset "9 bus system - simulation unstable" begin
    sys = System("test/data_tests/nine_bus_inv_gen.json")
    sys_flow_results = solve_powerflow(PowerFlows.ACPowerFlow(), sys)["flow_results"]


    for gfl in get_components(DynamicInverter{AverageConverter, OuterControl{ActivePowerPI, ReactivePowerPI}, CurrentModeControl, FixedDCSource, KauraPLL, LCLFilter}, sys)
        PLL = get_freq_estimator(gfl)
        set_ω_lp!(PLL, get_ω_lp(PLL) * 1000) 
        set_kp_pll!(PLL, get_kp_pll(PLL) * 1000) 
        set_ki_pll!(PLL, get_ki_pll(PLL) * 1000) 
    end 
    perturbations =
        Vector{Union{PowerSimulationsDynamics.Perturbation, SurrogatePerturbation}}[[
            PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0)),
        ]]
    operating_points = [
            RandomOperatingPointXiao(
                generator_voltage_range = (0.96, 1.04), 
                generator_power_range = (0.0, 1.0),
                load_multiplier_range = (0.5, 1.5),
            ),
        ]
    Random.seed!(2)
    generate_data_params = GenerateDataParams(
        all_lines_dynamic = true,
        tspan = (0.0, 2.0),
        tstops = 0:0.001:2.0,
        tsave = 0:0.001:2.0,
        solver = "Rodas4",
        formulation = "MassMatrix",
    )
    dataset_aux = generate_surrogate_data(
        TerminalData,
        sys,
        sys,
        perturbations,
        operating_points,
        Dict{String, Dict{Symbol, Symbol}}(
            "Bus_7 -> Bus_8" => Dict{Symbol, Symbol}(:direction => :in, :side => :from),
        ),
        generate_data_params,
    )
    Random.seed!(2)
    dataset = generate_surrogate_data(
        TerminalData,
        sys,
        sys,
        perturbations,
        operating_points,
        Dict{String, Dict{Symbol, Symbol}}(
            "Bus_7 -> Bus_8" => Dict{Symbol, Symbol}(:direction => :in, :side => :from),
        ),
        generate_data_params,
        dataset_aux = dataset_aux
    )
    @test dataset[1].stable == true
    @test dataset_aux[1].stable == true  
end =#
