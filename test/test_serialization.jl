
rng = Random.default_rng()
@testset "Test serialization with surrogates in the system" begin
    path = (joinpath(pwd(), "test_serialization_with_surrogates"))
    !isdir(path) && mkdir(path)
    try
        input = [0.2; 0.4]
        model = Lux.Chain(
            WrappedFunction(x -> x .* sin(1.5)),
            Lux.Dense(2, 5, NNlib.relu),
            Lux.Dense(5, 2),
        )
        Random.seed!(rng, 0)
        p1, st1 = Lux.setup(rng, model)
        sys = System("test/data_tests/9BusSystem.json")

        #Replace all dynamic models with a surrogate. 
        for g in get_components(ThermalStandard, sys)
            dyn_gen = remove_component!(sys, get_dynamic_injector(g))
            s = TerminalDataSurrogate(
                name = get_name(g),
                τ = 0.1,
                window_size = 5,
                ext = Dict{String, Any}(
                    "model_path" => "",
                    "model" => model,
                    "ps" => p1,
                    "st" => st1,
                ),
            )
            add_component!(sys, s, g)
        end

        to_json_with_surrogates(sys, joinpath(path, "test.json"))
        sys_2 = deserialize_with_surrogates(joinpath(path, "test.json"))
        g_2 = collect(get_components(ThermalStandard, sys_2))[1]

        model2 = get_ext(get_dynamic_injector(g_2))["model"]
        p2 = get_ext(get_dynamic_injector(g_2))["ps"]
        st2 = get_ext(get_dynamic_injector(g_2))["st"]

        y_1, st_1 = Lux.apply(model, input, p1, st1)
        y2, st_2 = Lux.apply(model2, input, p2, st2)
        @test isapprox(y_1, y2)
        @test length(readdir(joinpath(path, "surrogate_models"))) == 3
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
