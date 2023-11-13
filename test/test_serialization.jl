
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
        p1, st1 = Lux.setup(rng, model)
        sys = System("test/data_tests/9BusSystem.json")
        g = collect(get_components(ThermalStandard, sys))[1]
        set_ext!(
            g,
            Dict{String, Any}(
                "model_path" => "",
                "model" => model,
                "ps" => p1,
                "st" => st1,
            ),
        )  #Store the Flux model in the ext field of a component

        #TODO - make these functions general and move to the package, then test them here. 
        # Might require some modification into type structure of surrogate components (static or dynamic injectors?).
        function to_json_with_surrogates(sys, full_path)
            g = collect(get_components(ThermalStandard, sys))[1]    #Note: loop over abstract surrogate type to search for models.
            model = get_ext(g)["model"]
            p = get_ext(g)["ps"]
            st = get_ext(g)["st"]
            dir = dirname(full_path)
            println(dir)
            mkpath(joinpath(dir, "surrogate_models"))
            @save joinpath(dir, "surrogate_models", get_name(g)) model p st
            set_ext!(
                g,
                Dict{String, Any}(
                    "model_path" => joinpath(dir, "surrogate_models", get_name(g)),
                    "model" => nothing,
                ),
            )
            to_json(sys, full_path; force = true)
        end

        function deserialize_with_surrogates(full_path)
            sys = System(full_path)
            g = collect(get_components(ThermalStandard, sys))[1]
            @load get_ext(g)["model_path"] model p st
            set_ext!(
                g,
                Dict{String, Any}(
                    "model_path" => nothing,
                    "model" => model,
                    "ps" => p,
                    "st" => st,
                ),
            )
            return sys
        end

        to_json_with_surrogates(sys, joinpath(path, "test.json"))
        sys_2 = deserialize_with_surrogates(joinpath(path, "test.json"))
        g_2 = collect(get_components(ThermalStandard, sys_2))[1]
        model2 = get_ext(g_2)["model"]
        p2 = get_ext(g_2)["ps"]
        st2 = get_ext(g_2)["st"]

        y_1, st_1 = Lux.apply(model, input, p1, st1)
        y2, st_2 = Lux.apply(model2, input, p2, st2)
        @test isapprox(y_1, y2)
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
