statistics_filepath = joinpath(pwd(), "test", "data_tests", "statistics.csv")

@testset "MinMaxScaler" begin
    scaler = PSIDS.parse_min_max(statistics_filepath)
    Random.seed!(1)
    input = rand(1, 8)
    scaled = PSIDS._input_scale(scaler, input, 1:8)
    output = PSIDS._target_scale_inverse(scaler, scaled')
    @test isapprox(vec(input), vec(output))
end

@testset "StandardScaler" begin
    scaler = PSIDS.parse_standard(statistics_filepath)
    Random.seed!(1)
    input = rand(1, 8)
    scaled = PSIDS._input_scale(scaler, input, 1:8)
    output = PSIDS._target_scale_inverse(scaler, scaled')
    @test isapprox(vec(input), vec(output))
end
