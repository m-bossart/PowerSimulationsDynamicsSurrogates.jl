#Ensure the models are built correctly and match python implementation which was used in training.
#Copy over model.h5 and results.json file from PowerSystemSolutionSurrogates
#Note - this does not test the data scaling; only the core data-driven model

function parse_test_input_output(results_file_path)
    df = copy(JSON3.read(results_file_path))
    dynamic_input = df[Symbol("Input Sample")][Symbol("0")][1]
    dynamic_input = reduce(vcat, dynamic_input')
    static_input = df[Symbol("Input Sample")][Symbol("0")][2]
    static_input = reshape(static_input, (1, length(static_input)))
    output = df[Symbol("Output Sample")][Symbol("0")][1]

    return dynamic_input, static_input, output
end

@testset "Test FFNN" begin
    h5_file_path = joinpath(pwd(), "test", "data_tests", "model_ffnn.h5")
    results_file_path = joinpath(pwd(), "test", "data_tests", "results_ffnn.json")

    dynamic_input, static_input, output_python = parse_test_input_output(results_file_path)
    model_architecture, model_parameters = parse_h5(h5_file_path)
    ext = Dict{String, Any}()

    PSIDS._allocate_model_parameters!(ext, model_architecture, model_parameters)
    @show [size(x) for x in ext["nn"]]
    output_julia =
        PSIDS._call_model(ext["nn"], model_architecture, static_input, dynamic_input)
    @show output_julia, output_python
    @test isapprox(output_python, output_julia, atol = 1e-7)
end

@testset "Test RNN" begin
    h5_file_path = joinpath(pwd(), "test", "data_tests", "model_rnn.h5")
    results_file_path = joinpath(pwd(), "test", "data_tests", "results_rnn.json")

    dynamic_input, static_input, output_python = parse_test_input_output(results_file_path)
    model_architecture, model_parameters = parse_h5(h5_file_path)
    ext = Dict{String, Any}()

    PSIDS._allocate_model_parameters!(ext, model_architecture, model_parameters)
    @show [size(x) for x in ext["nn"]]
    output_julia =
        PSIDS._call_model(ext["nn"], model_architecture, static_input, dynamic_input)
    @show output_julia, output_python
    @test isapprox(output_python, output_julia, atol = 1e-6)
end

@testset "Test LSTM" begin
    h5_file_path = joinpath(pwd(), "test", "data_tests", "model_lstm.h5")
    results_file_path = joinpath(pwd(), "test", "data_tests", "results_lstm.json")

    dynamic_input, static_input, output_python = parse_test_input_output(results_file_path)
    model_architecture, model_parameters = parse_h5(h5_file_path)
    ext = Dict{String, Any}()

    PSIDS._allocate_model_parameters!(ext, model_architecture, model_parameters)
    @show [size(x) for x in ext["nn"]]
    output_julia =
        PSIDS._call_model(ext["nn"], model_architecture, static_input, dynamic_input)
    @show output_julia, output_python
    @test isapprox(output_python, output_julia, atol = 1e-7)
end
