using Revise
using Random
using ComponentArrays
using DataFrames
using JSON3
using LinearAlgebra
using DelayDiffEq
using OrdinaryDiffEq
using PowerFlows
using PowerSystems
using PowerSimulationsDynamics
using PowerSimulationsDynamicsSurrogates
using Sundials
using PowerSystemCaseBuilder
using SciMLSensitivity
using Lux
using Test
using Random
using Logging
using PlotlyJS
using Zygote
using Enzyme
Enzyme.API.runtimeActivity!(true)  #Needed for "activity unstable" code: https://enzymead.github.io/Enzyme.jl/stable/faq/
Enzyme.API.looseTypeAnalysis!(true)  #Required for using component arrays with Enzyme
const PSY = PowerSystems
const PSID = PowerSimulationsDynamics
const PSIDS = PowerSimulationsDynamicsSurrogates

#include data for tests
include("data_tests/data_utils.jl")
include("data_tests/dynamic_test_data.jl")

const DISABLED_TEST_FILES = []
macro includetests(testarg...)
    if length(testarg) == 0
        tests = []
    elseif length(testarg) == 1
        tests = testarg[1]
    else
        error("@includetests takes zero or one argument")
    end

    quote
        tests = $tests
        rootfile = @__FILE__
        if length(tests) == 0
            tests = readdir(dirname(rootfile))
            tests = filter(
                f ->
                    startswith(f, "test_") && endswith(f, ".jl") && f != basename(rootfile),
                tests,
            )
        else
            tests = map(f -> string(f, ".jl"), tests)
        end
        println()
        for test in tests
            test ∈ DISABLED_TEST_FILES && continue
            print(splitext(test)[1], ": ")
            include(test)
            println()
        end
    end
end
##
logger = PSY.configure_logging(;
    console_level = Logging.Warn,  # Logging.Error, Logging.Warn, Logging.Debug
    file_level = Logging.Error,
)
#with_logger(logger) do
@time @testset "Begin PowerSimulationsDynamicsSurrogates tests" begin
    @includetests ARGS
end
#run tests
#include("test_CurrentPlayback.jl")
#include("test_TerminalDataSurrogate.jl")
#include("test_serialization.jl")
#include("test_SourceLoad.jl")
#include("test_onebus.jl")
#include("test_SteadyStateNODE.jl")  
#include("test_type_format.jl")
#include("test_data_generation.jl")
#include("test_dataset_aux.jl")  
#include("test_ChirpVariableSource.jl")
#include("test_FrequencyChirpVariableSource.jl")
#end
flush(logger)
close(logger)
