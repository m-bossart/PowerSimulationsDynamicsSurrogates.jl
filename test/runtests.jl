using Revise
using BSON: @save, @load
using NLsolve
using Random
using DataFrames
using JSON3
using LinearAlgebra
using PowerFlows
using PowerSystems
using PowerSimulationsDynamics
using PowerSimulationsDynamicsSurrogates
using Sundials
using Lux
using Flux
using Test
using Random
using Logging
const PSY = PowerSystems
const PSID = PowerSimulationsDynamics
const PSIDS = PowerSimulationsDynamicsSurrogates

#include data for tests
include("data_tests/data_utils.jl")
include("data_tests/dynamic_test_data.jl")

logger = PSY.configure_logging(;
    console_level = Logging.Warn,  # Logging.Error, Logging.Debug
    file_level = Logging.Error,
)
with_logger(logger) do
    #run tests
    include("test_TerminalDataSurrogate.jl")
    include("test_serialization.jl")
    include("test_SourceLoad.jl")
    include("test_build_systems.jl")
    include("test_onebus.jl")
    include("test_SteadyStateNODE.jl")
    include("test_SteadyStateNODEObs.jl")
    include("test_type_format.jl")
    include("test_data_generation.jl")
    include("test_dataset_aux.jl")
    include("test_ChirpVariableSource.jl")
    include("test_FrequencyChirpVariableSource.jl")
end
flush(logger)
close(logger)
