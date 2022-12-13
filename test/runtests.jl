using Revise
using NLsolve
using Plots
using PowerFlows
using PowerSystems
using PowerSimulationsDynamics
using PowerSimulationsDynamicsSurrogates
using Sundials
using Flux
using Test
using Random

const PSY = PowerSystems
const PSID = PowerSimulationsDynamics
const PSIDS = PowerSimulationsDynamicsSurrogates

#include data for tests
include("data_tests/data_utils.jl")
include("data_tests/dynamic_test_data.jl")

#run tests
include("test_build_systems.jl")
include("test_onebus.jl")
include("test_SteadyStateNODE.jl")
include("test_SteadyStateNODEObs.jl")
include("test_type_format.jl")
include("test_data_generation.jl")
include("test_ChirpVariableSource.jl")
include("test_FrequencyChirpVariableSource.jl")
