using Revise 
using PowerSimulationsDynamicsSurrogates
using Flux
using PowerSystems
using PowerSimulationsDynamics
using Sundials
PSY = PowerSystems
PSID = PowerSimulationsDynamics
using Test

include("test_SteadyStateNODE.jl")
include("test_ChirpSource.jl")