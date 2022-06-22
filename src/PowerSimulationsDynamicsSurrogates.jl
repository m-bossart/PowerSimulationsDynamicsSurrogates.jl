module PowerSimulationsDynamicsSurrogates

export SteadyStateNODE
export get_SteadyStateNODE_states
export initializer_psids
export get_name #TODO -  remove some getter/setter functions that won't be used
export get_initializer_structure
export get_initializer_parameters
export get_node_structure_exogenous
export get_node_structure_states
export get_node_structure_common
export get_node_parameters
export get_observer_structure
export get_observer_parameters
export get_x_scale
export get_x_bias
export get_exogenous_scale
export get_exogenous_bias
export get_base_power
export get_states
export get_n_states
export get_ext
export get_internal
export set_initializer_structure!
export set_initializer_parameters!
export set_node_structure_exogenous!
export set_node_structure_states!
export set_node_structure_common!
export set_node_parameters!
export set_observer_structure!
export set_observer_parameters!
export set_x_scale!
export set_x_bias!
export set_exogenous_scale!
export set_exogenous_bias!
export set_base_power!
export set_ext

#export SurrogatePerturbations generation exports 
export SurrogatePerturbation
export PVS
export VStep

#export SurrogateOperatingPoints
export SurrogateOperatingPoint
export GenerationLoadScale

#export SurrogateTrainDatasets
export SurrogateTrainDataset
export SteadyStateNODEData

export generate_train_data
export GenerateDataParams

import InfrastructureSystems
import OrdinaryDiffEq
import PowerFlows
import PowerSystems
import PowerSimulationsDynamics
import NLsolve
const IS = InfrastructureSystems
const PSY = PowerSystems
const PSID = PowerSimulationsDynamics

# Write your package code here.
include("SteadyStateNODE/SteadyStateNODE.jl")
include("SteadyStateNODE/utils.jl")
include("generate_data/Perturbations.jl")
include("generate_data/OperatingPointChanges.jl")
include("generate_data/Datasets.jl")

end
