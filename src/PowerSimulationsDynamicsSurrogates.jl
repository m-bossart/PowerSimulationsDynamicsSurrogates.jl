module PowerSimulationsDynamicsSurrogates

export TGFixedAlt
export TGTypeIAlt
export SteamTurbineGov1Alt
export ChirpVariableSource
export FrequencyChirpVariableSource
export CurrentPlayback
export SourceLoad
export SteadyStateNODE
export get_SteadyStateNODE_states
export SolutionPredictionSurrogate
export TerminalDataSurrogate
export get_name
export get_initializer_structure
export get_initializer_parameters
export get_node_structure
export get_node_parameters
export get_observer_structure
export get_observer_parameters
export get_input_min
export get_input_max
export get_target_min
export get_target_max
export get_base_power
export get_states
export get_n_states
export get_ext
export get_internal
export set_initializer_parameters!
export set_node_parameters!
export set_observer_parameters!
export set_base_power!
export set_fc!
export set_ext!

export get_component_parameters
export set_parameter_relative!
export set_parameter_absolute!
export add_aggregate_model!
export generate_empty_plot
export add_data_trace!

# SurrogateParams
export NODEParams
export SteadyStateNODEParams
export SourceParams
export SourceLoadParams
export TerminalDataSurrogateParams
export ClassicGenParams
export GFLParams
export GFMParams
export ZIPParams
export MultiDeviceParams
export MultiDeviceLineParams

# SurrogatePerturbations   
export SurrogatePerturbation
export LineTrip
export PVS
export Chirp
export VStep
export RandomSourceVoltageChange
export RandomSourceFrequencyChange
export RandomLoadTrip
export RandomBranchTrip
export RandomLoadChange

# SurrogateOperatingPoints
export SurrogateOperatingPoint
export ScaleSource
export GenerationLoadScale
export RandomOperatingPointXiao

# SurrogateDatasets
export SurrogateDataset
export AllStatesData
export BusData
export FullSolutionData
export TerminalData

export GenerateDataParams

export generate_surrogate_data
export create_validation_system_from_buses
export create_train_system_from_buses

export to_json_with_surrogates
export deserialize_with_surrogates

import BSON
import DataFrames
import Distributions
import InfrastructureSystems
import JSON3
import SciMLBase
import NonlinearSolve
import OrdinaryDiffEq   #should remove this dependency: pass solver to data generation function
import DelayDiffEq      #should remove this dependency: pass solver to data generation function
import PlotlyJS
import PowerFlows
import PowerSimulationsDynamics
import PowerSystems
import Sundials     #should remove this dependency 
import Random
import ComponentArrays
import Zygote

const IS = InfrastructureSystems
const PSY = PowerSystems
const PSID = PowerSimulationsDynamics

include("utils.jl")
include("components/abstract_component_types.jl")
include("components/TerminalDataSurrogate/psy.jl")
include("components/TerminalDataSurrogate/psid.jl")
include("components/SteadyStateNODE/psy.jl")
include("components/SteadyStateNODE/psid.jl")
include("components/CurrentPlayback/psy.jl")
include("components/CurrentPlayback/psid.jl")
include("components/ChirpVariableSource/psy.jl")
include("components/ChirpVariableSource/psid.jl")
include("components/SourceLoad/psy.jl")
include("components/SourceLoad/psid.jl")
include("components/FrequencyChirpVariableSource/psy.jl")
include("components/FrequencyChirpVariableSource/psid.jl")
include("components/TGFixedAlt/psy.jl")
include("components/TGFixedAlt/psid.jl")
include("components/TGTypeIAlt/psy.jl")
include("components/TGTypeIAlt/psid.jl")
include("components/SteamTurbineGov1Alt/psy.jl")
include("components/SteamTurbineGov1Alt/psid.jl")
include("components/zygote_buffer_compatibility.jl")
include("generate_data/Perturbations.jl")
include("generate_data/OperatingPointChanges.jl")
include("generate_data/datasets/generate.jl")
include("generate_data/datasets/TerminalData.jl")
include("generate_data/datasets/AllStatesData.jl")
include("generate_data/datasets/BusData.jl")
include("generate_data/datasets/FullSolutionData.jl")
include("manipulations/system_manipulations.jl")
include("manipulations/parameter_manipulations.jl")
end
