abstract type LearnedDynamicsSurrogate <: PSY.DynamicInjection end
abstract type LearnedSolutionSurrogate <: PSY.StaticInjection end

abstract type SurrogateModelParams end

mutable struct NODEParams <: SurrogateModelParams
    type::String
    name::String
    n_ports::Int64
    initializer_layer_type::String
    initializer_n_layer::Int64
    initializer_width_layers_relative_input::Int64
    initializer_activation::String
    dynamic_layer_type::String
    dynamic_hidden_states::Int64
    dynamic_n_layer::Int64
    dynamic_width_layers_relative_input::Int64
    dynamic_activation::String
    dynamic_σ2_initialization::Float64
    dynamic_last_layer_bias::Bool
end

function NODEParams(;
    type = "NODEParams",
    name = "surrogate-NODE",
    n_ports = 1,
    initializer_layer_type = "dense",
    initializer_n_layer = 0,
    initializer_width_layers_relative_input = 0,
    initializer_activation = "tanh",
    dynamic_layer_type = "dense",
    dynamic_hidden_states = 5,
    dynamic_n_layer = 1,
    dynamic_width_layers_relative_input = 4,
    dynamic_activation = "tanh",
    dynamic_σ2_initialization = 0.0,
    dynamic_last_layer_bias = false,
)
    NODEParams(
        type,
        name,
        n_ports,
        initializer_layer_type,
        initializer_n_layer,
        initializer_width_layers_relative_input,
        initializer_activation,
        dynamic_layer_type,
        dynamic_hidden_states,
        dynamic_n_layer,
        dynamic_width_layers_relative_input,
        dynamic_activation,
        dynamic_σ2_initialization,
        dynamic_last_layer_bias,
    )
end

mutable struct SteadyStateNODEParams <: SurrogateModelParams
    type::String
    name::String
    n_ports::Int64
    initializer_layer_type::String
    initializer_n_layer::Int64
    initializer_width_layers_relative_input::Int64
    initializer_activation::String
    dynamic_layer_type::String
    dynamic_hidden_states::Int64
    dynamic_n_layer::Int64
    dynamic_width_layers_relative_input::Int64
    dynamic_activation::String
    dynamic_σ2_initialization::Float64
    dynamic_last_layer_bias::Bool
end

function SteadyStateNODEParams(;
    type = "SteadyStateNODEParams",
    name = "surrogate-SteadyStateNODE",
    n_ports = 1,
    initializer_layer_type = "dense",
    initializer_n_layer = 0,
    initializer_width_layers_relative_input = 0,
    initializer_activation = "tanh",
    dynamic_layer_type = "dense",
    dynamic_hidden_states = 5,
    dynamic_n_layer = 1,
    dynamic_width_layers_relative_input = 4,
    dynamic_activation = "tanh",
    dynamic_σ2_initialization = 0.0,
    dynamic_last_layer_bias = false,
)
    SteadyStateNODEParams(
        type,
        name,
        n_ports,
        initializer_layer_type,
        initializer_n_layer,
        initializer_width_layers_relative_input,
        initializer_activation,
        dynamic_layer_type,
        dynamic_hidden_states,
        dynamic_n_layer,
        dynamic_width_layers_relative_input,
        dynamic_activation,
        dynamic_σ2_initialization,
        dynamic_last_layer_bias,
    )
end

mutable struct SteadyStateNODEObsParams <: SurrogateModelParams
    type::String
    name::String
    n_ports::Int64
    initializer_layer_type::String
    initializer_n_layer::Int64
    initializer_width_layers_relative_input::Int64
    initializer_activation::String
    dynamic_layer_type::String
    dynamic_hidden_states::Int64
    dynamic_n_layer::Int64
    dynamic_width_layers_relative_input::Int64
    dynamic_activation::String
    dynamic_σ2_initialization::Float64
    dynamic_last_layer_bias::Bool
    observation_layer_type::String
    observation_n_layer::Int64
    observation_width_layers_relative_input::Int64
    observation_activation::String
end

function SteadyStateNODEObsParams(;
    type = "SteadyStateNODEObsParams",
    name = "surrogate-SteadyStateNODEObs",
    n_ports = 1,
    initializer_layer_type = "dense",
    initializer_n_layer = 0,
    initializer_width_layers_relative_input = 0,
    initializer_activation = "tanh",
    dynamic_layer_type = "dense",
    dynamic_hidden_states = 5,
    dynamic_n_layer = 1,
    dynamic_width_layers_relative_input = 4,
    dynamic_activation = "tanh",
    dynamic_σ2_initialization = 0.0,
    dynamic_last_layer_bias = false,
    observation_layer_type = "dense",
    observation_n_layer = 0,
    observation_width_layers_relative_input = 0,
    observation_activation = "tanh",
)
    SteadyStateNODEObsParams(
        type,
        name,
        n_ports,
        initializer_layer_type,
        initializer_n_layer,
        initializer_width_layers_relative_input,
        initializer_activation,
        dynamic_layer_type,
        dynamic_hidden_states,
        dynamic_n_layer,
        dynamic_width_layers_relative_input,
        dynamic_activation,
        dynamic_σ2_initialization,
        dynamic_last_layer_bias,
        observation_layer_type,
        observation_n_layer,
        observation_width_layers_relative_input,
        observation_activation,
    )
end

mutable struct ClassicGenParams <: SurrogateModelParams
    type::String
    name::String
end

function ClassicGenParams(; type = "ClassicGenParams", name = "surrogate-ClassicGen")
    ClassicGenParams(type, name)
end

mutable struct GFLParams <: SurrogateModelParams
    type::String
    name::String
end

function GFLParams(; type = "GFLParams", name = "surrogate-GFLParams")
    GFLParams(type, name)
end

mutable struct GFMParams <: SurrogateModelParams
    type::String
    name::String
end

function GFMParams(; type = "GFMParams", name = "surrogate-GFMParams")
    GFMParams(type, name)
end

mutable struct ZIPParams <: SurrogateModelParams
    type::String
    name::String
end

function ZIPParams(; type = "ZIPParams", name = "surrogate-ZIPParams")
    ZIPParams(type, name)
end

mutable struct MultiDeviceParams <: SurrogateModelParams
    type::String
    name::String
    static_devices::Vector{SurrogateModelParams}
    dynamic_devices::Vector{SurrogateModelParams}
end

function MultiDeviceParams(;
    type = "MultiDeviceParams",
    name = "surrogate-MultiDeviceParams",
    static_devices = SurrogateModelParams[ZIPParams()],
    dynamic_devices = SurrogateModelParams[GFLParams(), GFMParams()],
)
    MultiDeviceParams(type, name, static_devices, dynamic_devices)
end

mutable struct MultiDeviceLineParams <: SurrogateModelParams
    type::String
    name::String
    static_devices::Vector{SurrogateModelParams}
    dynamic_devices::Vector{SurrogateModelParams}
end

function MultiDeviceLineParams(;
    type = "MultiDeviceLineParams",
    name = "surrogate-MultiDeviceLineParams",
    static_devices = SurrogateModelParams[ZIPParams()],
    dynamic_devices = SurrogateModelParams[GFLParams(), GFMParams()],
)
    MultiDeviceLineParams(type, name, static_devices, dynamic_devices)
end
