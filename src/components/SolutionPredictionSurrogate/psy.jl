"""
    mutable struct SolutionPredictionSurrogate <: PSY.StaticInjection
        name::String
        available::Bool
        bus::Bus
        active_power::Float64
        reactive_power::Float64
        nn_structure::Vector{Tuple{Int64, Int64, Bool, String}}
        nn_parameters::Vector{Float64}
        input_min::Vector{Float64}
        input_max::Vector{Float64}
        input_lims::Tuple{Float64, Float64}
        target_min::Vector{Float64}
        target_max::Vector{Float64}
        target_lims::Tuple{Float64, Float64}
        services::Vector{Service}
        ext::Dict{String, Any}
        internal::IS.InfrastructureSystemsInternal
    end

Experimental surrogate

# Arguments
- `name::String`
- `available::Bool`: 
- `bus::Bus`: 
- `active_power::Float64`: 
- `reactive_power::Float64`:
- `nn_structure::Vector{Tuple{Int64, Int64, Bool, String}}`: structure of the neural network layers 
- `nn_parameters::Vector{Float64}`: parameters of the neural network 
- `input_min::Vector{Float64}`: minimum values of inputs
- `input_max::Vector{Float64}`: maximum values of inputs
- `input_lims::Tuple{Float64, Float64}`: limits for inputs
- `target_min::Vector{Float64}`: minimum values of targets
- `target_max::Vector{Float64}`: maximum values of targets
- `target_lims::Tuplse{Float64, Float64`: limits for targets
- `services::Vector{Service}`
- `ext::Dict{String, Any}`
- `internal::InfrastructureSystemsInternal`: power system internal reference, do not modify
"""
mutable struct SolutionPredictionSurrogate <: LearnedSolutionSurrogate
    name::String
    available::Bool
    bus::PSY.Bus
    active_power::Float64
    reactive_power::Float64
    active_power_limits::PSY.MinMax
    reactive_power_limits::Union{Nothing, PSY.MinMax}
    internal_voltage::Float64   #need to be able to place this component at any type of bus
    internal_angle::Float64     #need to be able to place this component at any type of bus
    nn_structure::Vector{Tuple{Int64, Int64, Bool, String}}
    nn_parameters::Vector{Float64}
    nn_type::Symbol #new
    length_cache::Int64 #new
    nn_features::Symbol #new
    input_min::Vector{Float64}
    input_max::Vector{Float64}
    input_lims::Tuple{Float64, Float64}
    target_min::Vector{Float64}
    target_max::Vector{Float64}
    target_lims::Tuple{Float64, Float64}
    services::Vector{PSY.Service}
    ext::Dict{String, Any}
    internal::IS.InfrastructureSystemsInternal
end

function SolutionPredictionSurrogate(
    name,
    available,
    bus,
    active_power,
    reactive_power,
    active_power_limits,
    reactive_power_limits,
    internal_voltage,
    internal_angle,
    nn_structure,
    nn_parameters,
    nn_type,
    length_cache,
    nn_features,
    input_min,
    input_max,
    input_lims,
    target_min,
    target_max,
    target_lims,
    services = PSY.Device[],
    ext = Dict{String, Any}(),
)
    SolutionPredictionSurrogate(
        name,
        available,
        bus,
        active_power,
        reactive_power,
        active_power_limits,
        reactive_power_limits,
        internal_voltage,
        internal_angle,
        nn_structure,
        nn_parameters,
        nn_type,
        length_cache,
        nn_features,
        input_min,
        input_max,
        input_lims,
        target_min,
        target_max,
        target_lims,
        services,
        ext,
        IS.InfrastructureSystemsInternal(),
    )
end

function SolutionPredictionSurrogate(;
    name,
    available,
    bus,
    active_power,
    reactive_power,
    active_power_limits,
    reactive_power_limits,
    internal_voltage,
    internal_angle,
    nn_structure,
    nn_parameters,
    nn_type,
    length_cache,
    nn_features,
    input_min,
    input_max,
    input_lims,
    target_min,
    target_max,
    target_lims,
    services = PSY.Device[],
    ext = Dict{String, Any}(),
)
    SolutionPredictionSurrogate(
        name,
        available,
        bus,
        active_power,
        reactive_power,
        active_power_limits,
        reactive_power_limits,
        internal_voltage,
        internal_angle,
        nn_structure,
        nn_parameters,
        nn_type,
        length_cache,
        nn_features,
        input_min,
        input_max,
        input_lims,
        target_min,
        target_max,
        target_lims,
        services,
        ext,
        IS.InfrastructureSystemsInternal(),
    )
end

# Constructor for demo purposes; non-functional.
function SolutionPredictionSurrogate(::Nothing)
    SolutionPredictionSurrogate(;
        name = "init",
        available = false,
        bus = PSY.Bus(nothing),
        active_power = 0.0,
        reactive_power = 0.0,
        active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = 0.0, max = 1.0),
        internal_voltage = 0.0,
        internal_angle = 0.0,
        nn_structure = [(0, 0, true, "init")],
        nn_parameters = [0],
        nn_type = :dense,
        length_cache = 0,
        nn_features = :direct,
        input_min = [0.0],
        input_max = [0.0],
        input_lims = (-1, 1),
        target_min = [0.0],
        target_max = [0.0],
        target_lims = (-1, 1),
        ext = Dict{String, Any}(),
    )
end

#If function exists in PSY, overload it here. 
"""Get [`SolutionPredictionSurrogate`](@ref) `name`."""
PSY.get_name(value::SolutionPredictionSurrogate) = value.name
"""Get [`SolutionPredictionSurrogate`](@ref) `available`."""
PSY.get_available(value::SolutionPredictionSurrogate) = value.available
"""Get [`SolutionPredictionSurrogate`](@ref) `bus`."""
PSY.get_bus(value::SolutionPredictionSurrogate) = value.bus
"""Get [`SolutionPredictionSurrogate`](@ref) `active_power`."""
PSY.get_active_power(value::SolutionPredictionSurrogate) = value.active_power
"""Get [`SolutionPredictionSurrogate`](@ref) `reactive_power`."""
PSY.get_reactive_power(value::SolutionPredictionSurrogate) = value.reactive_power
"""Get [`SolutionPredictionSurrogate`](@ref) `active_power_limits`."""
PSY.get_active_power_limits(value::SolutionPredictionSurrogate) = value.active_power_limits
"""Get [`SolutionPredictionSurrogate`](@ref) `reactive_power_limits`."""
PSY.get_reactive_power_limits(value::SolutionPredictionSurrogate) =
    value.reactive_power_limits
"""Get [`SolutionPredictionSurrogate`](@ref) `internal_voltage`."""
PSY.get_internal_voltage(value::SolutionPredictionSurrogate) = value.internal_voltage
"""Get [`SolutionPredictionSurrogate`](@ref) `internal_angle`."""
PSY.get_internal_angle(value::SolutionPredictionSurrogate) = value.internal_angle
"""Get [`SolutionPredictionSurrogate`](@ref) `nn_structure`."""
get_nn_structure(value::SolutionPredictionSurrogate) = value.nn_structure
"""Get [`SolutionPredictionSurrogate`](@ref) `nn_parameters`."""
get_nn_parameters(value::SolutionPredictionSurrogate) = value.nn_parameters
"""Get [`SolutionPredictionSurrogate`](@ref) `nn_type`."""
get_nn_types(value::SolutionPredictionSurrogate) = value.nn_type
"""Get [`SolutionPredictionSurrogate`](@ref) `length_cache`."""
get_length_cache(value::SolutionPredictionSurrogate) = value.length_cache
"""Get [`SolutionPredictionSurrogate`](@ref) `nn_features`."""
get_nn_features(value::SolutionPredictionSurrogate) = value.nn_features
"""Get [`SolutionPredictionSurrogate`](@ref) `input_min`."""
get_input_min(value::SolutionPredictionSurrogate) = value.input_min
"""Get [`SolutionPredictionSurrogate`](@ref) `input_max`."""
get_input_max(value::SolutionPredictionSurrogate) = value.input_max
"""Get [`SolutionPredictionSurrogate`](@ref) `input_lims`."""
get_input_lims(value::SolutionPredictionSurrogate) = value.input_lims
"""Get [`SolutionPredictionSurrogate`](@ref) `target_min`."""
get_target_min(value::SolutionPredictionSurrogate) = value.target_min
"""Get [`SolutionPredictionSurrogate`](@ref) `target_max`."""
get_target_max(value::SolutionPredictionSurrogate) = value.target_max
"""Get [`SolutionPredictionSurrogate`](@ref) `target_lims`."""
get_target_lims(value::SolutionPredictionSurrogate) = value.target_lims
"""Get [`SolutionPredictionSurrogate`](@ref) `ext`."""
PSY.get_ext(value::SolutionPredictionSurrogate) = value.ext
"""Get [`SolutionPredictionSurrogate`](@ref) `internal`."""
PSY.get_internal(value::SolutionPredictionSurrogate) = value.internal

"""Set [`SolutionPredictionSurrogate`](@ref) `available`."""
PSY.set_available!(value::SolutionPredictionSurrogate, val) = value.available = val
"""Set [`SolutionPredictionSurrogate`](@ref) `bus`."""
PSY.set_bus!(value::SolutionPredictionSurrogate, val) = value.bus = val
"""Set [`SolutionPredictionSurrogate`](@ref) `active_power`."""
PSY.set_active_power!(value::SolutionPredictionSurrogate, val) = value.active_power = val
"""Set [`SolutionPredictionSurrogate`](@ref) `reactive_power`."""
PSY.set_reactive_power!(value::SolutionPredictionSurrogate, val) =
    value.reactive_power = val
"""Set [`SolutionPredictionSurrogate`](@ref) `active_power_limits`."""
PSY.set_active_power_limits!(value::SolutionPredictionSurrogate, val) =
    value.active_power_limits = val
"""Set [`SolutionPredictionSurrogate`](@ref) `reactive_power_limits`."""
PSY.set_reactive_power_limits!(value::SolutionPredictionSurrogate, val) =
    value.reactive_power_limits = val
"""Set [`SolutionPredictionSurrogate`](@ref) `internal_voltage`."""
PSY.set_internal_voltage!(value::SolutionPredictionSurrogate, val) =
    value.internal_voltage = val
"""Set [`SolutionPredictionSurrogate`](@ref) `internal_angle`."""
PSY.set_internal_angle!(value::SolutionPredictionSurrogate, val) =
    value.internal_angle = val
"""Set [`SolutionPredictionSurrogate`](@ref) `ext`."""
PSY.set_ext!(value::SolutionPredictionSurrogate, val) = value.ext = val
"""Set [`SolutionPredictionSurrogate`](@ref) `initializer_parameters`."""
set_nn_parameters!(value::SolutionPredictionSurrogate, val) =
    value.initializer_parameters = val
