"""
    mutable struct SteadyStateNODE <: DynamicInjection
        name::String
        initializer_structure::Vector{Tuple{Int64, Int64, Bool, String}}
        initializer_parameters::Vector{Float64}
        node_structure_exogenous::Vector{Tuple{Int64, Int64, Bool, String}}
        node_structure_states::Vector{Tuple{Int64, Int64, Bool, String}}
        node_structure_common::Vector{Tuple{Int64, Int64, Bool, String}}
        node_parameters::Vector{Float64}
        observer_structure::Vector{Tuple{Int64, Int64, Bool, String}}
        observer_parameters::Vector{Float64}
        x_scale::Vector{Float64}
        x_bias::Vector{Float64}
        exogenous_scale::Vector{Float64}
        exogenous_bias::Vector{Float64}
        base_power::Float64
        states::Vector{Symbol}
        n_states::Int
        ext::Dict{String, Any}
        internal::InfrastructureSystemsInternal
    end

Experimental surrogate

# Arguments
- `name::String`
- `initializer_structure::Vector{Tuple{Int64, Int64, Bool, String}}`: layers of the initializer
- `initializer_parameters::Vector{Float64}`: parameters of the initializer
- `node_structure_exogenous::Vector{Tuple{Int64, Int64, Bool, String}}`: layers of the initializer
- `node_structure_states::Vector{Tuple{Int64, Int64, Bool, String}}`: layers of the initializer
- `node_structure_common::Vector{Tuple{Int64, Int64, Bool, String}}`: layers of the initializer
- `node_parameters::Vector{Float64}`: parameters of the node
- `observer_structure::Vector{Tuple{Int64, Int64, Bool, String}}`: layers of the initializer
- `observer_parameters::Vector{Float64}`: parameters of the node
- `x_scale::Vector{Float64}`: scale for x input
- `x_bias::Vector{Float64}`: bias for x input
- `exogenous_scale::Vector{Float64}`: scale for exogenous input
- `exogenous_bias::Vector{Float64}`: bias for exogenous input
- `base_power::Float64`: Base power
- `states::Vector{Symbol}`: The states of GenericDER depend on the Flags
- `n_states::Int`: The states of GenericDER depend on the Flags
- `ext::Dict{String, Any}`
- `internal::InfrastructureSystemsInternal`: power system internal reference, do not modify
"""
mutable struct SteadyStateNODE <: PSY.DynamicInjection
    name::String
    "layers of the initializer"
    initializer_structure::Vector{Tuple{Int64, Int64, Bool, String}}
    "parameters of the initializer"
    initializer_parameters::Vector{Float64}
    "layers of the initializer"
    node_structure::Vector{Tuple{Int64, Int64, Bool, String}}
    "layers of the initializer"
    node_parameters::Vector{Float64}
    "layers of the initializer"
    observer_structure::Vector{Tuple{Int64, Int64, Bool, String}}
    "parameters of the node"
    observer_parameters::Vector{Float64}
    "scale for x input"
    x_scale::Vector{Float64}
    "bias for x input"
    x_bias::Vector{Float64}
    "scale for exogenous input"
    exogenous_scale::Vector{Float64}
    "bias for exogenous input"
    exogenous_bias::Vector{Float64}
    "Base power"
    base_power::Float64
    "The states of GenericDER depend on the Flags"
    states::Vector{Symbol}
    "The states of GenericDER depend on the Flags"
    n_states::Int
    ext::Dict{String, Any}
    "power system internal reference, do not modify"
    internal::IS.InfrastructureSystemsInternal
end

function SteadyStateNODE(
    name,
    initializer_structure = [(0, 0, true, "init")],
    initializer_parameters = [0.0],
    node_structure = [(0, 0, true, "init")],
    node_parameters = [0.0],
    observer_structure = [(0, 0, true, "init")],
    observer_parameters = [0.0],
    x_scale = [0.0],
    x_bias = [0.0],
    exogenous_scale = [0.0],
    exogenous_bias = [0.0],
    base_power = 100.0,
    ext = Dict{String, Any}(),
)
    SteadyStateNODE(
        name,
        initializer_structure,
        initializer_parameters,
        node_structure,
        node_parameters,
        observer_structure,
        observer_parameters,
        x_scale,
        x_bias,
        exogenous_scale,
        exogenous_bias,
        base_power,
        ext,
        get_SteadyStateNODE_states(initializer_structure[end][2] - 2)[1],
        get_SteadyStateNODE_states(initializer_structure[end][2] - 2)[2],
        IS.InfrastructureSystemsInternal(),
    )
end

function SteadyStateNODE(;
    name,
    initializer_structure = [(0, 0, true, "init")],
    initializer_parameters = [0.0],
    node_structure = [(0, 0, true, "init")],
    node_parameters = [0.0],
    observer_structure = [(0, 0, true, "init")],
    observer_parameters = [0.0],
    x_scale = [0.0],
    x_bias = [0.0],
    exogenous_scale = [0.0],
    exogenous_bias = [0.0],
    base_power = 100.0,
    states = get_SteadyStateNODE_states(initializer_structure[end][2] - 2)[1],
    n_states = get_SteadyStateNODE_states(initializer_structure[end][2] - 2)[2],
    ext = Dict{String, Any}(),
    internal = IS.InfrastructureSystemsInternal(),
)
    SteadyStateNODE(
        name,
        initializer_structure,
        initializer_parameters,
        node_structure,
        node_parameters,
        observer_structure,
        observer_parameters,
        x_scale,
        x_bias,
        exogenous_scale,
        exogenous_bias,
        base_power,
        states,
        n_states,
        ext,
        internal,
    )
end

# Constructor for demo purposes; non-functional.
function SteadyStateNODE(::Nothing)
    SteadyStateNODE(;
        name = "init",
        initializer_structure = [(0, 0, true, "init")],
        initializer_parameters = [0],
        node_structure = [(0, 0, true, "init")],
        node_parameters = [0],
        observer_structure = [(0, 0, true, "init")],
        observer_parameters = Any[0],
        x_scale = [0],
        x_bias = [0],
        exogenous_scale = [0],
        exogenous_bias = [0],
        base_power = 0,
        ext = Dict{String, Any}(),
    )
end

#If function exists in PSY, overload it here. 
"""Get [`SteadyStateNODE`](@ref) `name`."""
PSY.get_name(value::SteadyStateNODE) = value.name
"""Get [`SteadyStateNODE`](@ref) `initializer_structure`."""
get_initializer_structure(value::SteadyStateNODE) = value.initializer_structure
"""Get [`SteadyStateNODE`](@ref) `initializer_parameters`."""
get_initializer_parameters(value::SteadyStateNODE) = value.initializer_parameters
"""Get [`SteadyStateNODE`](@ref) `node_structure_exogenous`."""
get_node_structure(value::SteadyStateNODE) = value.node_structure
"""Get [`SteadyStateNODE`](@ref) `node_structure_states`."""
get_node_parameters(value::SteadyStateNODE) = value.node_parameters
"""Get [`SteadyStateNODE`](@ref) `observer_structure`."""
get_observer_structure(value::SteadyStateNODE) = value.observer_structure
"""Get [`SteadyStateNODE`](@ref) `observer_parameters`."""
get_observer_parameters(value::SteadyStateNODE) = value.observer_parameters
"""Get [`SteadyStateNODE`](@ref) `x_scale`."""
get_x_scale(value::SteadyStateNODE) = value.x_scale
"""Get [`SteadyStateNODE`](@ref) `x_bias`."""
get_x_bias(value::SteadyStateNODE) = value.x_bias
"""Get [`SteadyStateNODE`](@ref) `exogenous_scale`."""
get_exogenous_scale(value::SteadyStateNODE) = value.exogenous_scale
"""Get [`SteadyStateNODE`](@ref) `exogenous_bias`."""
get_exogenous_bias(value::SteadyStateNODE) = value.exogenous_bias
"""Get [`SteadyStateNODE`](@ref) `base_power`."""
PSY.get_base_power(value::SteadyStateNODE) = value.base_power
"""Get [`SteadyStateNODE`](@ref) `states`."""
PSY.get_states(value::SteadyStateNODE) = value.states
"""Get [`SteadyStateNODE`](@ref) `n_states`."""
PSY.get_n_states(value::SteadyStateNODE) = value.n_states
"""Get [`SteadyStateNODE`](@ref) `ext`."""
PSY.get_ext(value::SteadyStateNODE) = value.ext
"""Get [`SteadyStateNODE`](@ref) `internal`."""
PSY.get_internal(value::SteadyStateNODE) = value.internal

"""Set [`SteadyStateNODE`](@ref) `base_power`."""
PSY.set_base_power!(value::SteadyStateNODE, val) = value.base_power = val
"""Set [`SteadyStateNODE`](@ref) `ext`."""
PSY.set_ext!(value::SteadyStateNODE, val) = value.ext = val
"""Set [`SteadyStateNODE`](@ref) `initializer_parameters`."""
set_initializer_parameters!(value::SteadyStateNODE, val) = value.initializer_parameters = val
"""Set [`SteadyStateNODE`](@ref) `node_parameters`."""
set_node_parameters!(value::SteadyStateNODE, val) = value.node_parameters = val
"""Set [`SteadyStateNODE`](@ref) `observer_parameters`."""
set_observer_parameters!(value::SteadyStateNODE, val) = value.observer_parameters = val
