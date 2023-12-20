"""
    mutable struct SteadyStateNODEObs <: DynamicInjection
        name::String
        initializer_structure::Vector{Tuple{Int64, Int64, Bool, String}}
        initializer_parameters::Vector{Float64}
        node_structure::Vector{Tuple{Int64, Int64, Bool, String}}
        node_parameters::Vector{Float64}
        observer_structure::Vector{Tuple{Int64, Int64, Bool, String}}
        observer_parameters::Vector{Float64}
        input_min::Vector{Float64}
        input_max::Vector{Float64}
        input_lims::Tuple{Float64, Float64}
        target_min::Vector{Float64}
        target_max::Vector{Float64}
        target_lims::Tuple{Float64, Float64}
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
- `node_structure::Vector{Tuple{Int64, Int64, Bool, String}}`: layers of the node 
- `node_parameters::Vector{Float64}`: parameters of the node
- `observer_structure::Vector{Tuple{Int64, Int64, Bool, String}}`: layers of the initializer
- `observer_parameters::Vector{Float64}`: parameters of the node
- `input_min::Vector{Float64}`: minimum values of inputs
- `input_max::Vector{Float64}`: maximum values of inputs
- `input_lims::Tuple{Float64, Float64}`: limits for inputs
- `target_min::Vector{Float64}`: minimum values of targets
- `target_max::Vector{Float64}`: maximum values of targets
- `target_lims::Tuplse{Float64, Float64`: limits for targets
- `base_power::Float64`: Base power
- `states::Vector{Symbol}`: The states of GenericDER depend on the Flags
- `n_states::Int`: The states of GenericDER depend on the Flags
- `ext::Dict{String, Any}`
- `internal::InfrastructureSystemsInternal`: power system internal reference, do not modify
"""
mutable struct SteadyStateNODEObs <: TimeSteppingSurrogate
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
    "minimum values of inputs"
    input_min::Vector{Float64}
    "maximum values of inputs"
    input_max::Vector{Float64}
    "limits for inputs"
    input_lims::Tuple{Float64, Float64}
    "minimum values of targets"
    target_min::Vector{Float64}
    "maximum values of targets"
    target_max::Vector{Float64}
    "limits for targets"
    target_lims::Tuple{Float64, Float64}
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

function SteadyStateNODEObs(
    name,
    initializer_structure = [(0, 0, true, "init")], #input dim, output dim, bias?, activation 
    initializer_parameters = [0.0],
    node_structure = [(0, 0, true, "init")],
    node_parameters = [0.0],
    observer_structure = [(0, 0, true, "init")],
    observer_parameters = [0.0],
    input_min = [0.0],
    input_max = [0.0],
    input_lims = (-1, 1),
    target_min = [0.0],
    target_max = [0.0],
    target_lims = (-1, 1),
    base_power = 100.0,
    ext = Dict{String, Any}(),
)
    SteadyStateNODEObs(
        name,
        initializer_structure,
        initializer_parameters,
        node_structure,
        node_parameters,
        observer_structure,
        observer_parameters,
        input_min,
        input_max,
        input_lims,
        target_min,
        target_max,
        target_lims,
        base_power,
        ext,
        get_SteadyStateNODEObs_states(initializer_structure[end][2] - 2)[1],
        get_SteadyStateNODEObs_states(initializer_structure[end][2] - 2)[2],
        IS.InfrastructureSystemsInternal(),
    )
end

function SteadyStateNODEObs(;
    name,
    initializer_structure = [(0, 0, true, "init")],
    initializer_parameters = [0.0],
    node_structure = [(0, 0, true, "init")],
    node_parameters = [0.0],
    observer_structure = [(0, 0, true, "init")],
    observer_parameters = [0.0],
    input_min = [0.0],
    input_max = [0.0],
    input_lims = (-1, 1),
    target_min = [0.0],
    target_max = [0.0],
    target_lims = (-1, 1),
    base_power = 100.0,
    states = get_SteadyStateNODEObs_states(initializer_structure[end][2] - 2)[1],
    n_states = get_SteadyStateNODEObs_states(initializer_structure[end][2] - 2)[2],
    ext = Dict{String, Any}(),
    internal = IS.InfrastructureSystemsInternal(),
)
    SteadyStateNODEObs(
        name,
        initializer_structure,
        initializer_parameters,
        node_structure,
        node_parameters,
        observer_structure,
        observer_parameters,
        input_min,
        input_max,
        input_lims,
        target_min,
        target_max,
        target_lims,
        base_power,
        states,
        n_states,
        ext,
        internal,
    )
end

# Constructor for demo purposes; non-functional.
function SteadyStateNODEObs(::Nothing)
    SteadyStateNODEObs(;
        name = "init",
        initializer_structure = [(0, 0, true, "init")],
        initializer_parameters = [0],
        node_structure = [(0, 0, true, "init")],
        node_parameters = [0],
        observer_structure = [(0, 0, true, "init")],
        observer_parameters = Any[0],
        input_min = [0.0],
        input_max = [0.0],
        input_lims = (-1, 1),
        target_min = [0.0],
        target_max = [0.0],
        target_lims = (-1, 1),
        base_power = 0,
        ext = Dict{String, Any}(),
    )
end

#If function exists in PSY, overload it here. 
"""Get [`SteadyStateNODEObs`](@ref) `name`."""
PSY.get_name(value::SteadyStateNODEObs) = value.name
"""Get [`SteadyStateNODEObs`](@ref) `initializer_structure`."""
get_initializer_structure(value::SteadyStateNODEObs) = value.initializer_structure
"""Get [`SteadyStateNODEObs`](@ref) `initializer_parameters`."""
get_initializer_parameters(value::SteadyStateNODEObs) = value.initializer_parameters
"""Get [`SteadyStateNODEObs`](@ref) `node_structure`."""
get_node_structure(value::SteadyStateNODEObs) = value.node_structure
"""Get [`SteadyStateNODEObs`](@ref) `node_structure_states`."""
get_node_parameters(value::SteadyStateNODEObs) = value.node_parameters
"""Get [`SteadyStateNODEObs`](@ref) `observer_structure`."""
get_observer_structure(value::SteadyStateNODEObs) = value.observer_structure
"""Get [`SteadyStateNODEObs`](@ref) `observer_parameters`."""
get_observer_parameters(value::SteadyStateNODEObs) = value.observer_parameters
"""Get [`SteadyStateNODEObs`](@ref) `input_min`."""
get_input_min(value::SteadyStateNODEObs) = value.input_min
"""Get [`SteadyStateNODEObs`](@ref) `input_max`."""
get_input_max(value::SteadyStateNODEObs) = value.input_max
"""Get [`SteadyStateNODEObs`](@ref) `input_lims`."""
get_input_lims(value::SteadyStateNODEObs) = value.input_lims
"""Get [`SteadyStateNODEObs`](@ref) `target_min`."""
get_target_min(value::SteadyStateNODEObs) = value.target_min
"""Get [`SteadyStateNODEObs`](@ref) `target_max`."""
get_target_max(value::SteadyStateNODEObs) = value.target_max
"""Get [`SteadyStateNODEObs`](@ref) `target_lims`."""
get_target_lims(value::SteadyStateNODEObs) = value.target_lims
"""Get [`SteadyStateNODEObs`](@ref) `base_power`."""
PSY.get_base_power(value::SteadyStateNODEObs) = value.base_power
"""Get [`SteadyStateNODEObs`](@ref) `states`."""
PSY.get_states(value::SteadyStateNODEObs) = value.states
"""Get [`SteadyStateNODEObs`](@ref) `n_states`."""
PSY.get_n_states(value::SteadyStateNODEObs) = value.n_states
"""Get [`SteadyStateNODEObs`](@ref) `ext`."""
PSY.get_ext(value::SteadyStateNODEObs) = value.ext
"""Get [`SteadyStateNODEObs`](@ref) `internal`."""
PSY.get_internal(value::SteadyStateNODEObs) = value.internal

"""Set [`SteadyStateNODEObs`](@ref) `base_power`."""
PSY.set_base_power!(value::SteadyStateNODEObs, val) = value.base_power = val
"""Set [`SteadyStateNODEObs`](@ref) `ext`."""
PSY.set_ext!(value::SteadyStateNODEObs, val) = value.ext = val
"""Set [`SteadyStateNODEObs`](@ref) `initializer_parameters`."""
set_initializer_parameters!(value::SteadyStateNODEObs, val) =
    value.initializer_parameters = val
"""Set [`SteadyStateNODEObs`](@ref) `node_parameters`."""
set_node_parameters!(value::SteadyStateNODEObs, val) = value.node_parameters = val
"""Set [`SteadyStateNODEObs`](@ref) `observer_parameters`."""
set_observer_parameters!(value::SteadyStateNODEObs, val) = value.observer_parameters = val

function get_SteadyStateNODEObs_states(dim_r::Int64)
    states = Symbol[]
    for i in 1:dim_r
        push!(states, Symbol(string("r", i)))
    end
    return states, dim_r
end
