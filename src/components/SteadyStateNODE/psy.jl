
mutable struct SteadyStateNODE <: TimeSteppingSurrogate
    name::String
    base_power::Float64
    states::Vector{Symbol}
    n_states::Int
    states_types::Vector{PSY.StateTypes}
    ext::Dict{String, Any}
    internal::IS.InfrastructureSystemsInternal
end

function SteadyStateNODE(name, n_states, ext = Dict{String, Any}())
    SteadyStateNODE(
        name,
        100.0,
        get_SteadyStateNODE_states(n_states),
        n_states,
        fill(PSY.StateTypes.Differential, n_states),
        ext,
        IS.InfrastructureSystemsInternal(),
    )
end

function SteadyStateNODE(;
    name,
    base_power = 100.0,
    states = get_SteadyStateNODE_states(1),
    n_states = 1,
    states_types = [PSY.StateTypes.Differential],
    ext = Dict{String, Any}(),
    internal = IS.InfrastructureSystemsInternal(),
)
    SteadyStateNODE(name, base_power, states, n_states, states_types, ext, internal)
end

# Constructor for demo purposes; non-functional.
function SteadyStateNODE(::Nothing)
    SteadyStateNODE(;
        name = "init",
        base_power = 100.0,
        states = [:r1],
        n_states = 1,
        states_types = [PSY.StateTypes.Differential],
        ext = Dict{String, Any}(),
        internal = IS.InfrastructureSystemsInternal(),
    )
end

#If function exists in PSY, overload it here. 
PSY.get_name(value::SteadyStateNODE) = value.name
PSY.get_base_power(value::SteadyStateNODE) = value.base_power
PSY.get_states(value::SteadyStateNODE) = value.states
PSY.get_n_states(value::SteadyStateNODE) = value.n_states
PSY.get_states_types(value::SteadyStateNODE) = value.states_types
PSY.get_ext(value::SteadyStateNODE) = value.ext
PSY.get_internal(value::SteadyStateNODE) = value.internal

PSY.set_base_power!(value::SteadyStateNODE, val) = value.base_power = val
PSY.set_ext!(value::SteadyStateNODE, val) = value.ext = val

function get_SteadyStateNODE_states(dim_r::Int64)
    states = Symbol[]
    for i in 1:dim_r
        push!(states, Symbol(string("r", i)))
    end
    return states
end
