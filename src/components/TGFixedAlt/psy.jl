
mutable struct TGFixedAlt <: PSY.TurbineGov
    efficiency::Float64
    P_ref::Float64
    ext::Dict{String, Any}
    states::Vector{Symbol}
    n_states::Int
    internal::IS.InfrastructureSystemsInternal
end

function TGFixedAlt(efficiency, P_ref = 1.0, ext = Dict{String, Any}())
    TGFixedAlt(efficiency, P_ref, ext, [:P_ref], 1, IS.InfrastructureSystemsInternal())
end

function TGFixedAlt(;
    efficiency,
    P_ref = 1.0,
    ext = Dict{String, Any}(),
    states = [:P_ref],
    n_states = 1,
    internal = IS.InfrastructureSystemsInternal(),
)
    TGFixedAlt(efficiency, P_ref, ext, states, n_states, internal)
end

# Constructor for demo purposes; non-functional.
function TGFixedAlt(::Nothing)
    TGFixedAlt(; efficiency = 0, P_ref = 0, ext = Dict{String, Any}())
end

PSY.get_efficiency(value::TGFixedAlt) = value.efficiency
PSY.get_P_ref(value::TGFixedAlt) = value.P_ref
PSY.get_ext(value::TGFixedAlt) = value.ext
PSY.get_states(value::TGFixedAlt) = value.states
PSY.get_n_states(value::TGFixedAlt) = value.n_states
PSY.get_internal(value::TGFixedAlt) = value.internal
PSY.set_efficiency!(value::TGFixedAlt, val) = value.efficiency = val
PSY.set_P_ref!(value::TGFixedAlt, val) = value.P_ref = val
PSY.set_ext!(value::TGFixedAlt, val) = value.ext = val
