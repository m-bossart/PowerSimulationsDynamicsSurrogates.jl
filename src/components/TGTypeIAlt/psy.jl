mutable struct TGTypeIAlt <: PSY.TurbineGov
    R::Float64
    Ts::Float64
    Tc::Float64
    T3::Float64
    T4::Float64
    T5::Float64
    valve_position_limits::PSY.MinMax
    P_ref::Float64
    ext::Dict{String, Any}
    states::Vector{Symbol}
    n_states::Int
    internal::IS.InfrastructureSystemsInternal
end

function TGTypeIAlt(
    R,
    Ts,
    Tc,
    T3,
    T4,
    T5,
    valve_position_limits,
    P_ref = 1.0,
    ext = Dict{String, Any}(),
)
    TGTypeIAlt(
        R,
        Ts,
        Tc,
        T3,
        T4,
        T5,
        valve_position_limits,
        P_ref,
        ext,
        [:x_g1, :x_g2, :x_g3, :P_ref],
        4,
        IS.InfrastructureSystemsInternal(),
    )
end

function TGTypeIAlt(;
    R,
    Ts,
    Tc,
    T3,
    T4,
    T5,
    valve_position_limits,
    P_ref = 1.0,
    ext = Dict{String, Any}(),
    states = [:x_g1, :x_g2, :x_g3, :P_ref],
    n_states = 4,
    internal = IS.InfrastructureSystemsInternal(),
)
    TGTypeIAlt(
        R,
        Ts,
        Tc,
        T3,
        T4,
        T5,
        valve_position_limits,
        P_ref,
        ext,
        states,
        n_states,
        internal,
    )
end

# Constructor for demo purposes; non-functional.
function TGTypeIAlt(::Nothing)
    TGTypeIAlt(;
        R = 0,
        Ts = 0,
        Tc = 0,
        T3 = 0,
        T4 = 0,
        T5 = 0,
        valve_position_limits = (min = 0.0, max = 0.0),
        P_ref = 0,
        ext = Dict{String, Any}(),
    )
end

PSY.get_R(value::TGTypeIAlt) = value.R
PSY.get_Ts(value::TGTypeIAlt) = value.Ts
PSY.get_Tc(value::TGTypeIAlt) = value.Tc
PSY.get_T3(value::TGTypeIAlt) = value.T3
PSY.get_T4(value::TGTypeIAlt) = value.T4
PSY.get_T5(value::TGTypeIAlt) = value.T5
PSY.get_valve_position_limits(value::TGTypeIAlt) = value.valve_position_limits
PSY.get_P_ref(value::TGTypeIAlt) = value.P_ref
PSY.get_ext(value::TGTypeIAlt) = value.ext
PSY.get_states(value::TGTypeIAlt) = value.states
PSY.get_n_states(value::TGTypeIAlt) = value.n_states
PSY.get_internal(value::TGTypeIAlt) = value.internal

PSY.set_R!(value::TGTypeIAlt, val) = value.R = val
PSY.set_Ts!(value::TGTypeIAlt, val) = value.Ts = val
PSY.set_Tc!(value::TGTypeIAlt, val) = value.Tc = val
PSY.set_T3!(value::TGTypeIAlt, val) = value.T3 = val
PSY.set_T4!(value::TGTypeIAlt, val) = value.T4 = val
PSY.set_T5!(value::TGTypeIAlt, val) = value.T5 = val
PSY.set_valve_position_limits!(value::TGTypeIAlt, val) = value.valve_position_limits = val
PSY.set_P_ref!(value::TGTypeIAlt, val) = value.P_ref = val
PSY.set_ext!(value::TGTypeIAlt, val) = value.ext = val
