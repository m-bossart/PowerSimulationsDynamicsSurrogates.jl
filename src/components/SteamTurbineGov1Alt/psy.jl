#This model should behave identically to SteamTurbineGov1Alt
#It is implemented to be able to perturb the power reference with a state perturbation in order to be compatible with ReverseDiffAdjoint()
mutable struct SteamTurbineGov1Alt <: PSY.TurbineGov
    "Droop parameter"
    R::Float64
    "Governor time constant"
    T1::Float64
    "Valve position limits"
    valve_position_limits::PSY.MinMax
    "Lead Lag Lead Time constant "
    T2::Float64
    "Lead Lag Lag Time constant "
    T3::Float64
    "Turbine Damping"
    D_T::Float64
    "Deadband for overspeed"
    DB_h::Float64
    "Deadband for underspeed"
    DB_l::Float64
    "Turbine Rate (MW). If zero, generator base is used."
    T_rate::Float64
    "Reference Power Set-point"
    P_ref::Float64
    ext::Dict{String, Any}
    "The states of the SteamTurbineGov1Alt model are:
 x_g1: Valve Opening,
 x_g2: Lead-lag state
    P_ref: Auxiliary state for P_ref"
    states::Vector{Symbol}
    "TGOV1 has 3 states"
    n_states::Int
    "TGOV1 has 2 differential states"
    states_types::Vector{PSY.StateTypes}
    "power system internal reference, do not modify"
    internal::IS.InfrastructureSystemsInternal
end

function SteamTurbineGov1Alt(
    R,
    T1,
    valve_position_limits,
    T2,
    T3,
    D_T,
    DB_h,
    DB_l,
    T_rate,
    P_ref = 1.0,
    ext = Dict{String, Any}(),
)
    SteamTurbineGov1Alt(
        R,
        T1,
        valve_position_limits,
        T2,
        T3,
        D_T,
        DB_h,
        DB_l,
        T_rate,
        P_ref,
        ext,
        [:x_g1, :x_g2, :P_ref],
        3,
        [
            PSY.StateTypes.Differential,
            PSY.StateTypes.Differential,
            PSY.StateTypes.Differential,
        ],
        IS.InfrastructureSystemsInternal(),
    )
end

function SteamTurbineGov1Alt(;
    R,
    T1,
    valve_position_limits,
    T2,
    T3,
    D_T,
    DB_h,
    DB_l,
    T_rate,
    P_ref = 1.0,
    ext = Dict{String, Any}(),
    states = [:x_g1, :x_g2, :P_ref],
    n_states = 3,
    states_types = [
        PSY.StateTypes.Differential,
        PSY.StateTypes.Differential,
        PSY.StateTypes.Differential,
    ],
    internal = IS.InfrastructureSystemsInternal(),
)
    SteamTurbineGov1Alt(
        R,
        T1,
        valve_position_limits,
        T2,
        T3,
        D_T,
        DB_h,
        DB_l,
        T_rate,
        P_ref,
        ext,
        states,
        n_states,
        states_types,
        internal,
    )
end

# Constructor for demo purposes; non-functional.
function SteamTurbineGov1Alt(::Nothing)
    SteamTurbineGov1Alt(;
        R = 0,
        T1 = 0,
        valve_position_limits = (min = 0.0, max = 0.0),
        T2 = 0,
        T3 = 0,
        D_T = 0,
        DB_h = 0,
        DB_l = 0,
        T_rate = 0,
        P_ref = 0,
        ext = Dict{String, Any}(),
    )
end

PSY.get_R(value::SteamTurbineGov1Alt) = value.R
PSY.get_T1(value::SteamTurbineGov1Alt) = value.T1
PSY.get_valve_position_limits(value::SteamTurbineGov1Alt) = value.valve_position_limits
PSY.get_T2(value::SteamTurbineGov1Alt) = value.T2
PSY.get_T3(value::SteamTurbineGov1Alt) = value.T3
PSY.get_D_T(value::SteamTurbineGov1Alt) = value.D_T
PSY.get_DB_h(value::SteamTurbineGov1Alt) = value.DB_h
PSY.get_DB_l(value::SteamTurbineGov1Alt) = value.DB_l
PSY.get_T_rate(value::SteamTurbineGov1Alt) = value.T_rate
PSY.get_P_ref(value::SteamTurbineGov1Alt) = value.P_ref
PSY.get_ext(value::SteamTurbineGov1Alt) = value.ext
PSY.get_states(value::SteamTurbineGov1Alt) = value.states
PSY.get_n_states(value::SteamTurbineGov1Alt) = value.n_states
PSY.get_states_types(value::SteamTurbineGov1Alt) = value.states_types
PSY.get_internal(value::SteamTurbineGov1Alt) = value.internal

PSY.set_R!(value::SteamTurbineGov1Alt, val) = value.R = val
PSY.set_T1!(value::SteamTurbineGov1Alt, val) = value.T1 = val
PSY.set_valve_position_limits!(value::SteamTurbineGov1Alt, val) =
    value.valve_position_limits = val
PSY.set_T2!(value::SteamTurbineGov1Alt, val) = value.T2 = val
PSY.set_T3!(value::SteamTurbineGov1Alt, val) = value.T3 = val
PSY.set_D_T!(value::SteamTurbineGov1Alt, val) = value.D_T = val
PSY.set_DB_h!(value::SteamTurbineGov1Alt, val) = value.DB_h = val
PSY.set_DB_l!(value::SteamTurbineGov1Alt, val) = value.DB_l = val
PSY.set_T_rate!(value::SteamTurbineGov1Alt, val) = value.T_rate = val
PSY.set_P_ref!(value::SteamTurbineGov1Alt, val) = value.P_ref = val
PSY.set_ext!(value::SteamTurbineGov1Alt, val) = value.ext = val
PSY.set_states_types!(value::SteamTurbineGov1Alt, val) = value.states_types = val
