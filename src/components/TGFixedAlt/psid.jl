function PSID.initialize_tg!(
    device_states,
    p,
    ::PSY.StaticInjection,
    dynamic_device::PSID.DynamicWrapper{PSY.DynamicGenerator{M, S, A, TGFixedAlt, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}
    tg = PSY.get_prime_mover(dynamic_device)
    eff = p[:params][:TurbineGov][:efficiency]
    τm0 = inner_vars[PSID.τm_var]

    P_ref = τm0 / eff
    #Update Control Refs
    PSY.set_P_ref!(tg, P_ref)
    PSID.set_P_ref!(p, P_ref)
    tg = PSY.get_prime_mover(dynamic_device)
    tg_ix = PSID.get_local_state_ix(dynamic_device, typeof(tg))
    tg_states = @view device_states[tg_ix]
    tg_states[1] = P_ref
    return
end

function PSID.mdl_tg_ode!(
    device_states::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    ω_sys::PSID.ACCEPTED_REAL_TYPES,
    device::PSID.DynamicWrapper{PSY.DynamicGenerator{M, S, A, TGFixedAlt, P}},
    h,
    t,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}
    tg = PSY.get_prime_mover(device)
    local_ix = PSID.get_local_state_ix(device, typeof(tg))
    efficiency = p[:params][:TurbineGov][:efficiency]
    P_ref = device_states[local_ix[1]] #p[:refs][:P_ref]
    inner_vars[PSID.τm_var] = P_ref * efficiency

    output_ode[local_ix[1]] = 0.0
    return
end

function PSID._mechanical_torque(
    tg::TGFixedAlt,
    name::String,
    res::PSID.SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
    unique_timestamps::Bool,
)
    # TODO: This will not plot correctly when changing P_ref in a callback
    ts, _ = PSID._post_proc_state_series(res.solution, 1, dt, unique_timestamps)
    setpoints = PSID.get_setpoints(res)
    P_ref = setpoints[name]["P_ref"]
    efficiency = PSY.get_efficiency(tg)
    τm0 = P_ref * efficiency
    return ts, τm0 * ones(length(ts))
end

PSID.get_params(x::TGFixedAlt) = (; efficiency = PSY.get_efficiency(x))
PSID.get_params_metadata(::TGFixedAlt) =
    (; efficiency = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true))
