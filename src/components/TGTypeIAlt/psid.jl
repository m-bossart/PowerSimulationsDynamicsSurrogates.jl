function PSID.initialize_tg!(
    device_states,
    p,
    static::PSY.StaticInjection,
    dynamic_device::PSID.DynamicWrapper{PSY.DynamicGenerator{M, S, A, TGTypeIAlt, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[PSID.τm_var]

    tg = PSY.get_prime_mover(dynamic_device)
    #Get Parameters
    params = p[:params][:TurbineGov]
    R = params[:R]
    Tc = params[:Tc]
    T3 = params[:T3]
    T4 = params[:T4]
    T5 = params[:T5]
    V_min = params[:valve_position_limits][:min]
    V_max = params[:valve_position_limits][:max]

    #Get References
    ω_ref = p[:refs][:ω_ref]
    ω0 = 1.0
    function f!(out, x, params)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]
        x_g3 = x[4]

        R = params[:R]
        Tc = params[:Tc]
        T3 = params[:T3]
        T4 = params[:T4]
        T5 = params[:T5]
        ω0 = 1.0
        #Compute auxiliary parameters
        inv_R = R < eps() ? 0.0 : (1.0 / R)
        P_in = P_ref + inv_R * (ω_ref - ω0)

        out[1] = P_in - x_g1
        out[2] = (1.0 - T3 / Tc) * x_g1 - x_g2
        out[3] = (1.0 - T4 / T5) * (x_g2 + (T3 / Tc) * x_g1) - x_g3
        out[4] = x_g3 + (T4 / T5) * (x_g2 + (T3 / Tc) * x_g1) - τm0
    end
    x0 = [
        τm0,
        τm0,
        (1.0 - T3 / Tc) * τm0,
        (1.0 - T4 / T5) * ((1.0 - T3 / Tc) * τm0 + (T3 / Tc) * τm0),
    ]
    prob = NonlinearSolve.NonlinearProblem{true}(f!, x0, params)
    sol = NonlinearSolve.solve(
        prob,
        NonlinearSolve.TrustRegion();
        reltol = PSID.STRICT_NLSOLVE_F_TOLERANCE,
        abstol = PSID.STRICT_NLSOLVE_F_TOLERANCE,
    )
    if !SciMLBase.successful_retcode(sol)
        @warn("Initialization of Turbine Governor $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.u
        #Update Control Refs
        PSY.set_P_ref!(tg, sol_x0[1])
        PSID.set_P_ref!(p, sol_x0[1])
        #Update states
        tg_ix = PSID.get_local_state_ix(dynamic_device, TGTypeIAlt)
        tg_states = @view device_states[tg_ix]
        if (sol_x0[2] > V_max) || (sol_x0[2] < V_min)
            @error(
                "Valve limits for TG in $(PSY.get_name(dynamic_device)) (x_g1 = $(sol_x0[2])), outside its limits V_max = $V_max, Vmin = $V_min. Consider updating the operating point."
            )
        end
        tg_states[1] = sol_x0[2]
        tg_states[2] = sol_x0[3]
        tg_states[3] = sol_x0[4]
        tg_states[4] = sol_x0[1]
    end
    return
end

function PSID.mdl_tg_ode!(
    device_states::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    ω_sys::PSID.ACCEPTED_REAL_TYPES,
    device::PSID.DynamicWrapper{PSY.DynamicGenerator{M, S, A, TGTypeIAlt, P}},
    h,
    t,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Obtain references
    ω_ref = p[:refs][:ω_ref]
    #P_ref = p[:refs][:P_ref]

    #Obtain indices for component w/r to device
    local_ix = PSID.get_local_state_ix(device, TGTypeIAlt)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    x_g1 = internal_states[1]
    x_g2 = internal_states[2]
    x_g3 = internal_states[3]
    P_ref = internal_states[4]
    #Obtain external states inputs for component
    external_ix = PSID.get_input_port_ix(device, TGTypeIAlt)
    ω = @view device_states[external_ix]

    #Get Parameters
    params = p[:params][:TurbineGov]
    R = params[:R]
    Ts = params[:Ts]
    Tc = params[:Tc]
    T3 = params[:T3]
    T4 = params[:T4]
    T5 = params[:T5]
    V_min = params[:valve_position_limits][:min]
    V_max = params[:valve_position_limits][:max]
    inv_R = R < eps() ? 0.0 : (1.0 / R)

    #Compute auxiliary parameters
    P_in_sat = PSID.clamp(P_ref + inv_R * (ω_ref - ω[1]), V_min, V_max)

    #Compute block derivatives
    _, dxg1_dt = PSID.low_pass(P_in_sat, x_g1, 1.0, Ts)
    y_ll, dxg2_dt = PSID.lead_lag(x_g1, x_g2, 1.0, T3, Tc)
    τ_m, dxg3_dt = PSID.lead_lag(y_ll, x_g3, 1.0, T4, T5)

    #Compute 3 States TG ODE:
    output_ode[local_ix[1]] = dxg1_dt
    output_ode[local_ix[2]] = dxg2_dt
    output_ode[local_ix[3]] = dxg3_dt
    output_ode[local_ix[4]] = 0.0

    #Update mechanical torque
    inner_vars[PSID.τm_var] = τ_m

    return
end

function PSID._mechanical_torque(
    tg::TGTypeIAlt,
    name::String,
    res::PSID.SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
    unique_timestamps::Bool,
)
    # Get params
    Tc = PSY.get_Tc(tg)
    T3 = PSY.get_T3(tg)
    T4 = PSY.get_T4(tg)
    T5 = PSY.get_T5(tg)
    # Get state results
    ts, x_g1 = PSID.post_proc_state_series(res, (name, :x_g1), dt, unique_timestamps)
    _, x_g2 = PSID.post_proc_state_series(res, (name, :x_g2), dt, unique_timestamps)
    _, x_g3 = PSID.post_proc_state_series(res, (name, :x_g3), dt, unique_timestamps)
    τm = zeros(length(ts))
    for ix in 1:length(ts)
        y_ll, _ = PSID.lead_lag(x_g1[ix], x_g2[ix], 1.0, T3, Tc)
        τ_out, _ = PSID.lead_lag(y_ll, x_g3[ix], 1.0, T4, T5)
        τm[ix] = τ_out
    end
    return ts, τm
end

PSID.get_params(x::TGTypeIAlt) = (
    R = PSY.get_R(x),
    Ts = PSY.get_Ts(x),
    Tc = PSY.get_Tc(x),
    T3 = PSY.get_T3(x),
    T4 = PSY.get_T4(x),
    T5 = PSY.get_T5(x),
    valve_position_limits = PSY.get_valve_position_limits(x),
)

PSID.get_params_metadata(::TGTypeIAlt) = (
    R = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
    Ts = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, false),
    Tc = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
    T3 = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
    T4 = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
    T5 = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
    valve_position_limits = (
        min = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
        max = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
    ),
)
