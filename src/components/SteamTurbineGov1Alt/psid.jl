function PSID.initialize_tg!(
    device_states,
    p,
    ::PSY.StaticInjection,
    dynamic_device::PSID.DynamicWrapper{
        PSY.DynamicGenerator{M, S, A, SteamTurbineGov1Alt, P},
    },
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[PSID.τm_var]
    Δω = 0.0

    tg = PSY.get_prime_mover(dynamic_device)
    #Get Parameters
    params = p[:params][:TurbineGov]
    R = params[:R]
    V_min = params[:valve_position_limits][:min]
    V_max = params[:valve_position_limits][:max]
    inv_R = R < eps() ? 0.0 : (1.0 / R)

    function f!(out, x, params)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]

        R = params[:R]
        T1 = params[:T1]
        V_min = params[:valve_position_limits][:min]
        V_max = params[:valve_position_limits][:max]
        T2 = params[:T2]
        T3 = params[:T3]
        D_T = params[:D_T]
        inv_R = R < eps() ? 0.0 : (1.0 / R)

        ref_in = inv_R * (P_ref - Δω)
        Pm = x_g2 + (T2 / T3) * x_g1

        out[1] = (1.0 / T1) * (ref_in - x_g1) #dx_g1/dt
        x_g1_sat = V_min < x_g1 < V_max ? x_g1 : max(V_min, min(V_max, x_g1))
        out[2] = (1.0 / T3) * (x_g1_sat * (1 - T2 / T3) - x_g2) #dx_g2/dt
        out[3] = (Pm - D_T * Δω) - τm0
    end
    x0 = [1.0 / inv_R, τm0, τm0]
    prob = NonlinearSolve.NonlinearProblem{true}(f!, x0, params)
    sol = NonlinearSolve.solve(
        prob,
        NonlinearSolve.TrustRegion();
        reltol = PSID.STRICT_NLSOLVE_F_TOLERANCE,
        abstol = PSID.STRICT_NLSOLVE_F_TOLERANCE,
    )
    if !SciMLBase.successful_retcode(sol)
        @warn("Initialization in TG failed")
    else
        sol_x0 = sol.u
        if (sol_x0[2] >= V_max + PSID.BOUNDS_TOLERANCE) ||
           (sol_x0[2] <= V_min - PSID.BOUNDS_TOLERANCE)
            @error(
                "Valve limits for TG in $(PSY.get_name(dynamic_device)) (x_g1 = $(sol_x0[2])), outside its limits V_max = $V_max, Vmin = $V_min.  Consider updating the operating point."
            )
        end
        #Update Control Refs
        PSY.set_P_ref!(tg, sol_x0[1])
        PSID.set_P_ref!(p, sol_x0[1])
        #Update states
        tg_ix = PSID.get_local_state_ix(dynamic_device, typeof(tg))
        tg_states = @view device_states[tg_ix]
        tg_states[1] = sol_x0[2]
        tg_states[2] = sol_x0[3]
        tg_states[3] = sol_x0[1]    #DIFF 
    end
    return
end

function PSID.mdl_tg_ode!(
    device_states::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    ω_sys::PSID.ACCEPTED_REAL_TYPES,
    device::PSID.DynamicWrapper{PSY.DynamicGenerator{M, S, A, SteamTurbineGov1Alt, P}},
    h,
    t,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Obtain TG
    tg = PSY.get_prime_mover(device)
    #Obtain references

    #Obtain indices for component w/r to device
    local_ix = PSID.get_local_state_ix(device, typeof(tg))

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    x_g1 = internal_states[1]
    x_g2 = internal_states[2]
    P_ref = internal_states[3]#p[:refs][:P_ref]
    #Obtain external states inputs for component
    external_ix = PSID.get_input_port_ix(device, typeof(tg))
    ω = @view device_states[external_ix]

    #Get Parameters
    params = @view(p[:params][:TurbineGov])
    R = params[:R]
    T1 = params[:T1]
    V_min = params[:valve_position_limits][:min]
    V_max = params[:valve_position_limits][:max]
    T2 = params[:T2]
    T3 = params[:T3]
    D_T = params[:D_T]
    inv_R = R < eps() ? 0.0 : (1.0 / R)

    #Compute auxiliary parameters
    ref_in = inv_R * (P_ref - (ω[1] - 1.0))

    #Compute block derivatives
    x_g1_sat, dxg1_dt = PSID.low_pass_nonwindup(ref_in, x_g1, 1.0, T1, V_min, V_max)
    y_ll, dxg2_dt = PSID.lead_lag(x_g1_sat, x_g2, 1.0, T2, T3)
    P_m = y_ll - D_T * (ω[1] - 1.0)

    #Compute 2 State TG ODE:
    output_ode[local_ix[1]] = dxg1_dt
    output_ode[local_ix[2]] = dxg2_dt
    output_ode[local_ix[3]] = 0.0

    #Update mechanical torque
    inner_vars[PSID.τm_var] = P_m / ω[1]

    return
end

function PSID._mechanical_torque(
    tg::SteamTurbineGov1Alt,
    name::String,
    res::PSID.SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
    unique_timestamps::Bool,
)
    # TODO: This will not plot correctly when changing P_ref in a callback
    # Get params
    setpoints = get_setpoints(res)
    P_ref = setpoints[name]["P_ref"]
    inv_R = PSY.get_R(tg) < eps() ? 0.0 : (1.0 / PSY.get_R(tg))
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)
    V_min, V_max = PSY.get_valve_position_limits(tg)
    T3 = PSY.get_T3(tg)
    D_T = PSY.get_D_T(tg)
    # Get state results
    ts, x_g1 = PSID.post_proc_state_series(res, (name, :x_g1), dt, unique_timestamps)
    _, x_g2 = PSID.post_proc_state_series(res, (name, :x_g2), dt, unique_timestamps)
    _, ω = PSID.post_proc_state_series(res, (name, :ω), dt, unique_timestamps)
    ref_in = inv_R * (P_ref .- (ω .- 1.0))
    τm = zeros(length(ts))
    for ix in 1:length(ts)
        x_g1_sat, _ = PSID.low_pass_nonwindup(ref_in[ix], x_g1[ix], 1.0, T1, V_min, V_max)
        y_ll, _ = PSID.lead_lag(x_g1_sat, x_g2[ix], 1.0, T2, T3)
        P_m = y_ll - D_T * (ω[ix] - 1.0)
        τm[ix] = P_m / ω[ix]
    end
    return ts, τm
end

PSID.get_params(x::SteamTurbineGov1Alt) = (
    R = PSY.get_R(x),
    T1 = PSY.get_T1(x),
    valve_position_limits = PSY.get_valve_position_limits(x),
    T2 = PSY.get_T2(x),
    T3 = PSY.get_T3(x),
    D_T = PSY.get_D_T(x),
)
PSID.get_params_metadata(::SteamTurbineGov1Alt) = (
    R = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
    T1 = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
    valve_position = (
        min = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
        max = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
    ),
    T2 = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
    T3 = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
    D_T = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
)
