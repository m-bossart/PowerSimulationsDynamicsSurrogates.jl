PSID.is_valid(::ChirpVariableSource) = nothing

function PSID._get_frequency_state(d::PSID.DynamicWrapper{ChirpVariableSource})
    return 0
end

function PSID.device_mass_matrix_entries!(
    mass_matrix::PSID.AbstractArray,
    dynamic_device::PSID.DynamicWrapper{ChirpVariableSource},
)
    global_index = PSID.get_global_index(dynamic_device)
    mass_matrix_chirp_entries!(mass_matrix, dynamic_device, global_index)
    return
end

function mass_matrix_chirp_entries!(
    mass_matrix,
    chirp::PSID.DynamicWrapper{ChirpVariableSource},
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    @debug "Using default mass matrix entries $chirp"
end

function PSID.device!(
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    ::AbstractArray{T},
    ::AbstractArray{T},
    dynamic_device::PSID.DynamicWrapper{ChirpVariableSource},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES}

    # Internal Voltage states
    V_R = device_states[1] * cos(device_states[2])
    V_I = device_states[1] * sin(device_states[2])

    tstart = get_tstart(PSID.get_device(dynamic_device))
    N = get_N(PSID.get_device(dynamic_device))
    ω1 = get_ω1(PSID.get_device(dynamic_device))
    ω2 = get_ω2(PSID.get_device(dynamic_device))
    V_amp = get_V_amp(PSID.get_device(dynamic_device))
    θ_amp = get_θ_amp(PSID.get_device(dynamic_device))

    if t >= tstart && t < N
        output_ode[1] =
            V_amp *
            cos(ω1 * (t - tstart) + (ω2 - ω1) * (t - tstart)^2 / (2 * N)) *
            (ω1 + (ω2 - ω1) * (t - tstart) / N)
        output_ode[2] =
            θ_amp *
            cos(ω1 * (t - tstart) + (ω2 - ω1) * (t - tstart)^2 / (2 * N)) *
            (ω1 + (ω2 - ω1) * (t - tstart) / N)
    else
        output_ode[1] = 0
        output_ode[2] = 0
    end

    #update current
    R_th = PSY.get_R_th(PSID.get_device(dynamic_device))
    X_th = PSY.get_X_th(PSID.get_device(dynamic_device))
    Zmag = R_th^2 + X_th^2
    current_r[1] += R_th * (V_R - voltage_r[1]) / Zmag + X_th * (V_I - voltage_i[1]) / Zmag #in system pu flowing out
    current_i[1] += R_th * (V_I - voltage_i[1]) / Zmag - X_th * (V_R - voltage_r[1]) / Zmag #in system pu flowing out

    return
end

function PSID.initialize_dynamic_device!(
    dynamic_device::PSID.DynamicWrapper{ChirpVariableSource},
    source::PSY.Source,
    ::AbstractVector,
)
    device_states = zeros(PSY.get_n_states(dynamic_device))

    #PowerFlow Data
    P0 = PSY.get_active_power(source)
    Q0 = PSY.get_reactive_power(source)
    Vm = PSY.get_magnitude(PSY.get_bus(source))
    θ = PSY.get_angle(PSY.get_bus(source))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)
    I_R = real(I)
    I_I = imag(I)
    R_th = PSY.get_R_th(source)
    X_th = PSY.get_X_th(source)
    Zmag = R_th^2 + X_th^2
    function f!(out, x)
        V_R_internal = x[1]
        V_I_internal = x[2]

        out[1] =
            R_th * (V_R_internal - V_R) / Zmag + X_th * (V_I_internal - V_I) / Zmag - I_R
        out[2] =
            R_th * (V_I_internal - V_I) / Zmag - X_th * (V_R_internal - V_R) / Zmag - I_I
    end
    x0 = [V_R, V_I]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Chirp Signal Source failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        V_internal = sqrt(sol_x0[1]^2 + sol_x0[2]^2)
        θ_internal = atan(sol_x0[2], sol_x0[1])

        device_states[1] = V_internal
        device_states[2] = θ_internal
    end
    return device_states
end

PSID.get_inner_vars_count(::ChirpVariableSource) = 0
