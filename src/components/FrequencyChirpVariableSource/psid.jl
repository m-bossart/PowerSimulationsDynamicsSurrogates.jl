PSID.is_valid(::FrequencyChirpVariableSource) = nothing

function PSID._get_frequency_state(d::PSID.DynamicWrapper{FrequencyChirpVariableSource})
    return PSID.get_global_index(d)[:ω]
end

function PSID.device_mass_matrix_entries!(
    mass_matrix::PSID.AbstractArray,
    dynamic_device::PSID.DynamicWrapper{FrequencyChirpVariableSource},
)
    global_index = PSID.get_global_index(dynamic_device)
    mass_matrix_chirp_entries!(mass_matrix, dynamic_device, global_index)
    return
end

function mass_matrix_chirp_entries!(
    mass_matrix,
    chirp::PSID.DynamicWrapper{FrequencyChirpVariableSource},
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
    dynamic_device::PSID.DynamicWrapper{FrequencyChirpVariableSource},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES}

    # Internal Voltage states
    V_R = device_states[1] * cos(device_states[2])
    V_I = device_states[1] * sin(device_states[2])
    ω = device_states[3]

    tstart = get_tstart(PSID.get_device(dynamic_device))
    N = get_N(PSID.get_device(dynamic_device))
    ω1 = get_ω1(PSID.get_device(dynamic_device))
    ω2 = get_ω2(PSID.get_device(dynamic_device))
    V_amp = get_V_amp(PSID.get_device(dynamic_device))
    ω_amp = get_ω_amp(PSID.get_device(dynamic_device))

    if t >= tstart && t < N
        output_ode[1] =
            V_amp *
            cos(ω1 * (t - tstart) + (ω2 - ω1) * (t - tstart)^2 / (2 * N)) *
            (ω1 + (ω2 - ω1) * (t - tstart) / N)
        output_ode[2] = ω - 1.0
        output_ode[3] =
            ω_amp *
            cos(ω1 * (t - tstart) + (ω2 - ω1) * (t - tstart)^2 / (2 * N)) *
            (ω1 + (ω2 - ω1) * (t - tstart) / N)
    else
        output_ode[1] = 0
        output_ode[2] = 0
        output_ode[3] = 0
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
    dynamic_device::PSID.DynamicWrapper{FrequencyChirpVariableSource},
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
        device_states[3] = 1.0
    end
    return device_states
end

PSID.get_inner_vars_count(::FrequencyChirpVariableSource) = 0

"""
Function to obtain the output current time series of a PeriodicVariableSource model out of the DAE Solution. It receives the simulation inputs,
the dynamic device and bus voltage. It is dispatched for device type to compute the specific current.
compute_output_current(::SimulationResults, ::PeriodicVariableSource, ::Vector{Float64}, ::Vector{Float64}, ::Nothing)
"""
function PSID.compute_output_current(
    res::PSID.SimulationResults,
    dynamic_device::FrequencyChirpVariableSource,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    dt::Union{Nothing, Float64},
)
    name = PSY.get_name(dynamic_device)
    ts, Vt_internal = PSID.post_proc_state_series(res, (name, :Vt), dt)
    _, θt_internal = PSID.post_proc_state_series(res, (name, :θt), dt)
    Vr_internal = Vt_internal .* cos.(θt_internal)
    Vi_internal = Vt_internal .* sin.(θt_internal)
    R_th = PSY.get_R_th(dynamic_device)
    X_th = PSY.get_X_th(dynamic_device)
    Z_sq = R_th^2 + X_th^2
    I_R = R_th * (Vr_internal .- V_R) / Z_sq + X_th * (Vi_internal .- V_I) / Z_sq
    I_I = R_th * (Vi_internal .- V_I) / Z_sq - X_th * (Vr_internal .- V_R) / Z_sq

    return ts, I_R, I_I
end
