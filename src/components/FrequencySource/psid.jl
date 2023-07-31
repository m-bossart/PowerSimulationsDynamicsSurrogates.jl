PSID.is_valid(::FrequencySource) = nothing

function PSID._get_frequency_state(d::PSID.DynamicWrapper{FrequencySource})
    return PSID.get_global_index(d)[:ω]
end

function PSID.device_mass_matrix_entries!(
    mass_matrix::PSID.AbstractArray,
    dynamic_device::PSID.DynamicWrapper{FrequencySource},
)
    global_index = PSID.get_global_index(dynamic_device)
    mass_matrix_device_entries!(mass_matrix, dynamic_device, global_index)
    return
end

function mass_matrix_device_entries!(
    mass_matrix,
    device::PSID.DynamicWrapper{FrequencySource},
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:Vt], global_index[:Vt]] = PSY.get_Tv(PSID.get_device(device))
    mass_matrix[global_index[:ω], global_index[:ω]] = get_Tω(PSID.get_device(device))
    mass_matrix[global_index[:θt], global_index[:θt]] = 1 / (2 * pi * 60)
    return
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
    dynamic_device::PSID.DynamicWrapper{FrequencySource},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES}

    # Internal Voltage states
    Vt = device_states[1]
    θt = device_states[2]
    ω = device_states[3]

    V_R = Vt * cos(θt)
    V_I = Vt * sin(θt)

    V_ref = PSID.get_V_ref(dynamic_device)
    ω_ref = PSID.get_ω_ref(dynamic_device)

    Tv = PSY.get_Tv(PSID.get_device(dynamic_device))
    Tω = get_Tω(PSID.get_device(dynamic_device))

    _, dVt_dt = PSID.low_pass_mass_matrix(V_ref, Vt, 1.0, Tv)
    _, dω_dt = PSID.low_pass_mass_matrix(ω_ref, ω, 1.0, Tω)
    _, dθt_dt = PSID.integrator_nonwindup_mass_matrix(
        ω - 1.0,
        θt,
        1.0,
        1 / (2 * pi * 60),
        -Inf,
        Inf,
    )

    output_ode[1] = dVt_dt
    output_ode[2] = dθt_dt
    output_ode[3] = dω_dt

    #update current
    R_th = PSY.get_R_th(PSID.get_device(dynamic_device))
    X_th = PSY.get_X_th(PSID.get_device(dynamic_device))
    Zmag = R_th^2 + X_th^2
    current_r[1] += R_th * (V_R - voltage_r[1]) / Zmag + X_th * (V_I - voltage_i[1]) / Zmag #in system pu flowing out
    current_i[1] += R_th * (V_I - voltage_i[1]) / Zmag - X_th * (V_R - voltage_r[1]) / Zmag #in system pu flowing out

    return
end

function PSID.initialize_dynamic_device!(
    dynamic_device::PSID.DynamicWrapper{FrequencySource},
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

        PSID.set_V_ref(dynamic_device, V_internal)
        PSID.set_ω_ref(dynamic_device, 1.0)
    end
    return device_states
end

PSID.get_inner_vars_count(::FrequencySource) = 0

"""
Function to obtain the output current time series of a PeriodicVariableSource model out of the DAE Solution. It receives the simulation inputs,
the dynamic device and bus voltage. It is dispatched for device type to compute the specific current.
compute_output_current(::SimulationResults, ::PeriodicVariableSource, ::Vector{Float64}, ::Vector{Float64}, ::Nothing)
"""
function PSID.compute_output_current(
    res::PSID.SimulationResults,
    dynamic_device::FrequencySource,
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
