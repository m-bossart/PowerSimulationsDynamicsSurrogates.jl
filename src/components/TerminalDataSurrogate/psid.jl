
function PSID.StaticWrapper(device::TerminalDataSurrogate, bus_ix::Int)
    bus = PSY.get_bus(device)
    ext_wrapper = Dict{String, Any}()
    return PSID.StaticWrapper{TerminalDataSurrogate, PSID.BUS_MAP[PSY.get_bustype(bus)]}(
        device,
        Base.Ref(1.0),
        Base.Ref(PSY.get_internal_voltage(device)),
        Base.Ref(PSY.get_internal_angle(device)),
        Base.Ref(PSY.get_active_power(device)),
        Base.Ref(PSY.get_reactive_power(device)),
        bus_ix,
        ext_wrapper,
    )
end

add_dim(x) = reshape(x, (size(x)..., 1))

function PSID.initialize_static_device!(
    device::PSID.StaticWrapper{TerminalDataSurrogate, T},
) where {T <: PSID.BusCategory}
    ext_wrapper = PSID.get_ext(device)
    ext_device = PSY.get_ext(device)

    P0 = PSY.get_active_power(device)
    Q0 = PSY.get_reactive_power(device)
    Vm = PSY.get_magnitude(PSY.get_bus(device))
    θ = PSY.get_angle(PSY.get_bus(device))
    S0 = P0 + Q0 * 1im

    VR0 = Vm * cos(θ)
    VI0 = Vm * sin(θ)
    V = VR0 + VI0 * 1im
    I = conj(S0 / V)
    IR0 = real(I)
    II0 = imag(I)

    ext_wrapper["v0"] = [VR0, VI0]
    ext_wrapper["i0"] = [IR0, II0]

    model = ext_device["model"]
    ps = ext_device["ps"]
    st = ext_device["st"]
    window_size = get_window_size(device.device)
    v_ss = add_dim(hcat(fill(VR0, window_size), fill(VI0, window_size)))
    i_ss = add_dim(hcat(fill(IR0, window_size), fill(II0, window_size)))
    ss_input = ([VR0, VI0], [IR0, II0], v_ss, i_ss)
    y, st = model(ss_input, ps, st)
    ext_wrapper["offset"] = y - [IR0, II0]

    #TODO -add a warning if offset is too large; offset should be zero if well trained at steady state operation
    PSID.set_V_ref(device, Vm)
    PSID.set_θ_ref(device, θ)
    return
end

function PSID.device!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    global_vars::AbstractArray{T},
    ::AbstractArray{T},
    device::PSID.StaticWrapper{TerminalDataSurrogate, U},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES, U <: PSID.BusCategory}
    mdl_solution_prediction_surrogate!(
        voltage_r,
        voltage_i,
        current_r,
        current_i,
        global_vars,
        device,
        t,
    )
    return
end

function mdl_solution_prediction_surrogate!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    global_vars, #::AbstractArray{T},
    static_device::PSID.StaticWrapper{TerminalDataSurrogate},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES}
    ext_wrapper = PSID.get_ext(static_device)
    ext_device = PSY.get_ext(static_device.device)
    window_size = get_window_size(static_device.device)
    τ = get_τ(static_device.device)
    v0 = ext_wrapper["v0"]
    i0 = ext_wrapper["i0"]
    offset = ext_wrapper["offset"]
    model = ext_device["model"]
    ps = ext_device["ps"]
    st = ext_device["st"]

    #TODO - build v, i based on window_size and τ and access to history function
    v_ss = add_dim(hcat(fill(v0[1], window_size), fill(v0[2], window_size)))
    i_ss = add_dim(hcat(fill(i0[1], window_size), fill(i0[2], window_size)))
    x = (v0, i0, v_ss, i_ss)
    y_pred, st = model(x, ps, st)

    current_r[1] += y_pred[1] - offset[1]
    current_i[1] += y_pred[2] - offset[2]

    return
end
