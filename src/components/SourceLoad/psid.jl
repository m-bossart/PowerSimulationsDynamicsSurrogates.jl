
function PSID.StaticWrapper(device::T, bus_ix::Int, p_range) where {T <: SourceLoad}
    bus = PSY.get_bus(device)
    return PSID.StaticWrapper{T, PSID.BUS_MAP[PSY.get_bustype(bus)]}(
        device,
        Base.Ref(1.0),
        Base.Ref(PSY.get_magnitude(bus)),
        Base.Ref(PSY.get_angle(bus)),
        Base.Ref(PSY.get_active_power(device)), #PSY.get_active_power(device))
        Base.Ref(PSY.get_reactive_power(device)), #PSY.get_reactive_power(device)
        p_range,
        bus_ix,
        Dict{String, Any}(),
    )
end

function PSID.initialize_static_device!(
    device::PSID.StaticWrapper{SourceLoad, T},
    local_parameters::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
) where {T <: PSID.BusCategory}
    return
end

function PSID.device!(
    device_parameters::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    ::AbstractArray{T},
    ::AbstractArray{T},
    device::PSID.StaticWrapper{SourceLoad, U},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES, U <: PSID.BusCategory}
    mdl_sourceload!(voltage_r, voltage_i, current_r, current_i, device)
    return
end

function mdl_sourceload!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    wrapper::PSID.StaticWrapper{SourceLoad},
) where {T <: PSID.ACCEPTED_REAL_TYPES}
    # Read power flow voltages
    #V0_mag_inv = 1.0 / get_V_ref(wrapper)
    V0_mag_inv = 1.0 / PSID.get_V_ref(wrapper)
    V0_mag_sq_inv = V0_mag_inv^2

    # Load device parameters
    P_impedance = PSY.get_active_power(wrapper)
    Q_impedance = PSY.get_reactive_power(wrapper)

    # Compute ZIP currents
    Iz_re = V0_mag_sq_inv * (voltage_r * P_impedance + voltage_i * Q_impedance)
    Iz_im = V0_mag_sq_inv * (voltage_i * P_impedance - voltage_r * Q_impedance)

    # Update current
    current_r[1] += Iz_re  #in system pu flowing out
    current_i[1] += Iz_im  #in system pu flowing out

    return
end
