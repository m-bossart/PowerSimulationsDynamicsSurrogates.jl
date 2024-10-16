PSID.get_delays(dynamic_injector::TerminalDataSurrogate) =
    [x * get_τ(dynamic_injector) for x in 1:(get_window_size(dynamic_injector) + 1)]
PSID.is_valid(::TerminalDataSurrogate) = nothing
PSID.get_inner_vars_count(::TerminalDataSurrogate) = 0
PSID._get_frequency_state(d::PSID.DynamicWrapper{TerminalDataSurrogate}) = 0

function PSID._get_refs(x::PSID.DynamicWrapper{TerminalDataSurrogate})
    return (ir_offset = 0.0, ii_offset = 0.0)
end

function PSID._get_refs_metadata(x::PSID.DynamicWrapper{TerminalDataSurrogate})
    return (
        ir_offset = PSID.ParamsMetadata(PSID.DEVICE_SETPOINT, false, true),
        ii_offset = PSID.ParamsMetadata(PSID.DEVICE_SETPOINT, false, true),
    )
end

PSID.get_params(x::TerminalDataSurrogate) = (; θ = PSY.get_ext(x)["ps"])
PSID.get_params_metadata(::TerminalDataSurrogate) =
    (; θ = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true))

function PSID.device_mass_matrix_entries!(
    mass_matrix::AbstractArray,
    dynamic_device::PSID.DynamicWrapper{TerminalDataSurrogate},
)
    global_index = PSID.get_global_index(dynamic_device)
    mass_matrix_entries!(mass_matrix, dynamic_device, global_index)
    return
end

function mass_matrix_entries!(
    mass_matrix,
    pvs::PSID.DynamicWrapper{TerminalDataSurrogate},
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    fc = get_fc(pvs.device)
    mass_matrix[global_index[:ir], global_index[:ir]] = fc
    mass_matrix[global_index[:ii], global_index[:ii]] = fc
    mass_matrix[global_index[:vr], global_index[:vr]] = fc
    mass_matrix[global_index[:vi], global_index[:vi]] = fc
end

function PSID.DynamicWrapper(
    static_device::PSY.Source,
    dynamic_device::TerminalDataSurrogate,
    bus_ix::Int,
    ix_range,
    ode_range,
    inner_var_range,
    sys_base_power,
    sys_base_freq,
)
    device_states = PSY.get_states(dynamic_device)

    return PSID.DynamicWrapper(
        dynamic_device,
        sys_base_power,
        sys_base_freq,
        static_device,
        PSID.BUS_MAP[PSY.get_bustype(PSY.get_bus(static_device))],
        1.0,
        collect(inner_var_range),
        collect(ix_range),
        collect(ode_range),
        bus_ix,
        Base.ImmutableDict(Dict(device_states .=> ix_range)...),
        Base.ImmutableDict{Int, Vector{Int}}(),
        Base.ImmutableDict{Int, Vector{Int}}(),
        Dict{String, Any}(),
    )
end

add_dim(x) = reshape(x, (size(x)..., 1))

function PSID.initialize_dynamic_device!(
    dynamic_device::PSID.DynamicWrapper{TerminalDataSurrogate},
    source::PSY.Source,
    ::AbstractVector,
    p::AbstractVector,
    device_states::AbstractVector,
)
    ext_wrapper = PSID.get_ext(dynamic_device)
    device = PSID.get_device(dynamic_device)
    ext_device = PSY.get_ext(device)

    basepower = PSY.get_base_power(dynamic_device)
    sys_Sbase = PSID.get_system_base_power(dynamic_device)
    P0 = PSY.get_active_power(source) * (sys_Sbase / basepower)
    Q0 = PSY.get_reactive_power(source) * (sys_Sbase / basepower)
    Vm = PSY.get_magnitude(PSY.get_bus(source))
    θ = PSY.get_angle(PSY.get_bus(source))
    S0 = P0 + Q0 * 1im

    VR0 = Vm * cos(θ)
    VI0 = Vm * sin(θ)
    V = VR0 + VI0 * 1im
    I = conj(S0 / V)
    IR0 = real(I)
    II0 = imag(I)

    ext_wrapper["v0"] = [VR0, VI0]
    ext_wrapper["i0"] = [IR0, II0]
    ext_wrapper["voltage_violation_warning"] = false

    model = ext_device["model"]
    ps = p[:params][:θ]
    st = ext_device["st"]
    window_size = get_window_size(device)
    v_ss = add_dim(vcat(fill(VR0, (1, window_size)), fill(VI0, (1, window_size))))
    i_ss = add_dim(vcat(fill(IR0, (1, window_size)), fill(II0, (1, window_size))))
    @assert size(v_ss)[1] == 2
    ss_input = (add_dim([VR0, VI0]), add_dim([IR0, II0]), v_ss, i_ss)
    y, st = model(ss_input, ps, st)
    if get_steadystate_offset_correction(device)
        device_states[1:4] = [IR0, II0, VR0, VI0]
        @warn "The surrogate model has a non-zero error at time zero which is corrected with an offset. Ir: $(y[1] - IR0), Ii: $(y[2] - II0)"
        @view(p[:refs])[:ir_offset] = y[1] - IR0
        @view(p[:refs])[:ii_offset] = y[2] - II0
    else
        @view(p[:refs])[:ir_offset] = 0.0
        @view(p[:refs])[:ii_offset] = 0.0
        device_states[1:4] = [y[1], y[2], VR0, VI0]
    end
    return
end

function PSID.device!(
    device_states::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{T},
    p::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    ::AbstractArray{T},
    ::AbstractArray{T},
    dynamic_device::PSID.DynamicWrapper{TerminalDataSurrogate},
    h,
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES}
    ext_wrapper = PSID.get_ext(dynamic_device)
    ext_device = PSY.get_ext(dynamic_device.device)
    window_size = get_window_size(dynamic_device.device)
    trained_voltage_range = get_trained_voltage_range(dynamic_device.device)
    τ = get_τ(dynamic_device.device)
    model = ext_device["model"]

    offset = [p[:refs][:ir_offset], p[:refs][:ii_offset]]
    ps = p[:params][:θ]
    st = ext_device["st"]
    vm = sqrt(voltage_r^2 + voltage_i^2)
    if ((vm < trained_voltage_range[1]) || (vm > trained_voltage_range[2])) &&
       t > 0.0 &&
       !ext_wrapper["voltage_violation_warning"]
        ext_wrapper["voltage_violation_warning"] = true
    end

    ix_range = dynamic_device.ix_range
    ir_ix = ix_range[1]
    ii_ix = ix_range[2]
    vr_ix = ix_range[3]
    vi_ix = ix_range[4]

    ir = device_states[1]
    ii = device_states[2]
    vr = device_states[3]
    vi = device_states[4]
    v0 = Array{typeof(ext_wrapper["v0"][1])}(undef, (2, 1))
    v0 .= ext_wrapper["v0"]
    i0 = Array{typeof(ext_wrapper["i0"][1])}(undef, (2, 1))
    i0 .= ext_wrapper["i0"]
    v = Array{typeof(vr)}(undef, (2, window_size, 1))
    i = Array{typeof(vr)}(undef, (2, window_size, 1))
    v[1, :, 1] =
        vcat([h(nothing, t - N * τ; idxs = vr_ix) for N in (window_size - 1):-1:1], [vr])
    v[2, :, 1] =
        vcat([h(nothing, t - N * τ; idxs = vi_ix) for N in (window_size - 1):-1:1], [vi])
    i[1, :, 1] = [h(nothing, t - N * τ; idxs = ir_ix) for N in window_size:-1:1]
    i[2, :, 1] = [h(nothing, t - N * τ; idxs = ii_ix) for N in window_size:-1:1]

    x = (v0, i0, v, i)
    y_pred, st = model(x, ps, st)

    output_ode[1] = y_pred[1] - offset[1] - ir
    output_ode[2] = y_pred[2] - offset[2] - ii
    output_ode[3] = voltage_r - vr
    output_ode[4] = voltage_i - vi
    basepower = PSY.get_base_power(dynamic_device)
    sys_Sbase = PSID.get_system_base_power(dynamic_device)
    current_r[1] += (y_pred[1] - offset[1]) * (basepower / sys_Sbase)
    current_i[1] += (y_pred[2] - offset[2]) * (basepower / sys_Sbase)
    return
end
