PSID.get_delays(dynamic_injector::IdealTerminalDataSurrogate) = [get_τ(dynamic_injector)]
PSID.is_valid(::IdealTerminalDataSurrogate) = nothing
PSID.get_inner_vars_count(::IdealTerminalDataSurrogate) = 0
PSID._get_frequency_state(d::PSID.DynamicWrapper{IdealTerminalDataSurrogate}) = 0

function PSID._get_refs(x::PSID.DynamicWrapper{IdealTerminalDataSurrogate})
    return (;)
end

function PSID._get_refs_metadata(x::PSID.DynamicWrapper{IdealTerminalDataSurrogate})
    return (;)
end

function PSID.device_mass_matrix_entries!(
    mass_matrix::AbstractArray,
    dynamic_device::PSID.DynamicWrapper{IdealTerminalDataSurrogate},
)
    global_index = PSID.get_global_index(dynamic_device)
    mass_matrix_entries!(mass_matrix, dynamic_device, global_index)

    return
end

function mass_matrix_entries!(
    mass_matrix,
    pvs::PSID.DynamicWrapper{IdealTerminalDataSurrogate},
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:ir], global_index[:ir]] = 0.0
    mass_matrix[global_index[:ii], global_index[:ii]] = 0.0
end

PSID.get_params(x::IdealTerminalDataSurrogate) = (;)
PSID.get_params_metadata(::IdealTerminalDataSurrogate) = (;)

function PSID.DynamicWrapper(
    static_device::PSY.Source,
    dynamic_device::IdealTerminalDataSurrogate,
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

function get_value_from_reference_solution(
    t,
    reference_solution,
    reference_solution_state_indices,
)
    if t < 0
        t_evaluate = 0.0
    else
        t_evaluate = t
    end
    value = 0.0
    for ix in reference_solution_state_indices
        value += reference_solution(t_evaluate, idxs = ix)
    end
    return value
end

function PSID.initialize_dynamic_device!(
    dynamic_device::PSID.DynamicWrapper{IdealTerminalDataSurrogate},
    source::PSY.Source,
    ::AbstractVector,
    p::AbstractVector,
    device_states::AbstractVector,
)
    device = PSID.get_device(dynamic_device)
    reference_solution = get_reference_solution(device)
    reference_solution_ir_indices = get_reference_solution_ir_indices(device)
    reference_solution_ii_indices = get_reference_solution_ii_indices(device)
    ir0 = get_value_from_reference_solution(
        0.0,
        reference_solution,
        reference_solution_ir_indices,
    )
    ii0 = get_value_from_reference_solution(
        0.0,
        reference_solution,
        reference_solution_ii_indices,
    )
    device_states[1:2] = [ir0, ii0]
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
    dynamic_device::PSID.DynamicWrapper{IdealTerminalDataSurrogate},
    h,
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES}
    device = PSID.get_device(dynamic_device)
    ext_device = PSY.get_ext(device)
    τ = get_τ(dynamic_device.device)

    if !isempty(ext_device)
        ix_range = dynamic_device.ix_range
        ir_ix = ix_range[1]
        ii_ix = ix_range[2]
        ir_delayed = h(nothing, t - τ; idxs = ir_ix)
        ii_delayed = h(nothing, t - τ; idxs = ii_ix)

        model = ext_device["model"]
        ps = ext_device["ps"]
        st = ext_device["st"]
        y, _ = model([ir_delayed, ii_delayed], ps, st)
    end
    reference_solution = get_reference_solution(dynamic_device.device)
    reference_solution_ir_indices = get_reference_solution_ir_indices(dynamic_device.device)
    reference_solution_ii_indices = get_reference_solution_ii_indices(dynamic_device.device)

    ix_range = dynamic_device.ix_range
    ir_ix = ix_range[1]
    ii_ix = ix_range[2]
    ir_delayed = h(nothing, t - τ; idxs = ir_ix)
    ii_delayed = h(nothing, t - τ; idxs = ii_ix)
    ir_ref_t = get_value_from_reference_solution(
        t,
        reference_solution,
        reference_solution_ir_indices,
    )
    ir_ref_τ = get_value_from_reference_solution(
        t - τ,
        reference_solution,
        reference_solution_ir_indices,
    )
    ii_ref_t = get_value_from_reference_solution(
        t,
        reference_solution,
        reference_solution_ii_indices,
    )
    ii_ref_τ = get_value_from_reference_solution(
        t - τ,
        reference_solution,
        reference_solution_ii_indices,
    )
    ir_pred = ir_ref_t + ir_delayed - ir_ref_τ
    ii_pred = ii_ref_t + ii_delayed - ii_ref_τ

    ir = device_states[1]
    ii = device_states[2]
    output_ode[1] = ir_pred - ir
    output_ode[2] = ii_pred - ii
    basepower = PSY.get_base_power(dynamic_device)
    sys_Sbase = PSID.get_system_base_power(dynamic_device)
    current_r[1] += ir_pred * (basepower / sys_Sbase)
    current_i[1] += ii_pred * (basepower / sys_Sbase)
    return
end
