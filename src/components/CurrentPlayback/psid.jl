
function PSID.StaticWrapper(device::T, bus_ix::Int) where {T <: CurrentPlayback}
    bus = PSY.get_bus(device)
    return PSID.StaticWrapper{T, PSID.BUS_MAP[PSY.get_bustype(bus)]}(
        device,
        1.0,
        bus_ix,
        Dict{String, Any}(),
    )
end

function PSID.initialize_static_device!(
    device::PSID.StaticWrapper{CurrentPlayback, T},
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
    device::PSID.StaticWrapper{CurrentPlayback, U},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES, U <: PSID.BusCategory}
    mdl_CurrentPlayback!(voltage_r, voltage_i, current_r, current_i, device, t)
    return
end

function mdl_CurrentPlayback!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    wrapper::PSID.StaticWrapper{CurrentPlayback},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES}
    device = PSID.get_device(wrapper)
    playback_name = get_playback_name(device)
    playback_type = get_playback_type(device)
    playback_result = get_playback_result(device)
    reverse_current_polarity = get_reverse_current_polarity(device)

    if t == 0.0
        if playback_type == :DynamicInjection
            ir = PSID.get_real_current_series(playback_result, playback_name)[2][1]
            ii = PSID.get_imaginary_current_series(playback_result, playback_name)[2][1]
        elseif playback_type == :Source
            ir = PSID.get_source_real_current_series(playback_result, playback_name)[2][1]
            ii =
                PSID.get_source_imaginary_current_series(playback_result, playback_name)[2][1]
        elseif playback_type == :Branch
            ir = PSID.get_real_current_branch_flow(playback_result, playback_name)[2][1]
            ii =
                PSID.get_imaginary_current_branch_flow(playback_result, playback_name)[2][1]
        else
            @error "Invalid playback_type, must be :DynamicInjection, :Source, or :Branch"
        end
    else
        if playback_type == :DynamicInjection
            ir = PSID.get_real_current_series(playback_result, playback_name, t)[2][1]
            ii = PSID.get_imaginary_current_series(playback_result, playback_name, t)[2][1]
        elseif playback_type == :Source
            ir =
                PSID.get_source_real_current_series(playback_result, playback_name, t)[2][1]
            ii = PSID.get_source_imaginary_current_series(
                playback_result,
                playback_name,
                t,
            )[2][1]
        elseif playback_type == :Branch
            ir =
                PSID.get_real_current_branch_flow(playback_result, playback_name; dt = t)[2][1]
            ii = PSID.get_imaginary_current_branch_flow(
                playback_result,
                playback_name;
                dt = t,
            )[2][1]
        else
            @error "Invalid playback_type, must be :DynamicInjection, :Source, or :Branch"
        end
    end

    current_r[1] += ir #in system pu flowing out
    current_i[1] += ii #in system pu flowing out

    return
end
