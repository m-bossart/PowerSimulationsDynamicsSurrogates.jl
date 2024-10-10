
mutable struct CurrentPlayback <: PSY.StaticInjection
    name::String
    available::Bool
    bus::PSY.ACBus
    active_power::Float64
    reacive_power::Float64
    base_power::Float64
    playback_name::String
    playback_type::Symbol   #Source, DynamicInjection, Branch -> impacts how you get the currents
    playback_result::Union{Nothing, PSID.SimulationResults}
    reverse_current_polarity::Bool
    dynamic_injector::Union{Nothing, PSY.DynamicInjection}
    services::Vector{PSY.Service}
    ext::Dict{String, Any}
    internal::PSY.InfrastructureSystemsInternal
end

function CurrentPlayback(
    name,
    available,
    bus,
    active_power,
    reactive_power,
    base_power,
    playback_name,
    playback_type,
    playback_result,
    reverse_current_polarity,
    dynamic_injector = nothing,
    services = PSY.Device[],
    ext = Dict{String, Any}(),
)
    CurrentPlayback(
        name,
        available,
        bus,
        active_power,
        reactive_power,
        base_power,
        playback_name,
        playback_type,
        playback_result,
        reverse_current_polarity,
        dynamic_injector,
        services,
        ext,
        PSY.InfrastructureSystemsInternal(),
    )
end

function CurrentPlayback(;
    name,
    available,
    bus,
    active_power,
    reactive_power,
    base_power,
    playback_name,
    playback_type,
    playback_result,
    reverse_current_polarity,
    dynamic_injector = nothing,
    services = PSY.Device[],
    ext = Dict{String, Any}(),
    internal = PSY.InfrastructureSystemsInternal(),
)
    CurrentPlayback(
        name,
        available,
        bus,
        active_power,
        reactive_power,
        base_power,
        playback_name,
        playback_type,
        playback_result,
        reverse_current_polarity,
        dynamic_injector,
        services,
        ext,
        internal,
    )
end

# Constructor for demo purposes; non-functional.
function CurrentPlayback(::Nothing)
    CurrentPlayback(;
        name = "init",
        available = false,
        bus = PSY.ACBus(nothing),
        active_power = 0.0,
        reactive_power = 0.0,
        base_power = 0.0,
        playback_name = "",
        playback_type = :none,
        playback_result = nothing,
        reverse_current_polarity = false,
        dynamic_injector = nothing,
        services = PSY.Device[],
        ext = Dict{String, Any}(),
    )
end

PSY.get_name(value::CurrentPlayback) = value.name
PSY.get_available(value::CurrentPlayback) = value.available
PSY.get_bus(value::CurrentPlayback) = value.bus
PSY.get_active_power(value::CurrentPlayback) = value.active_power
PSY.get_reactive_power(value::CurrentPlayback) = value.reacive_power
PSY.get_base_power(value::CurrentPlayback) = value.base_power
PSY.get_states(value::CurrentPlayback) = value.states
PSY.get_n_states(value::CurrentPlayback) = value.n_states
PSY.get_ext(value::CurrentPlayback) = value.ext
PSY.get_internal(value::CurrentPlayback) = value.internal
get_playback_name(value::CurrentPlayback) = value.playback_name
get_playback_type(value::CurrentPlayback) = value.playback_type
get_playback_result(value::CurrentPlayback) = value.playback_result
get_reverse_current_polarity(value::CurrentPlayback) = value.reverse_current_polarity

PSY.set_active_power!(value::CurrentPlayback, val) = value.active_power = val
PSY.set_reactive_power!(value::CurrentPlayback, val) = value.reacive_power = val

PSY.set_base_power!(value::CurrentPlayback, val) = value.base_power = val
PSY.set_ext!(value::CurrentPlayback, val) = value.ext = val

PSY.get_reactive_power_limits(::CurrentPlayback) = (min = -Inf, max = Inf)
