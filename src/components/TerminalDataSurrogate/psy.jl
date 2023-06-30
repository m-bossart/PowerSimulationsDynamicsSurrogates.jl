"""
    mutable struct TerminalDataSurrogate{
        M <: DataDrivenModelArchitecture,
        D <: PSY.DynamicInjection,
    } <: LearnedSolutionSurrogate
        name::String
        available::Bool
        bus::PSY.Bus
        active_power::Float64
        reactive_power::Float64
        active_power_limits::PSY.MinMax
        reactive_power_limits::Union{Nothing, PSY.MinMax}
        internal_voltage::Float64  
        internal_angle::Float64     
        model_architecture::M      
        underlying_dynamic_model::D 
        n_past_timesteps::Int64
        services::Vector{PSY.Service}
        ext::Dict{String, Any}
        internal::IS.InfrastructureSystemsInternal
    end

Experimental surrogate

# Arguments
- `name::String`
- `available::Bool`: 
- `bus::Bus`: 
- `active_power::Float64`: 
- `reactive_power::Float64`:
- `active_power_limits::PSY.MinMax`:
- `reactive_power_limits::Union{Nothing, PSY.MinMax}`:
- `internal_voltage::Float64`:
- `internal_angle::Float64`:
- `model_architecture::M` :     
- `underlying_dynamic_model::D`: 
- `n_past_timesteps::Int64`: 
- `services::Vector{Service}`
- `ext::Dict{String, Any}`
- `internal::InfrastructureSystemsInternal`: power system internal reference, do not modify
"""
mutable struct TerminalDataSurrogate{
    M <: DataDrivenModelArchitecture,
    D <: PSY.DynamicInjection,
} <: LearnedSolutionSurrogate
    name::String
    available::Bool
    bus::PSY.Bus
    active_power::Float64
    reactive_power::Float64
    active_power_limits::PSY.MinMax
    reactive_power_limits::Union{Nothing, PSY.MinMax}
    internal_voltage::Float64
    internal_angle::Float64
    model_architecture::M
    underlying_dynamic_model::D
    θ_ref_frame::Float64
    n_past_timesteps::Int64
    services::Vector{PSY.Service}
    ext::Dict{String, Any}
    internal::IS.InfrastructureSystemsInternal
end

function TerminalDataSurrogate(
    name,
    available,
    bus,
    active_power,
    reactive_power,
    active_power_limits,
    reactive_power_limits,
    internal_voltage,
    internal_angle,
    model_architecture,
    underlying_dynamic_model,
    θ_ref_frame,
    n_past_timesteps,
    services = PSY.Service[],
    ext = Dict{String, Any}(),
)
    TerminalDataSurrogate(
        name,
        available,
        bus,
        active_power,
        reactive_power,
        active_power_limits,
        reactive_power_limits,
        internal_voltage,
        internal_angle,
        model_architecture,
        underlying_dynamic_model,
        θ_ref_frame,
        n_past_timesteps,
        services,
        ext,
        IS.InfrastructureSystemsInternal(),
    )
end

function TerminalDataSurrogate(;
    name,
    available,
    bus,
    active_power,
    reactive_power,
    active_power_limits,
    reactive_power_limits,
    internal_voltage,
    internal_angle,
    model_architecture,
    underlying_dynamic_model,
    θ_ref_frame,
    n_past_timesteps,
    services = PSY.Service[],
    ext = Dict{String, Any}(),
)
    TerminalDataSurrogate(
        name,
        available,
        bus,
        active_power,
        reactive_power,
        active_power_limits,
        reactive_power_limits,
        internal_voltage,
        internal_angle,
        model_architecture,
        underlying_dynamic_model,
        θ_ref_frame,
        n_past_timesteps,
        services,
        ext,
        IS.InfrastructureSystemsInternal(),
    )
end

# Constructor for demo purposes; non-functional.
function TerminalDataSurrogate(::Nothing)
    TerminalDataSurrogate(;
        name = "init",
        available = false,
        bus = PSY.Bus(nothing),
        active_power = 0.0,
        reactive_power = 0.0,
        active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = 0.0, max = 1.0),
        internal_voltage = 0.0,
        internal_angle = 0.0,
        model_architecture = FullyConnected(nothing),
        underlying_dynamic_model = PSY.GenericDER(nothing),
        θ_ref_frame = 0.0,
        n_past_timesteps = 0,
        services = PSY.Service[],
        ext = Dict{String, Any}(),
    )
end

#If function exists in PSY, overload it here. 
"""Get [`TerminalDataSurrogate`](@ref) `name`."""
PSY.get_name(value::TerminalDataSurrogate) = value.name
"""Get [`TerminalDataSurrogate`](@ref) `available`."""
PSY.get_available(value::TerminalDataSurrogate) = value.available
"""Get [`TerminalDataSurrogate`](@ref) `bus`."""
PSY.get_bus(value::TerminalDataSurrogate) = value.bus
"""Get [`TerminalDataSurrogate`](@ref) `active_power`."""
PSY.get_active_power(value::TerminalDataSurrogate) = value.active_power
"""Get [`TerminalDataSurrogate`](@ref) `reactive_power`."""
PSY.get_reactive_power(value::TerminalDataSurrogate) = value.reactive_power
"""Get [`TerminalDataSurrogate`](@ref) `active_power_limits`."""
PSY.get_active_power_limits(value::TerminalDataSurrogate) = value.active_power_limits
"""Get [`TerminalDataSurrogate`](@ref) `reactive_power_limits`."""
PSY.get_reactive_power_limits(value::TerminalDataSurrogate) = value.reactive_power_limits
"""Get [`TerminalDataSurrogate`](@ref) `internal_voltage`."""
PSY.get_internal_voltage(value::TerminalDataSurrogate) = value.internal_voltage
"""Get [`TerminalDataSurrogate`](@ref) `internal_angle`."""
PSY.get_internal_angle(value::TerminalDataSurrogate) = value.internal_angle
"""Get [`FullyConnected`](@ref) `model_architecture`."""
get_model_architecture(value::TerminalDataSurrogate) = value.model_architecture
"""Get [`FullyConnected`](@ref) `underlying_dynamic_model`."""
get_underlying_dynamic_model(value::TerminalDataSurrogate) = value.underlying_dynamic_model
"""Get [`FullyConnected`](@ref) `θ_ref_frame`."""
get_θ_ref_frame(value::TerminalDataSurrogate) = value.θ_ref_frame
"""Get [`TerminalDataSurrogate`](@ref) `n_past_timesteps`."""
get_n_past_timesteps(value::TerminalDataSurrogate) = value.n_past_timesteps
"""Get [`TerminalDataSurrogate`](@ref) `ext`."""
PSY.get_ext(value::TerminalDataSurrogate) = value.ext
"""Get [`TerminalDataSurrogate`](@ref) `internal`."""
PSY.get_internal(value::TerminalDataSurrogate) = value.internal

"""Set [`TerminalDataSurrogate`](@ref) `available`."""
PSY.set_available!(value::TerminalDataSurrogate, val) = value.available = val
"""Set [`TerminalDataSurrogate`](@ref) `bus`."""
PSY.set_bus!(value::TerminalDataSurrogate, val) = value.bus = val
"""Set [`TerminalDataSurrogate`](@ref) `active_power`."""
PSY.set_active_power!(value::TerminalDataSurrogate, val) = value.active_power = val
"""Set [`TerminalDataSurrogate`](@ref) `reactive_power`."""
PSY.set_reactive_power!(value::TerminalDataSurrogate, val) = value.reactive_power = val
"""Set [`TerminalDataSurrogate`](@ref) `active_power_limits`."""
PSY.set_active_power_limits!(value::TerminalDataSurrogate, val) =
    value.active_power_limits = val
"""Set [`TerminalDataSurrogate`](@ref) `reactive_power_limits`."""
PSY.set_reactive_power_limits!(value::TerminalDataSurrogate, val) =
    value.reactive_power_limits = val
"""Set [`TerminalDataSurrogate`](@ref) `internal_voltage`."""
PSY.set_internal_voltage!(value::TerminalDataSurrogate, val) = value.internal_voltage = val
"""Set [`TerminalDataSurrogate`](@ref) `internal_angle`."""
PSY.set_internal_angle!(value::TerminalDataSurrogate, val) = value.internal_angle = val
"""Get [`FullyConnected`](@ref) `θ_ref_frame`."""
set_θ_ref_frame!(value::TerminalDataSurrogate, val) = value.θ_ref_frame = val
"""Set [`TerminalDataSurrogate`](@ref) `ext`."""
PSY.set_ext!(value::TerminalDataSurrogate, val) = value.ext = val
