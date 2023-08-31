mutable struct TerminalDataSurrogate{
    M <: MLLayer,
    D <: PSY.DynamicInjection,
    S <: DataScaler,
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
    model_architecture::Vector{M}
    model_parameters::Vector{Float64}
    underlying_dynamic_model::D
    data_scaler::S
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
    model_parameters,
    underlying_dynamic_model,
    data_scaler,
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
        model_parameters,
        underlying_dynamic_model,
        data_scaler,
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
    model_parameters,
    underlying_dynamic_model,
    data_scaler,
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
        model_parameters,
        underlying_dynamic_model,
        data_scaler,
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
        model_architecture = [FFNN(nothing)],
        model_parameters = [0.0],
        underlying_dynamic_model = PSY.GenericDER(nothing),
        data_scaler = MinMaxScaler(nothing),
        θ_ref_frame = 0.0,
        n_past_timesteps = 0,
        services = PSY.Service[],
        ext = Dict{String, Any}(),
    )
end

#If function exists in PSY, overload it here. 
PSY.get_name(value::TerminalDataSurrogate) = value.name
PSY.get_available(value::TerminalDataSurrogate) = value.available
PSY.get_bus(value::TerminalDataSurrogate) = value.bus
PSY.get_active_power(value::TerminalDataSurrogate) = value.active_power
PSY.get_reactive_power(value::TerminalDataSurrogate) = value.reactive_power
PSY.get_active_power_limits(value::TerminalDataSurrogate) = value.active_power_limits
PSY.get_reactive_power_limits(value::TerminalDataSurrogate) = value.reactive_power_limits
PSY.get_internal_voltage(value::TerminalDataSurrogate) = value.internal_voltage
PSY.get_internal_angle(value::TerminalDataSurrogate) = value.internal_angle
get_model_architecture(value::TerminalDataSurrogate) = value.model_architecture
get_model_parameters(value::TerminalDataSurrogate) = value.model_parameters
get_underlying_dynamic_model(value::TerminalDataSurrogate) = value.underlying_dynamic_model
get_data_scaler(value::TerminalDataSurrogate) = value.data_scaler
get_θ_ref_frame(value::TerminalDataSurrogate) = value.θ_ref_frame
get_n_past_timesteps(value::TerminalDataSurrogate) = value.n_past_timesteps
PSY.get_ext(value::TerminalDataSurrogate) = value.ext
PSY.get_internal(value::TerminalDataSurrogate) = value.internal

PSY.set_available!(value::TerminalDataSurrogate, val) = value.available = val
PSY.set_bus!(value::TerminalDataSurrogate, val) = value.bus = val
PSY.set_active_power!(value::TerminalDataSurrogate, val) = value.active_power = val
PSY.set_reactive_power!(value::TerminalDataSurrogate, val) = value.reactive_power = val
PSY.set_active_power_limits!(value::TerminalDataSurrogate, val) =
    value.active_power_limits = val
PSY.set_reactive_power_limits!(value::TerminalDataSurrogate, val) =
    value.reactive_power_limits = val
PSY.set_internal_voltage!(value::TerminalDataSurrogate, val) = value.internal_voltage = val
PSY.set_internal_angle!(value::TerminalDataSurrogate, val) = value.internal_angle = val
set_θ_ref_frame!(value::TerminalDataSurrogate, val) = value.θ_ref_frame = val
set_model_parameters!(value::TerminalDataSurrogate, val) = value.model_parameters = val
PSY.set_ext!(value::TerminalDataSurrogate, val) = value.ext = val
