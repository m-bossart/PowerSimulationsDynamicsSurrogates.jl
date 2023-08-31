mutable struct PhysicsInformedSurrogate{
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
    n_past_timesteps::Int64
    services::Vector{PSY.Service}
    ext::Dict{String, Any}
    internal::IS.InfrastructureSystemsInternal
end

function PhysicsInformedSurrogate(
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
    n_past_timesteps,
    services = PSY.Service[],
    ext = Dict{String, Any}(),
)
    PhysicsInformedSurrogate(
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
        n_past_timesteps,
        services,
        ext,
        IS.InfrastructureSystemsInternal(),
    )
end

function PhysicsInformedSurrogate(;
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
    n_past_timesteps,
    services = PSY.Service[],
    ext = Dict{String, Any}(),
)
    PhysicsInformedSurrogate(
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
        n_past_timesteps,
        services,
        ext,
        IS.InfrastructureSystemsInternal(),
    )
end

# Constructor for demo purposes; non-functional.
function PhysicsInformedSurrogate(::Nothing)
    PhysicsInformedSurrogate(;
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
        n_past_timesteps = 0,
        services = PSY.Service[],
        ext = Dict{String, Any}(),
    )
end

#If function exists in PSY, overload it here. 
PSY.get_name(value::PhysicsInformedSurrogate) = value.name
PSY.get_available(value::PhysicsInformedSurrogate) = value.available
PSY.get_bus(value::PhysicsInformedSurrogate) = value.bus
PSY.get_active_power(value::PhysicsInformedSurrogate) = value.active_power
PSY.get_reactive_power(value::PhysicsInformedSurrogate) = value.reactive_power
PSY.get_active_power_limits(value::PhysicsInformedSurrogate) = value.active_power_limits
PSY.get_reactive_power_limits(value::PhysicsInformedSurrogate) = value.reactive_power_limits
PSY.get_internal_voltage(value::PhysicsInformedSurrogate) = value.internal_voltage
PSY.get_internal_angle(value::PhysicsInformedSurrogate) = value.internal_angle
get_model_architecture(value::PhysicsInformedSurrogate) = value.model_architecture
get_model_parameters(value::PhysicsInformedSurrogate) = value.model_parameters
get_underlying_dynamic_model(value::PhysicsInformedSurrogate) =
    value.underlying_dynamic_model
get_data_scaler(value::PhysicsInformedSurrogate) = value.data_scaler
get_n_past_timesteps(value::PhysicsInformedSurrogate) = value.n_past_timesteps
PSY.get_ext(value::PhysicsInformedSurrogate) = value.ext
PSY.get_internal(value::PhysicsInformedSurrogate) = value.internal

PSY.set_available!(value::PhysicsInformedSurrogate, val) = value.available = val
PSY.set_bus!(value::PhysicsInformedSurrogate, val) = value.bus = val
PSY.set_active_power!(value::PhysicsInformedSurrogate, val) = value.active_power = val
PSY.set_reactive_power!(value::PhysicsInformedSurrogate, val) = value.reactive_power = val
PSY.set_active_power_limits!(value::PhysicsInformedSurrogate, val) =
    value.active_power_limits = val
PSY.set_reactive_power_limits!(value::PhysicsInformedSurrogate, val) =
    value.reactive_power_limits = val
PSY.set_internal_voltage!(value::PhysicsInformedSurrogate, val) =
    value.internal_voltage = val
PSY.set_internal_angle!(value::PhysicsInformedSurrogate, val) = value.internal_angle = val
set_model_parameters!(value::PhysicsInformedSurrogate, val) = value.model_parameters = val
PSY.set_ext!(value::PhysicsInformedSurrogate, val) = value.ext = val
