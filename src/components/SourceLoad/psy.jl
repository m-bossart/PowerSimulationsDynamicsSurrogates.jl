mutable struct SourceLoad <: PSY.StaticInjection
    name::String
    available::Bool
    bus::PSY.Bus
    base_power::Float64
    constant_active_power::Float64
    constant_reactive_power::Float64
    impedance_active_power::Float64
    impedance_reactive_power::Float64
    current_active_power::Float64
    current_reactive_power::Float64
    max_constant_active_power::Float64
    max_constant_reactive_power::Float64
    max_impedance_active_power::Float64
    max_impedance_reactive_power::Float64
    max_current_active_power::Float64
    max_current_reactive_power::Float64
    dynamic_injector::Union{Nothing, PSY.DynamicInjection}
    services::Vector{PSY.Service}
    ext::Dict{String, Any}
    internal::PSY.InfrastructureSystemsInternal
end

function SourceLoad(
    name,
    available,
    bus,
    base_power,
    constant_active_power,
    constant_reactive_power,
    impedance_active_power,
    impedance_reactive_power,
    current_active_power,
    current_reactive_power,
    max_constant_active_power,
    max_constant_reactive_power,
    max_impedance_active_power,
    max_impedance_reactive_power,
    max_current_active_power,
    max_current_reactive_power,
    dynamic_injector = nothing,
    services = PSY.Device[],
    ext = Dict{String, Any}(),
)
    SourceLoad(
        name,
        available,
        bus,
        base_power,
        constant_active_power,
        constant_reactive_power,
        impedance_active_power,
        impedance_reactive_power,
        current_active_power,
        current_reactive_power,
        max_constant_active_power,
        max_constant_reactive_power,
        max_impedance_active_power,
        max_impedance_reactive_power,
        max_current_active_power,
        max_current_reactive_power,
        dynamic_injector,
        services,
        ext,
        PSY.InfrastructureSystemsInternal(),
    )
end

function SourceLoad(;
    name,
    available,
    bus,
    base_power,
    constant_active_power,
    constant_reactive_power,
    impedance_active_power,
    impedance_reactive_power,
    current_active_power,
    current_reactive_power,
    max_constant_active_power,
    max_constant_reactive_power,
    max_impedance_active_power,
    max_impedance_reactive_power,
    max_current_active_power,
    max_current_reactive_power,
    dynamic_injector = nothing,
    services = PSY.Device[],
    ext = Dict{String, Any}(),
    internal = PSY.InfrastructureSystemsInternal(),
)
    SourceLoad(
        name,
        available,
        bus,
        base_power,
        constant_active_power,
        constant_reactive_power,
        impedance_active_power,
        impedance_reactive_power,
        current_active_power,
        current_reactive_power,
        max_constant_active_power,
        max_constant_reactive_power,
        max_impedance_active_power,
        max_impedance_reactive_power,
        max_current_active_power,
        max_current_reactive_power,
        dynamic_injector,
        services,
        ext,
        internal,
    )
end

# Constructor for demo purposes; non-functional.
function SourceLoad(::Nothing)
    SourceLoad(;
        name = "init",
        available = false,
        bus = PSY.Bus(nothing),
        base_poewr = 0.0,
        constant_active_power = 0.0,
        constant_reactive_power = 0.0,
        impedance_active_power = 0.0,
        impedance_reactive_power = 0.0,
        current_active_power = 0.0,
        current_reactive_power = 0.0,
        max_constant_active_power = 0.0,
        max_constant_reactive_power = 0.0,
        max_impedance_active_power = 0.0,
        max_impedance_reactive_power = 0.0,
        max_current_active_power = 0.0,
        max_current_reactive_power = 0.0,
        dynamic_injector = nothing,
        services = PSY.Device[],
        ext = Dict{String, Any}(),
    )
end

PSY.get_name(value::SourceLoad) = value.name
PSY.get_available(value::SourceLoad) = value.available
PSY.get_bus(value::SourceLoad) = value.bus
PSY.get_base_power(value::SourceLoad) = value.base_power
PSY.get_constant_active_power(value::SourceLoad) = value.constant_active_power
PSY.get_constant_reactive_power(value::SourceLoad) = value.constant_reactive_power
PSY.get_impedance_active_power(value::SourceLoad) = value.impedance_active_power
PSY.get_impedance_reactive_power(value::SourceLoad) = value.impedance_reactive_power
PSY.get_current_active_power(value::SourceLoad) = value.current_active_power
PSY.get_current_reactive_power(value::SourceLoad) = value.current_reactive_power
PSY.get_max_constant_active_power(value::SourceLoad) = value.max_constant_active_power
PSY.get_max_constant_reactive_power(value::SourceLoad) = value.max_constant_reactive_power
PSY.get_max_impedance_active_power(value::SourceLoad) = value.max_impedance_active_power
PSY.get_max_impedance_reactive_power(value::SourceLoad) = value.max_impedance_reactive_power
PSY.get_max_current_active_power(value::SourceLoad) = value.max_current_active_power
PSY.get_max_current_reactive_power(value::SourceLoad) = value.max_current_reactive_power
PSY.get_dynamic_injector(value::SourceLoad) = value.dynamic_injector
PSY.get_services(value::SourceLoad) = value.services
PSY.get_ext(value::SourceLoad) = value.ext
PSY.get_internal(value::SourceLoad) = value.internal

#ASSUME ONLY CONSTANT IMPEDANCE 
PSY.get_active_power(value::SourceLoad) = value.impedance_active_power
PSY.get_reactive_power(value::SourceLoad) = value.impedance_reactive_power

PSY.set_available!(value::SourceLoad, val) = value.available = val
PSY.set_bus!(value::SourceLoad, val) = value.bus = val
PSY.set_constant_active_power!(value::SourceLoad, val) = value.constant_active_power = val
PSY.set_constant_reactive_power!(value::SourceLoad, val) =
    value.constant_reactive_power = val
PSY.set_impedance_active_power!(value::SourceLoad, val) = value.impedance_active_power = val
PSY.set_impedance_reactive_power!(value::SourceLoad, val) =
    value.impedance_reactive_power = val
PSY.set_current_active_power!(value::SourceLoad, val) = value.current_active_power = val
PSY.set_current_reactive_power!(value::SourceLoad, val) = value.current_reactive_power = val
PSY.set_max_constant_active_power!(value::SourceLoad, val) =
    value.max_constant_active_power = val
PSY.set_max_constant_reactive_power!(value::SourceLoad, val) =
    value.max_constant_reactive_power = val
PSY.set_max_impedance_active_power!(value::SourceLoad, val) =
    value.max_impedance_active_power = val
PSY.set_max_impedance_reactive_power!(value::SourceLoad, val) =
    value.max_impedance_reactive_power = val
PSY.set_max_current_active_power!(value::SourceLoad, val) =
    value.max_current_active_power = val
PSY.set_max_current_reactive_power!(value::SourceLoad, val) =
    value.max_current_reactive_power = val
PSY.set_services!(value::SourceLoad, val) = value.services = val
PSY.set_ext!(value::SourceLoad, val) = value.ext = val

#ASSUME ONLY CONSTANT IMPEDANCE 
PSY.set_active_power!(value::SourceLoad, val) = value.impedance_active_power = val
PSY.set_reactive_power!(value::SourceLoad, val) = value.impedance_reactive_power = val
