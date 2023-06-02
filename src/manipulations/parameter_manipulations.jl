const PARAMETER_TYPES = [
    Float64,
    NamedTuple{(:min, :max), Tuple{Float64, Float64}},
    NamedTuple{(:in, :out), Tuple{Float64, Float64}},
]
import Base.:+
import Base.:/
import Base.:*
function Base.:+(
    A::NamedTuple{(:in, :out), Tuple{Float64, Float64}},
    B::NamedTuple{(:in, :out), Tuple{Float64, Float64}},
)
    return (in = (A.in + B.in), out = (A.out + B.out))
end
function Base.:+(
    A::NamedTuple{(:min, :max), Tuple{Float64, Float64}},
    B::NamedTuple{(:min, :max), Tuple{Float64, Float64}},
)
    return (min = (A.min + B.min), max = (A.max + B.max))
end
function Base.:/(
    A::NamedTuple{(:in, :out), Tuple{Float64, Float64}},
    B::R,
) where {R <: Real}
    return (in = (A.in / B), out = (A.out / B))
end
function Base.:/(
    A::NamedTuple{(:min, :max), Tuple{Float64, Float64}},
    B::R,
) where {R <: Real}
    return (min = (A.min / B), max = (A.max / B))
end
function Base.:*(
    A::NamedTuple{(:in, :out), Tuple{Float64, Float64}},
    B::R,
) where {R <: Real}
    return (in = (A.in * B), out = (A.out * B))
end
function Base.:*(
    A::NamedTuple{(:min, :max), Tuple{Float64, Float64}},
    B::R,
) where {R <: Real}
    return (min = (A.min * B), max = (A.max * B))
end

function build_default(T::Type{PowerSystems.GenericBattery})
    component = T(nothing)
    PowerSystems.set_available!(component, true)
    return component
end

function build_default(
    T::Type{PowerSystems.DynamicInverter{C, OC, IC, DC, P, F}},
) where {
    C <: PowerSystems.Converter,
    OC <: PowerSystems.OuterControl,
    IC <: PowerSystems.InnerControl,
    DC <: PowerSystems.DCSource,
    P <: PowerSystems.FrequencyEstimator,
    F <: PowerSystems.Filter,
}
    component = PowerSystems.DynamicInverter(;
        name = "init",
        ω_ref = 1.0,
        converter = build_default(C),
        outer_control = build_default(OC),
        inner_control = build_default(IC),
        dc_source = build_default(DC),
        freq_estimator = build_default(P),
        filter = build_default(F),
    )
    return component
end

build_default(T::Type{C}) where {C <: PowerSystems.Converter} = T(nothing)
build_default(
    T::Type{PowerSystems.OuterControl{A, B}},
) where {A <: PowerSystems.ActivePowerControl, B <: PowerSystems.ReactivePowerControl} =
    PowerSystems.OuterControl(A(nothing), B(nothing))
build_default(T::Type{C}) where {C <: PowerSystems.InnerControl} = T(nothing)
build_default(T::Type{C}) where {C <: PowerSystems.DCSource} = T(nothing)
build_default(T::Type{C}) where {C <: PowerSystems.FrequencyEstimator} = T(nothing)
build_default(T::Type{C}) where {C <: PowerSystems.Filter} = T(nothing)

function add_aggregate_model!(
    sys_add::PowerSystems.System,
    sys_derive::PowerSystems.System,
    devices::Vector{Tuple{DataType, DataType}},
    bus_number::Int64,
)
    for (static_injection_type, dynamic_injection_type) in devices
        df_static = DataFrames.DataFrame()
        df_dynamic = DataFrames.DataFrame()
        for static_injector in PowerSystems.get_components(
            x -> typeof(PowerSystems.get_dynamic_injector(x)) == dynamic_injection_type,
            static_injection_type,
            sys_derive,
        )
            dict_static = get_component_parameter_dictionary(static_injector)
            dynamic_injector = PowerSystems.get_dynamic_injector(static_injector)
            dict_dynamic = get_component_parameter_dictionary(dynamic_injector)
            if isempty(df_static)
                df_static = DataFrames.DataFrame(dict_static)
            else
                push!(df_static, dict_static)
            end
            if isempty(df_dynamic)
                df_dynamic = DataFrames.DataFrame(dict_dynamic)
            else
                push!(df_dynamic, dict_dynamic)
            end
        end

        println(df_static)
        println(df_dynamic)

        static_injector = build_default(static_injection_type)
        dynamic_injector = build_default(dynamic_injection_type)

        dict_static, dict_dynamic = calculate_aggregate_parameters(df_static, df_dynamic)

        println(dict_static)
        println(dict_dynamic)
        set_component_parameter_dictionary!(static_injector, dict_static)
        set_component_parameter_dictionary!(dynamic_injector, dict_dynamic)

        bus_to_add =
            collect(PowerSystems.get_components(x -> PSY.get_number(x) == bus_number, PSY.Bus, sys_add))[1]

        PowerSystems.set_bus!(static_injector, bus_to_add)
        PowerSystems.add_component!(sys_add, static_injector)
        PowerSystems.add_component!(sys_add, dynamic_injector, static_injector)
    end
end

function calculate_aggregate_parameters(
    df_static_injector::DataFrames.DataFrame,
    df_dynamic_injector::DataFrames.DataFrame,
)
    @assert "base_power" in names(df_static_injector)
    @assert "base_power" in names(df_dynamic_injector)
    @assert df_static_injector[!, :base_power] == df_dynamic_injector[!, :base_power]

    dict_static_injector = Dict()
    for col_name in names(df_static_injector)
        if col_name == "base_power"
            dict_static_injector[Symbol(col_name)] = sum(df_static_injector[!, col_name])
        else
            dict_static_injector[Symbol(col_name)] =
                sum(
                    df_static_injector[!, col_name] .* df_static_injector[!, "base_power"],
                ) / sum(df_static_injector[!, "base_power"])
        end
    end

    dict_dynamic_injector = Dict()
    for col_name in names(df_dynamic_injector)
        if col_name == "base_power"
            dict_dynamic_injector[Symbol(col_name)] = sum(df_dynamic_injector[!, col_name])
        else
            dict_dynamic_injector[Symbol(col_name)] =
                sum(
                    df_dynamic_injector[!, col_name] .*
                    df_dynamic_injector[!, "base_power"],
                ) / sum(df_dynamic_injector[!, "base_power"])
        end
    end

    return dict_static_injector, dict_dynamic_injector
end

function set_component_parameter_dictionary!(
    static_injector::SI,
    para_dict::Dict,
) where {SI <: PowerSystems.StaticInjection}
    println(para_dict)
    for (param, value) in para_dict
        println(param)
        println(value)
        _set_parameter!(static_injector, param, value)
    end
    return
end

function get_component_parameter_dictionary(
    static_injector::SI,
) where {SI <: PowerSystems.StaticInjection}
    d = Dict{Symbol, Any}()
    all_parameters = get_component_parameters(static_injector)
    for p in all_parameters
        val = get_parameter(static_injector, p)
        push!(d, p => val)
    end
    return d
end

function set_component_parameter_dictionary!(
    dynamic_inverter::D,
    para_dict::Dict,
) where {D <: PowerSystems.DynamicInverter}
    top_level_params = get_component_parameters(dynamic_inverter)
    for p in top_level_params
        _set_parameter!(dynamic_inverter, p, para_dict[p])
    end

    converter = PowerSystems.get_converter(dynamic_inverter)
    converter_params = get_component_parameters(converter)
    for p in converter_params
        _set_parameter!(converter, p, para_dict[p])
    end

    dc_source = PowerSystems.get_dc_source(dynamic_inverter)
    dc_source_params = get_component_parameters(dc_source)
    for p in dc_source_params
        _set_parameter!(dc_source, p, para_dict[p])
    end

    filter = PowerSystems.get_filter(dynamic_inverter)
    filter_params = get_component_parameters(filter)
    for p in filter_params
        _set_parameter!(filter, p, para_dict[p])
    end

    freq_estimator = PowerSystems.get_freq_estimator(dynamic_inverter)
    freq_estimator_params = get_component_parameters(freq_estimator)
    for p in freq_estimator_params
        _set_parameter!(freq_estimator, p, para_dict[p])
    end

    inner_control = PowerSystems.get_inner_control(dynamic_inverter)
    inner_control_params = get_component_parameters(inner_control)
    for p in inner_control_params
        _set_parameter!(inner_control, p, para_dict[p])
    end

    outer_control_active =
        PowerSystems.get_active_power_control(PowerSystems.get_outer_control(dynamic_inverter))
    outer_control_active_params = get_component_parameters(outer_control_active)
    for p in outer_control_active_params
        _set_parameter!(outer_control_active, p, para_dict[p])
    end

    outer_control_reactive =
        PowerSystems.get_reactive_power_control(PowerSystems.get_outer_control(dynamic_inverter))
    outer_control_reactive_params = get_component_parameters(outer_control_reactive)
    for p in outer_control_reactive_params
        _set_parameter!(outer_control_reactive, p, para_dict[p])
    end
    return
end

function get_component_parameter_dictionary(
    dynamic_inverter::D,
) where {D <: PowerSystems.DynamicInverter}
    d = Dict{Symbol, Any}()
    top_level_params = get_component_parameters(dynamic_inverter)
    for p in top_level_params
        val = get_parameter(dynamic_inverter, p)
        push!(d, p => val)
    end

    converter = PowerSystems.get_converter(dynamic_inverter)
    converter_params = get_component_parameters(converter)
    for p in converter_params
        val = get_parameter(converter, p)
        push!(d, p => val)
    end

    dc_source = PowerSystems.get_dc_source(dynamic_inverter)
    dc_source_params = get_component_parameters(dc_source)
    for p in dc_source_params
        val = get_parameter(dc_source, p)
        push!(d, p => val)
    end

    filter = PowerSystems.get_filter(dynamic_inverter)
    filter_params = get_component_parameters(filter)
    for p in filter_params
        val = get_parameter(filter, p)
        push!(d, p => val)
    end

    freq_estimator = PowerSystems.get_freq_estimator(dynamic_inverter)
    freq_estimator_params = get_component_parameters(freq_estimator)
    for p in freq_estimator_params
        val = get_parameter(freq_estimator, p)
        push!(d, p => val)
    end

    inner_control = PowerSystems.get_inner_control(dynamic_inverter)
    inner_control_params = get_component_parameters(inner_control)
    for p in inner_control_params
        val = get_parameter(inner_control, p)
        push!(d, p => val)
    end

    outer_control_active =
        PowerSystems.get_active_power_control(PowerSystems.get_outer_control(dynamic_inverter))
    outer_control_active_params = get_component_parameters(outer_control_active)
    for p in outer_control_active_params
        val = get_parameter(outer_control_active, p)
        push!(d, p => val)
    end

    outer_control_reactive =
        PowerSystems.get_reactive_power_control(PowerSystems.get_outer_control(dynamic_inverter))
    outer_control_reactive_params = get_component_parameters(outer_control_reactive)
    for p in outer_control_reactive_params
        val = get_parameter(outer_control_reactive, p)
        push!(d, p => val)
    end
    return d
end

function get_component_parameters(component)
    all_fields = collect(fieldnames(typeof(component)))
    return filter!(x -> typeof(getfield(component, x)) in PARAMETER_TYPES, all_fields)
end

function set_parameter_absolute!(
    component::Union{C, IS},
    param::Symbol,
    val::Float64,
) where {C <: PowerSystems.Component, IS <: InfrastructureSystems.DeviceParameter}
    _set_parameter!(component, param, val)
end

function set_parameter_absolute!(
    component::Union{C, IS},
    param::Symbol,
    dist::D,
) where {
    C <: PowerSystems.Component,
    IS <: InfrastructureSystems.DeviceParameter,
    D <: Distributions.Distribution,
}
    sampled_value = rand(dist, 1)[1]
    _set_parameter!(component, param, sampled_value)
end

function set_parameter_relative!(
    component::Union{C, IS},
    param::Symbol,
    mult::Float64,
) where {C <: PowerSystems.Component, IS <: InfrastructureSystems.DeviceParameter}
    new_value = get_parameter(component, param) * mult
    _set_parameter!(component, param, new_value)
end

function set_parameter_relative!(
    component::Union{C, IS},
    param::Symbol,
    dist::D,
) where {
    C <: PowerSystems.Component,
    IS <: InfrastructureSystems.DeviceParameter,
    D <: Distributions.Distribution,
}
    sampled_value = rand(dist, 1)[1]
    set_parameter_relative!(component, param, sampled_value)
end

function _set_parameter!(
    device::Union{C, IS},
    param::Symbol,
    value,
) where {C <: PowerSystems.Component, IS <: InfrastructureSystems.DeviceParameter}
    @assert typeof(value) in PARAMETER_TYPES
    setfield!(device, param, value)
end

function get_parameter(
    device::Union{C, IS},
    param::Symbol,
) where {C <: PowerSystems.Component, IS <: InfrastructureSystems.DeviceParameter}
    return getfield(device, param)
end
