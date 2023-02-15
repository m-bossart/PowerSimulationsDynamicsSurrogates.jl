
abstract type SurrogateOperatingPoint end

function update_operating_point!(
    sys::PSY.System,
    condition::T,
    sys_aux::PSY.System,
) where {T <: SurrogateOperatingPoint}
    @warn "update_operating_point! not implemented for this type of surrogate operation point"
end

###############################################################################
############################## ScaleSource ####################################
###############################################################################

struct ScaleSource <: SurrogateOperatingPoint
    type::String
    source_name::String
    V_scale::Float64
    θ_scale::Float64
    P_scale::Float64
    Q_scale::Float64
end

function ScaleSource(;
    type = "ScaleSource",
    source_name = "init",
    V_scale = 1.0,
    θ_scale = 1.0,
    P_scale = 1.0,
    Q_scale = 1.0,
)
    ScaleSource(type, source_name, V_scale, θ_scale, P_scale, Q_scale)
end

function update_operating_point!(
    sys::PSY.System,
    condition::ScaleSource,
    sys_aux::PSY.System,
)
    source = PSY.get_component(PSY.Source, sys, condition.source_name)
    if source === nothing
        @error "Source not found for ScaleSource operating point change"
    end
    PSY.set_active_power!(source, PSY.get_active_power(source) * condition.P_scale)
    PSY.set_reactive_power!(source, PSY.get_reactive_power(source) * condition.Q_scale)

    bus = PSY.get_bus(source)
    PSY.set_magnitude!(bus, PSY.get_magnitude(bus) * condition.V_scale)
    PSY.set_angle!(bus, PSY.get_angle(bus) * condition.θ_scale)
    return
end

###############################################################################
######################### GenerationLoadScale #################################
###############################################################################

struct GenerationLoadScale <: SurrogateOperatingPoint
    type::String
    generation_scale::Float64
    load_scale::Float64
end

function GenerationLoadScale(;
    type = "GenerationLoadScale",
    generation_scale = 0.0,
    load_scale = 0.0,
)
    GenerationLoadScale(type, generation_scale, load_scale)
end

function update_operating_point!(
    sys::PSY.System,
    condition::GenerationLoadScale,
    sys_aux::PSY.System,
)
    generation_scale = condition.generation_scale
    load_scale = condition.generation_scale
    for g in PSY.get_components(PSY.Generator, sys_aux) #Search in sys_aux, implement in sys
        g_name = PSY.get_name(g)
        g_new = PSY.get_component(typeof(g), sys, g_name)
        if g_new === nothing
            @error "Device from auxiliary system not found in main system"
        end
        PSY.set_active_power!(g_new, PSY.get_active_power(g_new) * generation_scale)
        if PSY.get_bustype(PSY.get_bus(g_new)) == PSY.BusTypes.PQ
            PSY.set_reactive_power!(g_new, PSY.get_reactive_power(g_new) * generation_scale)
        end
    end
    for g in PSY.get_components(PSY.Storage, sys_aux) #Search in sys_aux, implement in sys
        g_name = PSY.get_name(g)
        g_new = PSY.get_component(typeof(g), sys, g_name)
        if g_new === nothing
            @error "Device from auxiliary system not found in main system"
        end
        PSY.set_active_power!(g_new, PSY.get_active_power(g_new) * generation_scale)
        if PSY.get_bustype(PSY.get_bus(g_new)) == PSY.BusTypes.PQ
            PSY.set_reactive_power!(g_new, PSY.get_reactive_power(g_new) * generation_scale)
        end
    end
    for l in PSY.get_components(PSY.ElectricLoad, sys_aux)  #Search in sys_aux, implement in sys
        l_name = PSY.get_name(l)
        l_new = PSY.get_component(typeof(l), sys, l_name)
        if l_new === nothing
            @error "Device from auxiliary system not found in main system"
        end
        if typeof(l) <: PSY.StaticLoad
            PSY.set_active_power!(l_new, PSY.get_active_power(l_new) * load_scale)
            PSY.set_reactive_power!(l_new, PSY.get_reactive_power(l_new) * load_scale)
        elseif typeof(l) <: PSY.FixedAdmittance
            PSY.set_Y!(l_new, PSY.get_Y(l_new) * load_scale)
        else
            @error "Not sure how to scale load model"
        end
    end
end

###############################################################################
###################### RandomOperatingPointXiao ###############################
###############################################################################

struct RandomOperatingPointXiao <: SurrogateOperatingPoint
    type::String
    generator_voltage_range::Tuple{Float64, Float64}
    generator_power_range::Tuple{Float64, Float64}
    load_multiplier_range::Tuple{Float64, Float64}
end

function RandomOperatingPointXiao(;
    type = "RandomOperatingPointXiao",
    generator_voltage_range = (0.94, 1.06),
    generator_power_range = (0.0, 1.0),
    load_multiplier_range = (0.5, 1.5),
)
    RandomOperatingPointXiao(
        type,
        generator_voltage_range,
        generator_power_range,
        load_multiplier_range,
    )
end

function update_operating_point!(
    sys::PSY.System,
    condition::RandomOperatingPointXiao,
    sys_aux::PSY.System,
)
    generator_voltage_range = condition.generator_voltage_range
    generator_power_range = condition.generator_power_range
    load_multiplier_range = condition.load_multiplier_range

    for g in PSY.get_components(PSY.Generator, sys_aux) #Search in sys_aux, implement in sys
        g_name = PSY.get_name(g)
        g_new = PSY.get_component(typeof(g), sys, g_name)
        if g_new === nothing
            @error "Device from auxiliary system not found in main system"
        end
        PSY.set_active_power!(
            g_new,
            rand() * (generator_power_range[2] - generator_power_range[1]) +
            generator_power_range[1],
        )
        PSY.set_magnitude!(
            PSY.get_bus(g_new),
            rand() * (generator_voltage_range[2] - generator_voltage_range[1]) +
            generator_voltage_range[1],
        )
    end
    for l in PSY.get_components(PSY.ElectricLoad, sys_aux)  #Search in sys_aux, implement in sys
        l_name = PSY.get_name(l)
        l_new = PSY.get_component(typeof(l), sys, l_name)
        if l_new === nothing
            @error "Device from auxiliary system not found in main system"
        end
        if typeof(l) <: PSY.StaticLoad
            PSY.set_active_power!(
                l_new,
                PSY.get_active_power(l_new) * (
                    rand() * (load_multiplier_range[2] - load_multiplier_range[1]) +
                    load_multiplier_range[1]
                ),
            )
            PSY.set_reactive_power!(
                l_new,
                PSY.get_reactive_power(l_new) * (
                    rand() * (load_multiplier_range[2] - load_multiplier_range[1]) +
                    load_multiplier_range[1]
                ),
            )
        elseif typeof(l) <: PSY.FixedAdmittance
            PSY.set_Y!(
                l_new,
                PSY.get_Y(l_new) * (
                    rand() * (load_multiplier_range[2] - load_multiplier_range[1]) +
                    load_multiplier_range[1]
                ),
            )
        else
            @error "Not sure how to scale load model"
        end
    end
end
