
abstract type SurrogateOperatingPoint end

function update_operating_point!(
    sys::PSY.System,
    condition::T,
    sys_aux::PSY.System,
) where {T <: SurrogateOperatingPoint}
    @warn "update_operating_point! not implemented for this type of surrogate operation point"
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
        PSY.set_active_power!(g_new, PSY.get_active_power(g_new) * generation_scale)
    end
    for l in PSY.get_components(PSY.ElectricLoad, sys_aux)  #Search in sys_aux, implement in sys
        l_name = PSY.get_name(l)
        l_new = PSY.get_component(typeof(l), sys, l_name)
        PSY.set_active_power!(l_new, PSY.get_active_power(l_new) * load_scale)
        PSY.set_reactive_power!(l_new, PSY.get_reactive_power(l_new) * load_scale)
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
    end
end
