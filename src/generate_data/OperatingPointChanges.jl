
abstract type SurrogateOperatingPoint end

function update_operating_point!(
    sys::PSY.System,
    condition::T,
    sys_aux::PSY.System,
) where {T <: SurrogateOperatingPoint}
    @warn "update_operating_point! not implemented for this type of surrogate operation point"
end

###############################################################################
############################## FunctionCall ######################################
###############################################################################

struct FunctionCall <: SurrogateOperatingPoint
    type::String
    f::Any
end

function FunctionCall(; type = "FunctionCall", f = (sys, sys_aux) -> ())
    FunctionCall(type, f)
end

function update_operating_point!(
    sys::PSY.System,
    condition::FunctionCall,
    sys_aux::PSY.System,
)
    condition.f(sys, sys_aux)
    return
end

###############################################################################
############################## DoNothing ######################################
###############################################################################

struct DoNothing <: SurrogateOperatingPoint
    type::String
end

function DoNothing(; type = "DoNothing")
    DoNothing(type)
end

function update_operating_point!(sys::PSY.System, condition::DoNothing, sys_aux::PSY.System)
    return
end

###############################################################################
############################## SetVoltageSource ###############################
###############################################################################

struct SetVoltageSource <: SurrogateOperatingPoint
    type::String
    source_name::String
    internal_voltage::Float64
    internal_angle::Float64
end

function SetVoltageSource(;
    type = "SetVoltageSource",
    source_name = "test-load",
    internal_voltage = 0.0,
    internal_angle = 0.0,
)
    SetVoltageSource(type, source_name, internal_voltage, internal_angle)
end

function update_operating_point!(
    sys::PSY.System,
    condition::SetVoltageSource,
    sys_aux::PSY.System,
)
    source = PSY.get_component(PSY.Source, sys, condition.source_name)
    if source === nothing
        @error "Source not found for SetVoltageSource operating point change"
    end
    PSY.set_internal_voltage!(source, condition.internal_voltage)
    PSY.set_internal_angle!(source, condition.internal_angle)
    return
end

###############################################################################
############################## SetStandardLoad ################################
###############################################################################

struct SetStandardLoad <: SurrogateOperatingPoint
    type::String
    load_name::String
    constant_active_power::Float64
    constant_reactive_power::Float64
    impedance_active_power::Float64
    impedance_reactive_power::Float64
    current_active_power::Float64
    current_reactive_power::Float64
end

function SetStandardLoad(;
    type = "SetStandardLoad",
    load_name = "test-load",
    constant_active_power = 0.0,
    constant_reactive_power = 0.0,
    impedance_active_power = 0.0,
    impedance_reactive_power = 0.0,
    current_active_power = 0.0,
    current_reactive_power = 0.0,
)
    SetStandardLoad(
        type,
        load_name,
        constant_active_power,
        constant_reactive_power,
        impedance_active_power,
        impedance_reactive_power,
        current_active_power,
        current_reactive_power,
    )
end

function update_operating_point!(
    sys::PSY.System,
    condition::SetStandardLoad,
    sys_aux::PSY.System,
)
    load = PSY.get_component(PSY.StandardLoad, sys, condition.load_name)
    if load === nothing
        @error "Load not found for SetStandardLoad operating point change"
    end
    PSY.set_constant_active_power!(load, condition.constant_active_power)
    PSY.set_constant_reactive_power!(load, condition.constant_reactive_power)
    PSY.set_impedance_active_power!(load, condition.impedance_active_power)
    PSY.set_impedance_reactive_power!(load, condition.impedance_reactive_power)
    PSY.set_current_active_power!(load, condition.current_active_power)
    PSY.set_current_reactive_power!(load, condition.current_reactive_power)
    return
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
        if PSY.get_bustype(PSY.get_bus(g_new)) == PSY.ACBusTypes.PQ
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
        if PSY.get_bustype(PSY.get_bus(g_new)) == PSY.ACBusTypes.PQ
            PSY.set_reactive_power!(g_new, PSY.get_reactive_power(g_new) * generation_scale)
        end
    end
    for l in PSY.get_components(PSY.ElectricLoad, sys_aux)  #Search in sys_aux, implement in sys
        l_name = PSY.get_name(l)
        l_new = PSY.get_component(typeof(l), sys, l_name)
        if l_new === nothing
            @error "Device from auxiliary system not found in main system"
        end
        if typeof(l) == PSY.StandardLoad
            PSY.set_impedance_active_power!(
                l_new,
                PSY.get_impedance_active_power(l_new) * load_scale,
            )
            PSY.set_impedance_reactive_power!(
                l_new,
                PSY.get_impedance_reactive_power(l_new) * load_scale,
            )
        elseif typeof(l) == PSY.PowerLoad
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
        P =
            rand() * (generator_power_range[2] - generator_power_range[1]) +
            generator_power_range[1]
        PSY.set_active_power!(g_new, P)
        Vm =
            rand() * (generator_voltage_range[2] - generator_voltage_range[1]) +
            generator_voltage_range[1]
        PSY.set_magnitude!(PSY.get_bus(g_new), Vm)
    end
    for l in PSY.get_components(PSY.ElectricLoad, sys_aux)  #Search in sys_aux, implement in sys
        l_name = PSY.get_name(l)
        l_new = PSY.get_component(typeof(l), sys, l_name)
        if l_new === nothing
            @error "Device from auxiliary system not found in main system"
        end
        if typeof(l) <: PSY.StaticLoad
            Pload =
                PSY.get_impedance_active_power(l_new) * (
                    rand() * (load_multiplier_range[2] - load_multiplier_range[1]) +
                    load_multiplier_range[1]
                )
            PSY.set_impedance_active_power!(l_new, Pload)
            Qload =
                PSY.get_impedance_reactive_power(l_new) * (
                    rand() * (load_multiplier_range[2] - load_multiplier_range[1]) +
                    load_multiplier_range[1]
                )
            PSY.set_impedance_reactive_power!(l_new, Qload)
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
