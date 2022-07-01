
abstract type SurrogatePerturbation end

struct PVS <: SurrogatePerturbation
    type::String
    source_name::String
    internal_voltage_frequencies::Vector{Float64}
    internal_voltage_coefficients::Vector{Tuple{Float64, Float64}}
    internal_angle_frequencies::Vector{Float64}
    internal_angle_coefficients::Vector{Tuple{Float64, Float64}}
end

function PVS(;
    type = "PVS",
    source_name = "init",
    internal_voltage_frequencies = [0.0],
    internal_voltage_coefficients = [(0.0, 0.0)],
    internal_angle_frequencies = [0.0],
    internal_angle_coefficients = [(0.0, 0.0)],
)
    PVS(
        type,
        source_name,
        internal_voltage_frequencies,
        internal_voltage_coefficients,
        internal_angle_frequencies,
        internal_angle_coefficients,
    )
end

struct VStep <: SurrogatePerturbation
    type::String
    source_name::String
    t_step::Float64
    ΔV::Float64
end

function VStep(; type = "VStep", source_name = "init", t_step = 0.0, ΔV = 0.0)
    VStep(type, source_name, t_step, ΔV)
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::T,
) where {T <: SurrogatePerturbation}
    @warn "add_surrogate_perturbation not implemented for this type of surrogate perturbation"
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::T,
) where {T <: PSID.Perturbation}
    push!(psid_pertubations, perturbation)
end

function add_surrogate_perturbation!(sys::PSY.System, psid_perturbations, perturbation::PVS)
    source_name = perturbation.source_name
    source = PSY.get_component(PSY.Source, sys, source_name)
    if source === nothing
        @error "Source not found - check name!"
    end
    pvs = PSY.PeriodicVariableSource(
        name = PSY.get_name(source),
        R_th = PSY.get_R_th(source),
        X_th = PSY.get_X_th(source),
        internal_voltage_bias = 0.0,    #set in PSID initialization 
        internal_voltage_frequencies = perturbation.internal_voltage_frequencies,
        internal_voltage_coefficients = perturbation.internal_voltage_coefficients,
        internal_angle_bias = 0.0,      #set in PSID initialization 
        internal_angle_frequencies = perturbation.internal_angle_frequencies,
        internal_angle_coefficients = perturbation.internal_angle_coefficients,
    )
    PSY.add_component!(sys, pvs, source)
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::VStep,
)
    source_name = perturbation.source_name
    source = PSY.get_component(PSY.Source, sys, source_name)
    if source === nothing
        @error "Source not found - check name!"
    end
    V0_bus = PSY.get_magnitude(PSY.get_bus(source)) #TODO - update this to calculate the internal voltage, not the bus voltage
    push!(
        psid_perturbations,
        PSID.SourceBusVoltageChange(
            perturbation.t_step,
            source,
            :V_ref,
            V0_bus + perturbation.ΔV,
        ),
    )
end
