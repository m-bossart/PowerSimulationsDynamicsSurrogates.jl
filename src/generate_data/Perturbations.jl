
abstract type SurrogatePerturbation end

struct PVS <: SurrogatePerturbation
    type::String
    internal_voltage_frequencies::Vector{Float64}
    internal_voltage_coefficients::Vector{Tuple{Float64, Float64}}
    internal_angle_frequencies::Vector{Float64}
    internal_angle_coefficients::Vector{Tuple{Float64, Float64}}
end
function PVS(;
    type = "PVS",
    internal_voltage_frequencies = [0.0],
    internal_voltage_coefficients = [(0.0,0.0)],
    internal_angle_frequencies = [0.0],
    internal_angle_coefficients = [(0.0,0.0)],
)
    PVS(
        type,
        internal_voltage_frequencies,
        internal_voltage_coefficients,
        internal_angle_frequencies,
        internal_angle_coefficients,
    )
end

struct VStep <: SurrogatePerturbation
    type::String
    t_step::Float64
    ΔV::Float64
end

function VStep(; type = "VStep", t_step = 0.0, ΔV = 0.0)
    VStep(type, t_step, ΔV)
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    component::Tuple{String, T},
) where {T <: SurrogatePerturbation}
    @warn "add_surrogate_component! not implemented for this type of surrogate perturbation"
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    component::Tuple{String, PVS},
)
    source_name = component[1]
    perturbation = component[2]
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
        internal_angle_bias = 0.0,  #set in PSID initialization 
        internal_angle_frequencies = perturbation.internal_angle_frequencies,
        internal_angle_coefficients = perturbation.internal_angle_coefficients,
    )
    PSY.add_component!(sys, pvs, source)
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    component::Tuple{String, VStep},
)
    source_name = component[1]
    perturbation = component[2]

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
