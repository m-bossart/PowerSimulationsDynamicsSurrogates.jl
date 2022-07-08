
abstract type SurrogatePerturbation end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::T,
    sys_aux::PSY.System,
) where {T <: SurrogatePerturbation}
    @warn "add_surrogate_perturbation not implemented for this type of surrogate perturbation"
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::T,
    sys_aux::PSY.System,
) where {T <: PSID.Perturbation}
    push!(psid_pertubations, perturbation)
end

###############################################################################
################################## PVS ########################################
###############################################################################

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

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::PVS,
    sys_aux::PSY.System,
)
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

###############################################################################
################################# VStep #######################################
###############################################################################

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
    perturbation::VStep,
    sys_aux::PSY.System,
)
    source_name = perturbation.source_name
    source = PSY.get_component(PSY.Source, sys, source_name)
    if source === nothing
        @error "Source not found - check name!"
    end
    V0_source = PSY.get_internal_voltage(source)
    push!(
        psid_perturbations,
        PSID.SourceBusVoltageChange(
            perturbation.t_step,
            source,
            :V_ref,
            V0_source + perturbation.ΔV,
        ),
    )
end

###############################################################################
############################ RandomLoadTrip ###################################
###############################################################################

struct RandomLoadTrip <: SurrogatePerturbation
    type::String
    time::Float64
end

function RandomLoadTrip(; type = "RandomLoadTrip", time = 0.0)
    RandomLoadTrip(type, time)
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::RandomLoadTrip,
    sys_aux::PSY.System,
)
    time = perturbation.time
    electric_loads = collect(PSY.get_components(PSY.ElectricLoad, sys_aux))
    if length(electric_loads) === 0
        @error "Trying to trip a load but an electric load not found in system"
        return
    end
    l = rand(electric_loads)
    l_new = PSY.get_component(typeof(l), sys, PSY.get_name(l))
    push!(psid_perturbations, PSID.LoadTrip(time, l_new))
end

###############################################################################
############################ RandomBranchTrip #################################
###############################################################################

struct RandomBranchTrip <: SurrogatePerturbation
    type::String
    time::Float64
end

function RandomBranchTrip(; type = "RandomBranchTrip", time = 0.0)
    RandomBranchTrip(type, time)
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::RandomBranchTrip,
    sys_aux::PSY.System,
)
    time = perturbation.time
    ac_branches = collect(PSY.get_components(PSY.ACBranch, sys_aux))
    if length(ac_branches) === 0
        @error "Trying to trip an AC branch but an ACBranch not found in system"
        return
    end
    b = nothing
    for i in 1:10
        b = rand(ac_branches)
        if _check_if_connected_to_source(b, sys_aux) == false
            break
        else
            b = rand(ac_branches)
        end
    end
    push!(psid_perturbations, PSID.BranchTrip(time, typeof(b), PSY.get_name(b)))
end

function _check_if_connected_to_source(branch, sys)
    from_bus = PSY.get_from(PSY.get_arc(branch))
    to_bus = PSY.get_to(PSY.get_arc(branch))
    for s in PSY.get_components(PSY.Source, sys)
        b = PSY.get_bus(s)
        if b == from_bus || b == to_bus
            return true
        end
    end
    return false
end

###############################################################################
############################ RandomLoadChange #################################
###############################################################################

struct RandomLoadChange <: SurrogatePerturbation
    type::String
    time::Float64
    load_multiplier_range::Tuple{Float64, Float64}
end

function RandomLoadChange(;
    type = "RandomLoadChange",
    time = 0.0,
    load_multiplier_range = (1.0, 1.0),
)
    RandomLoadChange(type, time, load_multiplier_range)
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::RandomLoadChange,
    sys_aux::PSY.System,
)
    time = perturbation.time
    load_multiplier_range = perturbation.load_multiplier_range
    electric_loads = collect(PSY.get_components(PSY.ElectricLoad, sys_aux))
    if length(electric_loads) === 0
        @error "Trying to change a load but an electric load not found in system"
        return
    end
    l = rand(electric_loads)
    l_new = PSY.get_component(typeof(l), sys, PSY.get_name(l))
    multiplier =
        rand() * (load_multiplier_range[2] - load_multiplier_range[1]) +
        load_multiplier_range[1]
    Pnew = PSY.get_active_power(l_new) * multiplier
    Qnew = PSY.get_reactive_power(l_new) * multiplier
    push!(psid_perturbations, PSID.LoadChange(time, l_new, :P_ref, Pnew))
    push!(psid_perturbations, PSID.LoadChange(time, l_new, :Q_ref, Qnew))
end
