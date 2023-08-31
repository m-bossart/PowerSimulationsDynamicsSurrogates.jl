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
    push!(psid_perturbations, perturbation)
end

# Cannot trip a dynamic line in PSID - can remove this special case once this feature is added. 
# This workaround will not work for simulation with multiple perturbations
function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::PSID.BranchTrip,
    sys_aux::PSY.System,
)
    if PSY.get_component(PSY.DynamicBranch, sys, perturbation.branch_name) !== nothing
        dyn_branch = PSY.get_component(PSY.DynamicBranch, sys, perturbation.branch_name)
        branch = dyn_branch.branch
        @error "removing this dynamic line in order to trip it", perturbation.branch_name
        PSY.remove_component!(sys, dyn_branch)
        PSY.add_component!(sys, branch)
    end
    push!(psid_perturbations, perturbation)
end

#NOTE: Remove and use built in option after this issue is addressed: https://github.com/NREL-Sienna/PowerSimulationsDynamics.jl/issues/336
###############################################################################
################################## PerturbStateAbsolute #######################
###############################################################################
"""
    function PerturbState(
        time::Float64,
        index::Int,
        value::Float64,
    )
Allows the user to modify the state `index` by adding `value`. The user should modify dynamic states only, since algebraic state may require to do a reinitialization.
# Arguments:
- `time::Float64` : Defines when the modification of the state will happen. This time should be inside the time span considered in the Simulation.
- `index::Int` : Defines which state index you want to modify
- `value::Float64` : Defines how much the state will increase in value
"""
mutable struct PerturbStateAbsolute <: PSID.Perturbation
    time::Float64
    index::Int
    value::Float64
end

function PSID.get_affect(::PSID.SimulationInputs, ::PSY.System, pert::PerturbStateAbsolute)
    return (integrator) -> begin
        @debug "Modifying state"
        integrator.u[pert.index] = pert.value
        return
    end
end

###############################################################################
################################## ParameterChange ############################
###############################################################################

mutable struct ParameterChange <: PSID.Perturbation
    time::Float64
    device::PSY.DynamicInjection
    parameter::Symbol
    value::Float64
end

function PSID.get_affect(inputs::PSID.SimulationInputs, ::PSY.System, pert::ParameterChange)
    wrapped_device_ix = PSID._find_device_index(inputs, pert.device)
    return (integrator) -> begin
        wrapped_device = PSID.get_dynamic_injectors(integrator.p)[wrapped_device_ix]
        @debug "Changing $(PSY.get_name(wrapped_device)) $(pert.parameter) to $(pert.value)"
        machine = PSY.get_machine(PSID.get_device(wrapped_device))
        PSY.set_eq_p!(machine, pert.value)
        return
    end
end

###############################################################################
################################## BranchTrip #################################
###############################################################################

struct LineTrip <: SurrogatePerturbation
    type::String
    time::Float64
    line_name::String
end

function LineTrip(; type = "LineTrip", time = 0.0, line_name = "")
    LineTrip(type, time, line_name)
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::LineTrip,
    sys_aux::PSY.System,
)
    psid_perturbation = PSID.BranchTrip(perturbation.time, PSY.Line, perturbation.line_name)
    if PSY.get_component(PSY.DynamicBranch, sys, perturbation.line_name) !== nothing
        dyn_branch = PSY.get_component(PSY.DynamicBranch, sys, perturbation.line_name)
        branch = dyn_branch.branch
        @error "removing this dynamic line in order to trip it", perturbation.line_name
        PSY.remove_component!(sys, dyn_branch)
        PSY.add_component!(sys, branch)
    end
    push!(psid_perturbations, psid_perturbation)
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
################################## Chirp ######################################
###############################################################################
struct Chirp <: SurrogatePerturbation
    type::String
    source_name::String
    ω1::Float64
    ω2::Float64
    tstart::Float64
    N::Float64
    V_amp::Float64
    ω_amp::Float64
end

function Chirp(;
    type = "Chirp",
    source_name = "init",
    ω1 = 0.0,
    ω2 = 0.0,
    tstart = 0.0,
    N = 0.0,
    V_amp = 0.0,
    ω_amp = 0.0,
)
    Chirp(type, source_name, ω1, ω2, tstart, N, V_amp, ω_amp)
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::Chirp,
    sys_aux::PSY.System,
)
    source_name = perturbation.source_name
    source = PSY.get_component(PSY.Source, sys, source_name)
    if source === nothing
        @error "Source not found - check name!"
    end
    chirp = FrequencyChirpVariableSource(
        name = PSY.get_name(source),
        R_th = PSY.get_R_th(source),
        X_th = PSY.get_X_th(source),
        ω1 = perturbation.ω1,
        ω2 = perturbation.ω2,
        tstart = perturbation.tstart,
        N = perturbation.N,
        V_amp = perturbation.V_amp,
        ω_amp = perturbation.ω_amp,
    )
    PSY.add_component!(sys, chirp, source)
end

###############################################################################
################################# RandomSourceVoltageChange ###################
###############################################################################

struct RandomSourceVoltageChange <: SurrogatePerturbation
    type::String
    source_name::String
    t_step::Float64
    vm_bounds::Vector{Float64}
    vθ_bounds::Vector{Float64}
    vm_index::Int   #Remove once https://github.com/NREL-Sienna/PowerSimulationsDynamics.jl/issues/336 is implemented. 
    vθ_index::Int   #Remove once https://github.com/NREL-Sienna/PowerSimulationsDynamics.jl/issues/336 is implemented. 
end

function RandomSourceVoltageChange(;
    type = "VStep",
    source_name = "init",
    t_step = 0.0,
    vm_bounds = [0.0, 1.0],
    vθ_bounds = [0.0, 1.0],
    vm_index = 1,
    vθ_index = 2,
)
    RandomSourceVoltageChange(
        type,
        source_name,
        t_step,
        vm_bounds,
        vθ_bounds,
        vm_index,
        vθ_index,
    )
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::RandomSourceVoltageChange,
    sys_aux::PSY.System,
)
    source_name = perturbation.source_name
    source = PSY.get_component(PSY.Source, sys, source_name)
    if source === nothing
        @error "Source not found - check name!"
    end
    vm = rand(Distributions.Uniform(perturbation.vm_bounds[1], perturbation.vm_bounds[2]))
    vθ = rand(Distributions.Uniform(perturbation.vθ_bounds[1], perturbation.vθ_bounds[2]))
    dyn_source = PSY.get_component(FrequencySource, sys, source_name)
    if dyn_source === nothing
        push!(
            psid_perturbations,
            PSID.SourceBusVoltageChange(perturbation.t_step, source, :V_ref, vm),
        )
        push!(
            psid_perturbations,
            PSID.SourceBusVoltageChange(perturbation.t_step, source, :θ_ref, vθ),
        )
    else
        @error "Source has a dynamic model"
    end
end

###############################################################################
################################# RandomSourceFrequencyChange #################
###############################################################################

struct RandomSourceFrequencyChange <: SurrogatePerturbation
    type::String
    source_name::String
    t_step::Float64
    ω_bounds::Vector{Float64}
    ω_index::Int   #Remove once https://github.com/NREL-Sienna/PowerSimulationsDynamics.jl/issues/336 is implemented. 
end

function RandomSourceFrequencyChange(;
    type = "VStep",
    source_name = "init",
    t_step = 0.0,
    ω_bounds = [0.0, 1.0],
    ω_index = 1,
)
    RandomSourceFrequencyChange(type, source_name, t_step, ω_bounds, ω_index)
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::RandomSourceFrequencyChange,
    sys_aux::PSY.System,
)
    source_name = perturbation.source_name
    source = PSY.get_component(PSY.Source, sys, source_name)
    if source === nothing
        @error "Source not found - check name!"
    end
    ω = rand(Distributions.Uniform(perturbation.ω_bounds[1], perturbation.ω_bounds[2]))
    dyn_source = PSY.get_component(FrequencySource, sys, source_name)
    if dyn_source === nothing
        #TODO - implement a change in the Source frequency reference (need to add frequency ref to source and change PSID handling of system frequency. )
    else
        @error "FrequencySource component not found."
    end
end

###############################################################################
################################# VStep #######################################
###############################################################################

struct VStep <: SurrogatePerturbation
    type::String
    source_name::String
    t_step::Float64
    V_step::Float64
end

function VStep(; type = "VStep", source_name = "init", t_step = 0.0, V_step = 0.0)
    VStep(type, source_name, t_step, V_step)
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
    #V0_source = PSY.get_internal_voltage(source) 
    push!(
        psid_perturbations,
        PSID.SourceBusVoltageChange(
            perturbation.t_step,
            source,
            :V_ref,
            perturbation.V_step,
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
    #Cannot trip a dynamic branch in PSID, replace with Line
    if typeof(b) == PSY.DynamicBranch
        line = PSY.get_branch(deepcopy(b))
        PSY.remove_component!(sys, b)
        PSY.add_component!(sys, line)
        push!(psid_perturbations, PSID.BranchTrip(time, typeof(line), PSY.get_name(line)))
    else
        push!(psid_perturbations, PSID.BranchTrip(time, typeof(b), PSY.get_name(b)))
    end
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
    electric_loads = collect(PSY.get_components(PSY.StaticLoad, sys_aux))
    if length(electric_loads) === 0
        @error "Trying to change a load but an electric load not found in system"
        return
    end
    l = rand(electric_loads)
    l_new = PSY.get_component(typeof(l), sys, PSY.get_name(l))
    multiplier =
        rand() * (load_multiplier_range[2] - load_multiplier_range[1]) +
        load_multiplier_range[1]
    Pnew = PSY.get_impedance_active_power(l_new) * multiplier
    Qnew = PSY.get_impedance_reactive_power(l_new) * multiplier
    push!(psid_perturbations, PSID.LoadChange(time, l_new, :P_ref_impedance, Pnew))
    push!(psid_perturbations, PSID.LoadChange(time, l_new, :Q_ref_impedance, Qnew))
end

###############################################################################
############################ RandomGenerationChange ###########################
###############################################################################

struct RandomGenerationChange <: SurrogatePerturbation
    type::String
    time::Float64
    P_multiplier_range::Tuple{Float64, Float64}
end

function RandomGenerationChange(;
    type = "RandomLoadChange",
    time = 0.0,
    P_multiplier_range = (1.0, 1.0),
)
    RandomGenerationChange(type, time, P_multiplier_range)
end

function add_surrogate_perturbation!(
    sys::PSY.System,
    psid_perturbations,
    perturbation::RandomGenerationChange,
    sys_aux::PSY.System,
)
    time = perturbation.time
    P_multiplier_range = perturbation.P_multiplier_range
    static_injectors = collect(
        PSY.get_components(
            x ->
                (PSY.get_dynamic_injector(x) !== nothing) &&
                    PSY.get_available(x) == true,
            PSY.StaticInjection,
            sys_aux,
        ),
    )
    println(PSY.get_name.(static_injectors))
    if length(static_injectors) === 0
        @error "Trying to change a dynamic injector but a dynamic injector not found in system"
        return
    end
    s = rand(static_injectors)
    s_new = PSY.get_component(typeof(s), sys, PSY.get_name(s))
    multiplier =
        rand() * (P_multiplier_range[2] - P_multiplier_range[1]) + P_multiplier_range[1]
    Pnew = PSY.get_active_power(s_new) * multiplier
    #Qnew = PSY.get_impedance_reactive_power(l_new) * multiplier
    println(
        PSID.ControlReferenceChange(time, PSY.get_dynamic_injector(s_new), :P_ref, Pnew),
    )
    push!(
        psid_perturbations,
        PSID.ControlReferenceChange(time, PSY.get_dynamic_injector(s_new), :P_ref, Pnew),
    )
    #push!(psid_perturbations, PSID.LoadChange(time, d_new, :Q_ref_impedance, Qnew))
end
