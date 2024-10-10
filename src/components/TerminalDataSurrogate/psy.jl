"""
    mutable struct TerminalDataSurrogate <: DirectSolutionSurrogate <: PSID.DynamicInjection
        name::String
        τ::Float64  
        window_size::Int64  
        base_power::Float64
        states::Vector{Symbol}
        n_states::Int
        ext::Dict{String, Any}
        internal::InfrastructureSystemsInternal
    end

Experimental surrogate

# Arguments
- `name::String`
- `τ::Float64`: The fixed time interval between delayed inputs
- `window_size::Int64`: The number of delaye inputs
- `base_power::Float64`: Base power
- `states::Vector{Symbol}`: 
- `n_states::Int`:
- `n_states::Vector{StateTypes}`: 
- `ext::Dict{String, Any}`
- `internal::InfrastructureSystemsInternal`: power system internal reference, do not modify
"""
mutable struct TerminalDataSurrogate <: DirectSolutionSurrogate
    name::String
    τ::Float64
    window_size::Int64
    steadystate_offset_correction::Bool
    trained_voltage_range::Tuple{Float64, Float64}
    fc::Float64
    base_power::Float64
    states::Vector{Symbol}
    n_states::Int
    states_types::Vector{PSY.StateTypes}
    ext::Dict{String, Any}
    internal::IS.InfrastructureSystemsInternal
end

function TerminalDataSurrogate(
    name,
    τ = 0.0,
    window_size = 1,
    steadystate_offset_correction = true,
    trained_voltage_range = (0.0, 1.0),
    fc = 0.0,
    base_power = 100.0,
    ext = Dict{String, Any}(),
)
    TerminalDataSurrogate(
        name,
        τ,
        window_size,
        steadystate_offset_correction,
        trained_voltage_range,
        fc,
        base_power,
        [:ir, :ii, :vr, :vi],
        4,
        [
            PSY.StateTypes.Differential,
            PSY.StateTypes.Differential,
            PSY.StateTypes.Differential,
            PSY.StateTypes.Differential,
        ],
        ext,
        IS.InfrastructureSystemsInternal(),
    )
end

function TerminalDataSurrogate(;
    name,
    τ = 0.0,
    window_size = 1,
    steadystate_offset_correction = true,
    trained_voltage_range = (0.0, 1.0),
    fc = 0.0,
    base_power = 100.0,
    states = [:ir, :ii, :vr, :vi],
    n_states = 4,
    states_types = [
        PSY.StateTypes.Differential,
        PSY.StateTypes.Differential,
        PSY.StateTypes.Differential,
        PSY.StateTypes.Differential,
    ],
    ext = Dict{String, Any}(),
    internal = IS.InfrastructureSystemsInternal(),
)
    TerminalDataSurrogate(
        name,
        τ,
        window_size,
        steadystate_offset_correction,
        trained_voltage_range,
        fc,
        base_power,
        states,
        n_states,
        states_types,
        ext,
        internal,
    )
end

# Constructor for demo purposes; non-functional.
function TerminalDataSurrogate(::Nothing)
    TerminalDataSurrogate(;
        name = "init",
        τ = 0.0,
        window_size = 1,
        steadystate_offset_correction = true,
        trained_voltage_range = (0.0, 1.0),
        fc = 0.0,
        base_power = 100.0,
        ext = Dict{String, Any}(),
    )
end

#If function exists in PSY, overload it here. 
PSY.get_name(value::TerminalDataSurrogate) = value.name
get_τ(value::TerminalDataSurrogate) = value.τ
get_window_size(value::TerminalDataSurrogate) = value.window_size
get_steadystate_offset_correction(value::TerminalDataSurrogate) =
    value.steadystate_offset_correction
get_trained_voltage_range(value::TerminalDataSurrogate) = value.trained_voltage_range
get_fc(value::TerminalDataSurrogate) = value.fc
PSY.get_base_power(value::TerminalDataSurrogate) = value.base_power
PSY.get_states(value::TerminalDataSurrogate) = value.states
PSY.get_n_states(value::TerminalDataSurrogate) = value.n_states
PSY.get_states_types(value::TerminalDataSurrogate) = value.states_types
PSY.get_ext(value::TerminalDataSurrogate) = value.ext
PSY.get_internal(value::TerminalDataSurrogate) = value.internal

PSY.set_base_power!(value::TerminalDataSurrogate, val) = value.base_power = val
PSY.set_ext!(value::TerminalDataSurrogate, val) = value.ext = val
set_fc!(value::TerminalDataSurrogate, val) = value.fc = val
