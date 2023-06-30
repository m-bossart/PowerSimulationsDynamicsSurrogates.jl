"""
    mutable struct FrequencyChirpVariableSource <: DynamicInjection
        name::String
        R_th::Float64
        X_th::Float64
        ω1::Float64
        ω2::Float64
        tstart::Float64
        N::Float64
        V_amp::Float64
        ω_amp::Float64
        base_power::Float64
        states::Vector{Symbol}
        n_states::Int
        ext::Dict{String, Any}
        internal::InfrastructureSystemsInternal
    end
This struct acts as an infinity bus with time varying phasor values magnitude and angle V(t) 	heta(t). Time varying functions are represented using fourier series
# Arguments
- `name::String`
- `R_th::Float64`: Source Thevenin resistance, validation range: `(0, nothing)`
- `X_th::Float64`: Source Thevenin reactance, validation range: `(0, nothing)`
- `ω1::Float64`: beginning frequency
- `ω2::Float64`: final frequency
- `tstart::Float64`: chirp start time
- `N::Float64`: chirp duration
- `V_amp::Float64`: amplitude of magnitude perturbation
- `ω_amp::Float64}`: amplitude of frequency perturbation
- `base_power::Float64`: Base power
- `states::Vector{Symbol}`: State for time, voltage and angle
- `n_states::Int`
- `ext::Dict{String, Any}`
- `internal::InfrastructureSystemsInternal`: power system internal reference, do not modify
"""
mutable struct FrequencyChirpVariableSource <: PSY.DynamicInjection
    name::String
    "Source Thevenin resistance"
    R_th::Float64
    "Source Thevenin reactance"
    X_th::Float64
    "Beginning frequency in radians/s"
    ω1::Float64
    "Ending frequency in radians/s"
    ω2::Float64
    "chirp start time in seconds"
    tstart::Float64
    "chirp duration"
    N::Float64
    "Voltage magnitude ampltude perturbation"
    V_amp::Float64
    "Voltage angle ampltude perturbation"
    ω_amp::Float64
    "Base power"
    base_power::Float64
    "State for time, voltage and angle"
    states::Vector{Symbol}
    n_states::Int
    ext::Dict{String, Any}
    "power system internal reference, do not modify"
    internal::PSY.InfrastructureSystemsInternal
end

function FrequencyChirpVariableSource(
    name,
    R_th,
    X_th,
    ω1 = 1.0,
    ω2 = 2.0,
    tstart = 1000.0,
    N = 100.0,
    V_amp = 0.0,
    ω_amp = 0.0,
    base_power = 100.0,
    ext = Dict{String, Any}(),
)
    FrequencyChirpVariableSource(
        name,
        R_th,
        X_th,
        ω1,
        ω2,
        tstart,
        N,
        V_amp,
        ω_amp,
        base_power,
        ext,
        [:Vt, :θt, :ω],
        3,
        PSY.InfrastructureSystemsInternal(),
    )
end

function FrequencyChirpVariableSource(;
    name,
    R_th,
    X_th,
    ω1 = 1.0,
    ω2 = 2.0,
    tstart = 1000.0,
    N = 100.0,
    V_amp = 0.0,
    ω_amp = 0.0,
    base_power = 100.0,
    states = [:Vt, :θt, :ω],
    n_states = 3,
    ext = Dict{String, Any}(),
    internal = PSY.InfrastructureSystemsInternal(),
)
    FrequencyChirpVariableSource(
        name,
        R_th,
        X_th,
        ω1,
        ω2,
        tstart,
        N,
        V_amp,
        ω_amp,
        base_power,
        states,
        n_states,
        ext,
        internal,
    )
end

# Constructor for demo purposes; non-functional.
function FrequencyChirpVariableSource(::Nothing)
    FrequencyChirpVariableSource(;
        name = "init",
        R_th = 0,
        X_th = 0,
        ω1 = 1.0,
        ω2 = 2.0,
        tstart = 1,
        N = 10,
        V_amp = 0.1,
        ω_amp = 0.05,
        base_power = 0,
        ext = Dict{String, Any}(),
    )
end

"""Get [`FrequencyChirpVariableSource`](@ref) `name`."""
PSY.get_name(value::FrequencyChirpVariableSource) = value.name
"""Get [`FrequencyChirpVariableSource`](@ref) `R_th`."""
PSY.get_R_th(value::FrequencyChirpVariableSource) = value.R_th
"""Get [`FrequencyChirpVariableSource`](@ref) `X_th`."""
PSY.get_X_th(value::FrequencyChirpVariableSource) = value.X_th
"""Get [`FrequencyChirpVariableSource`](@ref) `internal_voltage_bias`."""
get_ω1(value::FrequencyChirpVariableSource) = value.ω1
"""Get [`FrequencyChirpVariableSource`](@ref) `internal_voltage_frequencies`."""
get_ω2(value::FrequencyChirpVariableSource) = value.ω2
"""Get [`FrequencyChirpVariableSource`](@ref) `internal_voltage_coefficients`."""
get_tstart(value::FrequencyChirpVariableSource) = value.tstart
"""Get [`FrequencyChirpVariableSource`](@ref) `internal_angle_bias`."""
get_N(value::FrequencyChirpVariableSource) = value.N
"""Get [`FrequencyChirpVariableSource`](@ref) `internal_angle_frequencies`."""
get_V_amp(value::FrequencyChirpVariableSource) = value.V_amp
"""Get [`FrequencyChirpVariableSource`](@ref) `internal_angle_coefficients`."""
get_ω_amp(value::FrequencyChirpVariableSource) = value.ω_amp
"""Get [`FrequencyChirpVariableSource`](@ref) `base_power`."""
PSY.get_base_power(value::FrequencyChirpVariableSource) = value.base_power
"""Get [`FrequencyChirpVariableSource`](@ref) `states`."""
PSY.get_states(value::FrequencyChirpVariableSource) = value.states
"""Get [`FrequencyChirpVariableSource`](@ref) `n_states`."""
PSY.get_n_states(value::FrequencyChirpVariableSource) = value.n_states
"""Get [`FrequencyChirpVariableSource`](@ref) `ext`."""
PSY.get_ext(value::FrequencyChirpVariableSource) = value.ext
"""Get [`FrequencyChirpVariableSource`](@ref) `internal`."""
PSY.get_internal(value::FrequencyChirpVariableSource) = value.internal

"""Set [`FrequencyChirpVariableSource`](@ref) `R_th`."""
PSY.set_R_th!(value::FrequencyChirpVariableSource, val) = value.R_th = val
"""Set [`FrequencyChirpVariableSource`](@ref) `X_th`."""
PSY.set_X_th!(value::FrequencyChirpVariableSource, val) = value.X_th = val
"""Set [`FrequencyChirpVariableSource`](@ref) `internal_voltage_bias`."""
set_ω1!(value::FrequencyChirpVariableSource, val) = value.ω1 = val
"""Set [`FrequencyChirpVariableSource`](@ref) `internal_voltage_frequencies`."""
set_ω2!(value::FrequencyChirpVariableSource, val) = value.ω2 = val
"""Set [`FrequencyChirpVariableSource`](@ref) `internal_voltage_coefficients`."""
set_tstart!(value::FrequencyChirpVariableSource, val) = value.tstart = val
"""Set [`FrequencyChirpVariableSource`](@ref) `internal_angle_bias`."""
set_N!(value::FrequencyChirpVariableSource, val) = value.N = val
"""Set [`FrequencyChirpVariableSource`](@ref) `internal_angle_frequencies`."""
set_V_amp!(value::FrequencyChirpVariableSource, val) = value.V_amp = val
"""Set [`FrequencyChirpVariableSource`](@ref) `internal_angle_coefficients`."""
set_ω_amp!(value::FrequencyChirpVariableSource, val) = value.ω_amp = val
"""Set [`FrequencyChirpVariableSource`](@ref) `base_power`."""
PSY.set_base_power!(value::FrequencyChirpVariableSource, val) = value.base_power = val
"""Set [`FrequencyChirpVariableSource`](@ref) `ext`."""
PSY.set_ext!(value::FrequencyChirpVariableSource, val) = value.ext = val
