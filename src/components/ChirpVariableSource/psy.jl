"""
    mutable struct ChirpVariableSource <: DynamicInjection
        name::String
        R_th::Float64
        X_th::Float64
        ω1::Float64
        ω2::Float64
        tstart::Float64
        N::Float64
        V_amp::Float64
        θ_amp::Float64
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
- `θ_amp::Float64}`: amplitude of angle perturbation
- `base_power::Float64`: Base power
- `states::Vector{Symbol}`: State for time, voltage and angle
- `n_states::Int`
- `ext::Dict{String, Any}`
- `internal::InfrastructureSystemsInternal`: power system internal reference, do not modify
"""
mutable struct ChirpVariableSource <: PSY.DynamicInjection
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
    θ_amp::Float64
    "Base power"
    base_power::Float64
    "State for time, voltage and angle"
    states::Vector{Symbol}
    n_states::Int
    ext::Dict{String, Any}
    "power system internal reference, do not modify"
    internal::PSY.InfrastructureSystemsInternal
end

function ChirpVariableSource(
    name,
    R_th,
    X_th,
    ω1 = 1.0,
    ω2 = 2.0,
    tstart = 1000.0,
    N = 100.0,
    V_amp = 0.0,
    θ_amp = 0.0,
    base_power = 100.0,
    ext = Dict{String, Any}(),
)
    ChirpVariableSource(
        name,
        R_th,
        X_th,
        ω1,
        ω2,
        tstart,
        N,
        V_amp,
        θ_amp,
        base_power,
        ext,
        [:Vt, :θt],
        2,
        PSY.InfrastructureSystemsInternal(),
    )
end

function ChirpVariableSource(;
    name,
    R_th,
    X_th,
    ω1 = 1.0,
    ω2 = 2.0,
    tstart = 1000.0,
    N = 100.0,
    V_amp = 0.0,
    θ_amp = 0.0,
    base_power = 100.0,
    states = [:Vt, :θt],
    n_states = 2,
    ext = Dict{String, Any}(),
    internal = PSY.InfrastructureSystemsInternal(),
)
    ChirpVariableSource(
        name,
        R_th,
        X_th,
        ω1,
        ω2,
        tstart,
        N,
        V_amp,
        θ_amp,
        base_power,
        states,
        n_states,
        ext,
        internal,
    )
end

# Constructor for demo purposes; non-functional.
function ChirpVariableSource(::Nothing)
    ChirpVariableSource(;
        name = "init",
        R_th = 0,
        X_th = 0,
        ω1 = 1.0,
        ω2 = 2.0,
        tstart = 1,
        N = 10,
        V_amp = 0.1,
        θ_amp = 0.05,
        base_power = 0,
        ext = Dict{String, Any}(),
    )
end

"""Get [`ChirpVariableSource`](@ref) `name`."""
PSY.get_name(value::ChirpVariableSource) = value.name
"""Get [`ChirpVariableSource`](@ref) `R_th`."""
PSY.get_R_th(value::ChirpVariableSource) = value.R_th
"""Get [`ChirpVariableSource`](@ref) `X_th`."""
PSY.get_X_th(value::ChirpVariableSource) = value.X_th
"""Get [`ChirpVariableSource`](@ref) `internal_voltage_bias`."""
get_ω1(value::ChirpVariableSource) = value.ω1
"""Get [`ChirpVariableSource`](@ref) `internal_voltage_frequencies`."""
get_ω2(value::ChirpVariableSource) = value.ω2
"""Get [`ChirpVariableSource`](@ref) `internal_voltage_coefficients`."""
get_tstart(value::ChirpVariableSource) = value.tstart
"""Get [`ChirpVariableSource`](@ref) `internal_angle_bias`."""
get_N(value::ChirpVariableSource) = value.N
"""Get [`ChirpVariableSource`](@ref) `internal_angle_frequencies`."""
get_V_amp(value::ChirpVariableSource) = value.V_amp
"""Get [`ChirpVariableSource`](@ref) `internal_angle_coefficients`."""
get_θ_amp(value::ChirpVariableSource) = value.θ_amp
"""Get [`ChirpVariableSource`](@ref) `base_power`."""
PSY.get_base_power(value::ChirpVariableSource) = value.base_power
"""Get [`ChirpVariableSource`](@ref) `states`."""
PSY.get_states(value::ChirpVariableSource) = value.states
"""Get [`ChirpVariableSource`](@ref) `n_states`."""
PSY.get_n_states(value::ChirpVariableSource) = value.n_states
"""Get [`ChirpVariableSource`](@ref) `ext`."""
PSY.get_ext(value::ChirpVariableSource) = value.ext
"""Get [`ChirpVariableSource`](@ref) `internal`."""
PSY.get_internal(value::ChirpVariableSource) = value.internal

"""Set [`ChirpVariableSource`](@ref) `R_th`."""
PSY.set_R_th!(value::ChirpVariableSource, val) = value.R_th = val
"""Set [`ChirpVariableSource`](@ref) `X_th`."""
PSY.set_X_th!(value::ChirpVariableSource, val) = value.X_th = val
"""Set [`ChirpVariableSource`](@ref) `internal_voltage_bias`."""
set_ω1!(value::ChirpVariableSource, val) = value.ω1 = val
"""Set [`ChirpVariableSource`](@ref) `internal_voltage_frequencies`."""
set_ω2!(value::ChirpVariableSource, val) = value.ω2 = val
"""Set [`ChirpVariableSource`](@ref) `internal_voltage_coefficients`."""
set_tstart!(value::ChirpVariableSource, val) = value.tstart = val
"""Set [`ChirpVariableSource`](@ref) `internal_angle_bias`."""
set_N!(value::ChirpVariableSource, val) = value.N = val
"""Set [`ChirpVariableSource`](@ref) `internal_angle_frequencies`."""
set_V_amp!(value::ChirpVariableSource, val) = value.V_amp = val
"""Set [`ChirpVariableSource`](@ref) `internal_angle_coefficients`."""
set_θ_amp!(value::ChirpVariableSource, val) = value.θ_amp = val
"""Set [`ChirpVariableSource`](@ref) `base_power`."""
PSY.set_base_power!(value::ChirpVariableSource, val) = value.base_power = val
"""Set [`ChirpVariableSource`](@ref) `ext`."""
PSY.set_ext!(value::ChirpVariableSource, val) = value.ext = val
