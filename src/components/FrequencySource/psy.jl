"""
    mutable struct FrequencySource <: DynamicInjection
        name::String
        R_th::Float64
        X_th::Float64
        Tv::Float64
        Tω::Float64
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
- `Tv::Float64`: time constant for voltage low pass filter
- `Tω::Float64`: time contstant for frequency low pass filter
- `base_power::Float64`: Base power
- `states::Vector{Symbol}`: State for time, voltage and angle
- `n_states::Int`
- `ext::Dict{String, Any}`
- `internal::InfrastructureSystemsInternal`: power system internal reference, do not modify
"""
mutable struct FrequencySource <: PSY.DynamicInjection
    name::String
    "Source Thevenin resistance"
    R_th::Float64
    "Source Thevenin reactance"
    X_th::Float64
    "Voltage reference"
    V_ref::Float64
    "Frequency reference"
    ω_ref::Float64
    "Time constant for voltage low pass filter"
    Tv::Float64
    "Time constant for frequency low pass filter"
    Tω::Float64
    "Base power"
    base_power::Float64
    "State for time, voltage and angle"
    states::Vector{Symbol}
    n_states::Int
    ext::Dict{String, Any}
    "power system internal reference, do not modify"
    internal::PSY.InfrastructureSystemsInternal
end

function FrequencySource(
    name,
    R_th,
    X_th,
    V_ref,
    ω_ref,
    Tv = 1.0,
    Tω = 1.0,
    base_power = 100.0,
    ext = Dict{String, Any}(),
)
    FrequencySource(
        name,
        R_th,
        X_th,
        V_ref,
        ω_ref,
        Tv,
        Tω,
        base_power,
        ext,
        [:Vt, :θt, :ω],
        3,
        PSY.InfrastructureSystemsInternal(),
    )
end

function FrequencySource(;
    name,
    R_th,
    X_th,
    V_ref,
    ω_ref,
    Tv = 1.0,
    Tω = 2.0,
    base_power = 100.0,
    states = [:Vt, :θt, :ω],
    n_states = 3,
    ext = Dict{String, Any}(),
    internal = PSY.InfrastructureSystemsInternal(),
)
    FrequencySource(
        name,
        R_th,
        X_th,
        V_ref,
        ω_ref,
        Tv,
        Tω,
        base_power,
        states,
        n_states,
        ext,
        internal,
    )
end

# Constructor for demo purposes; non-functional.
function FrequencySource(::Nothing)
    FrequencySource(;
        name = "init",
        R_th = 0,
        X_th = 0,
        V_ref = 1.0,
        ω_ref = 1.0,
        Tv = 1.0,
        Tω = 1.0,
        base_power = 0,
        ext = Dict{String, Any}(),
    )
end

"""Get [`FrequencySource`](@ref) `name`."""
PSY.get_name(value::FrequencySource) = value.name
"""Get [`FrequencySource`](@ref) `R_th`."""
PSY.get_R_th(value::FrequencySource) = value.R_th
"""Get [`FrequencySource`](@ref) `X_th`."""
PSY.get_X_th(value::FrequencySource) = value.X_th
"""Get [`FrequencySource`](@ref) `V_ref`."""
PSY.get_V_ref(value::FrequencySource) = value.V_ref
"""Get [`FrequencySource`](@ref) `ω_ref`."""
PSY.get_ω_ref(value::FrequencySource) = value.ω_ref
"""Get [`FrequencySource`](@ref) `Tv`."""
PSY.get_Tv(value::FrequencySource) = value.Tv
"""Get [`FrequencySource`](@ref) `Tω`."""
get_Tω(value::FrequencySource) = value.Tω
"""Get [`FrequencySource`](@ref) `base_power`."""
PSY.get_base_power(value::FrequencySource) = value.base_power
"""Get [`FrequencySource`](@ref) `states`."""
PSY.get_states(value::FrequencySource) = value.states
"""Get [`FrequencySource`](@ref) `n_states`."""
PSY.get_n_states(value::FrequencySource) = value.n_states
"""Get [`FrequencySource`](@ref) `ext`."""
PSY.get_ext(value::FrequencySource) = value.ext
"""Get [`FrequencySource`](@ref) `internal`."""
PSY.get_internal(value::FrequencySource) = value.internal

"""Set [`FrequencySource`](@ref) `R_th`."""
PSY.set_R_th!(value::FrequencySource, val) = value.R_th = val
"""Set [`FrequencySource`](@ref) `X_th`."""
PSY.set_X_th!(value::FrequencySource, val) = value.X_th = val
"""Set [`FrequencySource`](@ref) `V_ref`."""
PSY.set_V_ref!(value::FrequencySource, val) = value.V_ref = val
"""Set [`FrequencySource`](@ref) `ω_ref`."""
PSY.set_ω_ref!(value::FrequencySource, val) = value.ω_ref = val
"""Set [`FrequencySource`](@ref) `Tv`."""
PSY.set_Tv!(value::FrequencySource, val) = value.Tv = val
"""Set [`FrequencySource`](@ref) `Tω`."""
set_Tω!(value::FrequencySource, val) = value.Tω = val
"""Set [`FrequencySource`](@ref) `base_power`."""
PSY.set_base_power!(value::FrequencySource, val) = value.base_power = val
"""Set [`FrequencySource`](@ref) `ext`."""
PSY.set_ext!(value::FrequencySource, val) = value.ext = val
