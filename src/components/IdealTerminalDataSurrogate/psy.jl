mutable struct IdealTerminalDataSurrogate <: DirectSolutionSurrogate
    name::String
    τ::Float64
    reference_solution::Any
    reference_solution_ir_indices::Vector{Int64}
    reference_solution_ii_indices::Vector{Int64}
    base_power::Float64
    states::Vector{Symbol}
    n_states::Int
    states_types::Vector{PSY.StateTypes}
    ext::Dict{String, Any}
    internal::IS.InfrastructureSystemsInternal
end

function IdealTerminalDataSurrogate(
    name,
    τ = 0.0,
    reference_solution = nothing,
    reference_solution_ir_indices = Int64[],
    reference_solution_ii_indices = Int64[],
    base_power = 100.0,
    ext = Dict{String, Any}(),
)
    IdealTerminalDataSurrogate(
        name,
        τ,
        reference_solution,
        reference_solution_ir_indices,
        reference_solution_ii_indices,
        base_power,
        [:ir, :ii],
        2,
        [PSY.StateTypes.Differential, PSY.StateTypes.Differential],
        ext,
        IS.InfrastructureSystemsInternal(),
    )
end

function IdealTerminalDataSurrogate(;
    name,
    τ = 0.0,
    reference_solution = nothing,
    reference_solution_ir_indices = Int64[],
    reference_solution_ii_indices = Int64[],
    base_power = 100.0,
    states = [:ir, :ii],
    n_states = 2,
    states_types = [PSY.StateTypes.Differential, PSY.StateTypes.Differential],
    ext = Dict{String, Any}(),
    internal = IS.InfrastructureSystemsInternal(),
)
    IdealTerminalDataSurrogate(
        name,
        τ,
        reference_solution,
        reference_solution_ir_indices,
        reference_solution_ii_indices,
        base_power,
        states,
        n_states,
        states_types,
        ext,
        internal,
    )
end

# Constructor for demo purposes; non-functional.
function IdealTerminalDataSurrogate(::Nothing)
    IdealTerminalDataSurrogate(;
        name = "init",
        τ = 0.0,
        reference_solution = nothing,
        reference_solution_ir_indices = Int64[],
        reference_solution_ii_indices = Int64[],
        base_power = 100.0,
        ext = Dict{String, Any}(),
    )
end

#If function exists in PSY, overload it here. 
PSY.get_name(value::IdealTerminalDataSurrogate) = value.name
get_τ(value::IdealTerminalDataSurrogate) = value.τ
get_reference_solution(value::IdealTerminalDataSurrogate) = value.reference_solution
get_reference_solution_ir_indices(value::IdealTerminalDataSurrogate) =
    value.reference_solution_ir_indices
get_reference_solution_ii_indices(value::IdealTerminalDataSurrogate) =
    value.reference_solution_ii_indices
PSY.get_base_power(value::IdealTerminalDataSurrogate) = value.base_power
PSY.get_states(value::IdealTerminalDataSurrogate) = value.states
PSY.get_n_states(value::IdealTerminalDataSurrogate) = value.n_states
PSY.get_states_types(value::IdealTerminalDataSurrogate) = value.states_types
PSY.get_ext(value::IdealTerminalDataSurrogate) = value.ext
PSY.get_internal(value::IdealTerminalDataSurrogate) = value.internal

PSY.set_base_power!(value::IdealTerminalDataSurrogate, val) = value.base_power = val
PSY.set_ext!(value::IdealTerminalDataSurrogate, val) = value.ext = val
