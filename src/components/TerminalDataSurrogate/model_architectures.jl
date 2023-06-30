abstract type DataDrivenModelArchitecture end

mutable struct FullyConnected <: DataDrivenModelArchitecture
    nn_structure::Vector{Tuple{Int64, Int64, Bool, String}}
    nn_parameters::Vector{Float64}
    input_min::Vector{Float64}
    input_max::Vector{Float64}
    input_lims::Tuple{Float64, Float64}
    target_min::Vector{Float64}
    target_max::Vector{Float64}
    target_lims::Tuple{Float64, Float64}
end

function FullyConnected(
    nn_structure,
    nn_parameters,
    input_min,
    input_max,
    input_lims,
    target_min,
    target_max,
    target_lims,
)
    FullyConnected(
        nn_structure,
        nn_parameters,
        input_min,
        input_max,
        input_lims,
        target_min,
        target_max,
        target_lims,
    )
end
function FullyConnected(;
    nn_structure,
    nn_parameters,
    input_min,
    input_max,
    input_lims,
    target_min,
    target_max,
    target_lims,
)
    FullyConnected(
        nn_structure,
        nn_parameters,
        input_min,
        input_max,
        input_lims,
        target_min,
        target_max,
        target_lims,
    )
end

# Constructor for demo purposes; non-functional.
function FullyConnected(::Nothing)
    FullyConnected(;
        nn_structure = [(0, 0, false, "")],
        nn_parameters = Float64[],
        input_min = Float64[],
        input_max = Float64[],
        input_lims = (0.0, 0.0),
        target_min = Float64[],
        target_max = Float64[],
        target_lims = (0.0, 0.0),
    )
end

"""Get [`FullyConnected`](@ref) `nn_structure`."""
get_nn_structure(value::FullyConnected) = value.nn_structure
"""Get [`FullyConnected`](@ref) `nn_parameters`."""
get_nn_parameters(value::FullyConnected) = value.nn_parameters
"""Get [`FullyConnected`](@ref) `input_min`."""
get_input_min(value::FullyConnected) = value.input_min
"""Get [`FullyConnected`](@ref) `input_max`."""
get_input_max(value::FullyConnected) = value.input_max
"""Get [`FullyConnected`](@ref) `input_lims`."""
get_input_lims(value::FullyConnected) = value.input_lims
"""Get [`FullyConnected`](@ref) `target_min`."""
get_target_min(value::FullyConnected) = value.target_min
"""Get [`FullyConnected`](@ref) `target_max`."""
get_target_max(value::FullyConnected) = value.target_max
"""Get [`FullyConnected`](@ref) `target_lims`."""
get_target_lims(value::FullyConnected) = value.target_lims
"""Set [`FullyConnected`](@ref) `initializer_parameters`."""
set_nn_parameters!(value::FullyConnected, val) = value.initializer_parameters = val
