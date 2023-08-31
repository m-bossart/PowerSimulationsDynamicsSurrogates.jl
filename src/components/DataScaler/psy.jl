abstract type DataScaler end

#################################################
################### MinMaxScaler ################
#################################################

mutable struct MinMaxScaler <: DataScaler
    scale_input::Bool
    input_min::Vector{Float64}
    input_max::Vector{Float64}
    input_lims::Tuple{Float64, Float64}
    scale_target::Bool
    target_min::Vector{Float64}
    target_max::Vector{Float64}
    target_lims::Tuple{Float64, Float64}
end

function MinMaxScaler(;
    scale_input,
    input_min,
    input_max,
    input_lims,
    scale_target,
    target_min,
    target_max,
    target_lims,
)
    MinMaxScaler(
        scale_input,
        input_min,
        input_max,
        input_lims,
        scale_target,
        target_min,
        target_max,
        target_lims,
    )
end

# Constructor for demo purposes; non-functional.
function MinMaxScaler(::Nothing)
    MinMaxScaler(;
        scale_input = true,
        input_min = Float64[],
        input_max = Float64[],
        input_lims = (0.0, 0.0),
        scale_target = true,
        target_min = Float64[],
        target_max = Float64[],
        target_lims = (0.0, 0.0),
    )
end

get_scale_input(value::MinMaxScaler) = value.scale_input
get_input_min(value::MinMaxScaler) = value.input_min
get_input_max(value::MinMaxScaler) = value.input_max
get_input_lims(value::MinMaxScaler) = value.input_lims
get_scale_target(value::MinMaxScaler) = value.scale_target
get_target_min(value::MinMaxScaler) = value.target_min
get_target_max(value::MinMaxScaler) = value.target_max
get_target_lims(value::MinMaxScaler) = value.target_lims

#################################################
################ StandardScaler #################
#################################################

mutable struct StandardScaler <: DataScaler
    scale_input::Bool
    input_mean::Vector{Float64}
    input_std::Vector{Float64}
    scale_target::Bool
    target_mean::Vector{Float64}
    target_std::Vector{Float64}
end

function StandardScaler(;
    scale_input,
    input_mean,
    input_std,
    scale_target,
    target_mean,
    target_std,
)
    StandardScaler(
        scale_input,
        input_mean,
        input_std,
        scale_target,
        target_mean,
        target_std,
    )
end

# Constructor for demo purposes; non-functional.
function StandardScaler(::Nothing)
    StandardScaler(;
        scale_input = true,
        input_mean = Float64[],
        input_std = Float64[],
        scale_target = true,
        target_mean = Float64[],
        target_std = Float64[],
    )
end

get_scale_input(value::StandardScaler) = value.scale_input
get_input_mean(value::StandardScaler) = value.input_mean
get_input_std(value::StandardScaler) = value.input_std
get_scale_target(value::StandardScaler) = value.scale_target
get_target_mean(value::StandardScaler) = value.target_mean
get_target_std(value::StandardScaler) = value.target_std
