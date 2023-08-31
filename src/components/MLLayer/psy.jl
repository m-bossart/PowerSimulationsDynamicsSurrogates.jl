abstract type MLLayer end

#################################################
################### FFNN ########################
#################################################
mutable struct FFNN <: MLLayer
    input_dim::Int
    output_dim::Int
    bias::Bool
    activation::String
end

function FFNN(; input_dim, output_dim, bias, activation)
    FFNN(input_dim, output_dim, bias, activation)
end

# Constructor for demo purposes; non-functional.
function FFNN(::Nothing)
    FFNN(; input_dim = 0, output_dim = 0, bias = false, activation = "")
end

get_input_dim(value::FFNN) = value.input_dim
get_output_dim(value::FFNN) = value.output_dim
get_bias(value::FFNN) = value.bias
get_activation(value::FFNN) = value.activation

#################################################
################### RNN #########################
#################################################
mutable struct RNN <: MLLayer
    input_dim::Int
    output_dim::Int
    time_dim::Int
    output_sequences::Bool
    bias::Bool
    activation::String
end

function RNN(; input_dim, output_dim, time_dim, output_sequences, bias, activation)
    RNN(input_dim, output_dim, time_dim, output_sequences, bias, activation)
end

# Constructor for demo purposes; non-functional.
function RNN(::Nothing)
    RNN(;
        input_dim = 0,
        output_dim = 0,
        time_dim = 0,
        output_sequences = false,
        bias = false,
        activation = "",
    )
end

get_input_dim(value::RNN) = value.input_dim
get_output_dim(value::RNN) = value.output_dim
get_time_dim(value::RNN) = value.time_dim
get_output_sequences(value::RNN) = value.output_sequences
get_bias(value::RNN) = value.bias
get_activation(value::RNN) = value.activation

#################################################
################### LSTM #########################
#################################################
mutable struct LSTM <: MLLayer
    input_dim::Int
    output_dim::Int
    time_dim::Int
    output_sequences::Bool
    bias::Bool
    activation::String
end

function LSTM(; input_dim, output_dim, time_dim, output_sequences, bias, activation)
    LSTM(input_dim, output_dim, time_dim, output_sequences, bias, activation)
end

# Constructor for demo purposes; non-functional.
function LSTM(::Nothing)
    LSTM(;
        input_dim = 0,
        output_dim = 0,
        time_dim = 0,
        output_sequences = false,
        bias = false,
        activation = "",
    )
end

get_input_dim(value::LSTM) = value.input_dim
get_output_dim(value::LSTM) = value.output_dim
get_time_dim(value::LSTM) = value.time_dim
get_output_sequences(value::LSTM) = value.output_sequences
get_bias(value::LSTM) = value.bias
get_activation(value::LSTM) = value.activation
