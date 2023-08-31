#= 
NOTE: Need to be consistent in the order in which arrays are allocated to parameter_arrays
and the indexing in the _call_layer function
=#

function _allocate_model_parameters!(
    ext::Dict{String, Any},
    model_architecture,
    model_parameters,
)
    parameter_arrays = []
    p_index = 0
    for layer in model_architecture
        p_index =
            _allocate_layer_parameters!(parameter_arrays, layer, model_parameters, p_index)
    end
    @assert p_index == length(model_parameters)  #make sure all the parameters are allocated to an array 
    ext["nn"] = parameter_arrays
    return
end

function _allocate_layer_parameters!(
    parameter_arrays,
    layer::FFNN,
    model_parameters,
    p_index,
)
    input_dim = get_input_dim(layer)
    output_dim = get_output_dim(layer)
    bias = get_bias(layer)
    W = reshape(
        model_parameters[(p_index + 1):(p_index + input_dim * output_dim)],
        (output_dim, input_dim),
    )
    p_index += input_dim * output_dim
    push!(parameter_arrays, W)
    if bias
        b = model_parameters[(p_index + 1):(p_index + output_dim)]
        p_index += output_dim
        push!(parameter_arrays, b)
    end
    return p_index
end

#order: input kernel, recurrent kernel, bias
function _allocate_layer_parameters!(
    parameter_arrays,
    layer::RNN,
    model_parameters,
    p_index,
)
    input_dim = get_input_dim(layer)
    output_dim = get_output_dim(layer)
    bias = get_bias(layer)
    W = reshape(
        model_parameters[(p_index + 1):(p_index + input_dim * output_dim)],
        (output_dim, input_dim),
    )
    p_index += input_dim * output_dim
    push!(parameter_arrays, W)
    W_recurrent = reshape(
        model_parameters[(p_index + 1):(p_index + output_dim * output_dim)],
        (output_dim, output_dim),
    )
    p_index += output_dim * output_dim
    push!(parameter_arrays, W_recurrent)
    if bias
        b = model_parameters[(p_index + 1):(p_index + output_dim)]
        p_index += output_dim
        push!(parameter_arrays, b)
    end
    return p_index
end

#order: input kernel, recurrent kernel, bias
#forget (f) input (i), output(o), candidate(c)  #Which are the orders 
function _allocate_layer_parameters!(
    parameter_arrays,
    layer::LSTM,
    model_parameters,
    p_index,
)
    input_dim = get_input_dim(layer)
    output_dim = get_output_dim(layer)
    bias = get_bias(layer)
    #input kerlels (f,i,o,c)
    for i in 1:4
        W = reshape(
            model_parameters[(p_index + 1):(p_index + input_dim * output_dim)],
            (output_dim, input_dim),
        )
        p_index += input_dim * output_dim
        push!(parameter_arrays, W)
    end
    for i in 1:4
        W_recurrent = reshape(
            model_parameters[(p_index + 1):(p_index + output_dim * output_dim)],
            (output_dim, output_dim),
        )
        p_index += output_dim * output_dim
        push!(parameter_arrays, W_recurrent)
    end
    if bias
        for i in 1:4
            b = model_parameters[(p_index + 1):(p_index + output_dim)]
            p_index += output_dim
            push!(parameter_arrays, b)
        end
    end
    return p_index
end

function _call_model(parameter_arrays, model_architecture, static_input, dynamic_input)
    output = _reshape_input(static_input, dynamic_input, model_architecture[1])  #either concatonate or pad input based on the first layer type
    ix = 1
    for layer in model_architecture
        output, ix = _call_layer(output, layer, parameter_arrays, ix)
    end
    return output
end

function _call_layer(x, layer::FFNN, parameter_arrays, ix)
    starting_ix = ix
    x = parameter_arrays[starting_ix] * x
    starting_ix += 1
    if get_bias(layer) == true
        x += parameter_arrays[starting_ix]
        starting_ix += 1
    end
    x = activation_map(get_activation(layer)).(x)
    return x, starting_ix
end

function _call_layer(x, layer::RNN, parameter_arrays, ix)
    starting_ix = ix
    output_dim = get_output_dim(layer)
    time_dim = get_time_dim(layer)
    hidden_state = zeros(output_dim)
    for i in 1:time_dim
        hidden_state =
            parameter_arrays[starting_ix] * x[i, :] .+
            parameter_arrays[starting_ix + 1] * hidden_state
        if get_bias(layer) == true
            hidden_state += parameter_arrays[starting_ix + 2]
        end
        hidden_state = activation_map(get_activation(layer)).(hidden_state)
    end
    starting_ix += 3
    return hidden_state, starting_ix
end

#Follows notation here: https://en.wikipedia.org/wiki/Long_short-term_memory
#TODO - check the axis of concatonation of individual kernels
#Order in Keras is i,f,c,o
function _call_layer(x, layer::LSTM, parameter_arrays, ix)
    starting_ix = ix
    @assert get_bias(layer)
    starting_ix = ix
    output_dim = get_output_dim(layer)
    time_dim = get_time_dim(layer)
    hidden_state = zeros(output_dim)
    #Input kernel
    Wi = parameter_arrays[starting_ix]
    Wf = parameter_arrays[starting_ix + 1]
    Wc = parameter_arrays[starting_ix + 2]
    Wo = parameter_arrays[starting_ix + 3]
    #Recurrent kernel
    Ui = parameter_arrays[starting_ix + 4]
    Uf = parameter_arrays[starting_ix + 5]
    Uc = parameter_arrays[starting_ix + 6]
    Uo = parameter_arrays[starting_ix + 7]
    #Bias 
    bi = parameter_arrays[starting_ix + 8]
    bf = parameter_arrays[starting_ix + 9]
    bc = parameter_arrays[starting_ix + 10]
    bo = parameter_arrays[starting_ix + 11]
    starting_ix += 12

    h = zeros(output_dim)
    c = zeros(output_dim)
    for t in 1:time_dim
        f = sigmoid.(Wf * x[t, :] + Uf * h + bf)
        i = sigmoid.(Wi * x[t, :] + Ui * h + bi)
        o = sigmoid.(Wo * x[t, :] + Uo * h + bo)
        c_tilde = activation_map(get_activation(layer)).(Wc * x[t, :] + Uc * h + bc)
        c = f .* c .+ i .* c_tilde
        h = o .* activation_map(get_activation(layer)).(c)
    end
    return h, starting_ix
end

function _reshape_input(static_input, dynamic_input, layer::FFNN)
    static_input = reshape(static_input, length(static_input))
    dynamic_input = reshape(dynamic_input', length(dynamic_input))  #TODO - better implementation of reshape? 
    return vcat(dynamic_input, static_input)
end

function _reshape_input(static_input, dynamic_input, layer::Union{RNN, LSTM})
    dim_time = size(dynamic_input)[1]
    dim_static_features = size(static_input)[2]
    dim_dynamic_features = size(dynamic_input)[2]
    input = Array{Any}(undef, dim_time, dim_static_features + dim_dynamic_features)
    for i in 1:dim_time
        input[i, :] = vcat(dynamic_input[i, :], static_input[1, :])
    end
    return input
end

function parse_h5(h5_file_path)
    model_architecture = []
    model_parameters = Vector{Float64}[]

    fid = HDF5.h5open(h5_file_path, "r")
    layer_names = HDF5.read_attribute(fid["model_weights"], "layer_names")
    layer_configurations =
        JSON3.read(HDF5.read_attribute(fid, "model_config"))[:config][:layers]
    layer_weights = HDF5.h5open(h5_file_path, "r") do file
        read(file, "model_weights")
    end

    delete!(layer_weights, "top_level_model_weights")

    for (name, config, weights) in zip(layer_names, layer_configurations, layer_weights)
        if config["class_name"] == "Dense"
            input_size = size(layer_weights[name][name]["kernel:0"])[2]
            output_size = size(layer_weights[name][name]["kernel:0"])[1]
            @assert size(layer_weights[name][name]["bias:0"])[1] != 0
            model_parameters =
                vcat(model_parameters, vec(layer_weights[name][name]["kernel:0"]))
            model_parameters =
                vcat(model_parameters, vec(layer_weights[name][name]["bias:0"]))
            push!(
                model_architecture,
                FFNN(
                    input_dim = input_size,
                    output_dim = output_size,
                    bias = true,
                    activation = config[:config]["activation"],
                ),
            )
        elseif config["class_name"] == "SimpleRNN"   #Need to find both time dimension and return sequences 
            input_size = size(layer_weights[name][name]["simple_rnn_cell"]["kernel:0"])[2]
            output_size = size(layer_weights[name][name]["simple_rnn_cell"]["kernel:0"])[1]
            @assert size(layer_weights[name][name]["simple_rnn_cell"]["bias:0"])[1] != 0
            model_parameters = vcat(
                model_parameters,
                vec(layer_weights[name][name]["simple_rnn_cell"]["kernel:0"]),
            )
            model_parameters = vcat(
                model_parameters,
                vec(layer_weights[name][name]["simple_rnn_cell"]["recurrent_kernel:0"]),
            )
            model_parameters = vcat(
                model_parameters,
                vec(layer_weights[name][name]["simple_rnn_cell"]["bias:0"]),
            )
            time_dim = layer_configurations[1][:config]["batch_input_shape"][2]
            push!(
                model_architecture,
                RNN(
                    input_dim = input_size,
                    output_dim = output_size,
                    time_dim = time_dim,
                    output_sequences = config[:config]["return_sequences"],
                    bias = true,
                    activation = config[:config]["activation"],
                ),
            )
        elseif config["class_name"] == "LSTM"
            input_size = size(layer_weights[name][name]["lstm_cell"]["kernel:0"])[2]
            output_size = size(layer_weights[name][name]["lstm_cell"]["kernel:0"])[1] / 4       #four different kernels in 1 
            @assert size(layer_weights[name][name]["lstm_cell"]["bias:0"])[1] != 0
            for i in 1:4
                model_parameters = vcat(
                    model_parameters,
                    vec(
                        layer_weights[name][name]["lstm_cell"]["kernel:0"][
                            Int((i - 1) * output_size + 1):Int(i * output_size),
                            :,
                        ],
                    ),
                )
            end
            for i in 1:4
                model_parameters = vcat(
                    model_parameters,
                    vec(
                        layer_weights[name][name]["lstm_cell"]["recurrent_kernel:0"][
                            Int((i - 1) * output_size + 1):Int(i * output_size),
                            :,
                        ],
                    ),
                )    #NOTE: This has 4 independent recurrent kernels in it...
            end
            for i in 1:4
                model_parameters = vcat(
                    model_parameters,
                    vec(
                        layer_weights[name][name]["lstm_cell"]["bias:0"][Int(
                            (i - 1) * output_size + 1,
                        ):Int(i * output_size)],
                    ),
                )
            end
            time_dim = layer_configurations[1][:config]["batch_input_shape"][2]
            push!(
                model_architecture,
                LSTM(
                    input_dim = input_size,
                    output_dim = output_size,
                    time_dim = time_dim,
                    output_sequences = config[:config]["return_sequences"],
                    bias = true,
                    activation = config[:config]["activation"],
                ),
            )
        end
    end
    return model_architecture, model_parameters
end
