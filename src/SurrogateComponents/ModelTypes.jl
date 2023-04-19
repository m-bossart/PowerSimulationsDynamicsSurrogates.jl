mutable struct DataDrivenParams 
    initializer_layer_type::String
    initializer_n_layer::Int64
    initializer_width_layers_relative_input::Int64
    initializer_activation::String
    dynamic_layer_type::String
    dynamic_hidden_states::Int64
    dynamic_n_layer::Int64
    dynamic_width_layers_relative_input::Int64
    dynamic_activation::String
    dynamic_last_layer_bias::Bool
    observation_layer_type::String
    observation_n_layer::Int64
    observation_width_layers_relative_input::Int64
    observation_activation::String
end

function DataDrivenParams(;
    initializer_layer_type = "dense",
    initializer_n_layer = 0,
    initializer_width_layers_relative_input = 0,
    initializer_activation = "tanh",
    dynamic_layer_type = "dense",
    dynamic_hidden_states = 5,
    dynamic_n_layer = 1,
    dynamic_width_layers_relative_input = 4,
    dynamic_activation = "tanh",
    dynamic_last_layer_bias = false,
    observation_layer_type = "dense",
    observation_n_layer = 0,
    observation_width_layers_relative_input = 0,
    observation_activation = "tanh",
    input_normalization = true,  #depends on the surrounding components, should be user specified
    target_normalization = true, 
    input_ref_frame = true, 
    target_ref_frame = true, 
    ) 

    DataDrivenParams(
        initializer_layer_type,
        initializer_n_layer,
        initializer_width_layers_relative_input,
        initializer_activation,
        dynamic_layer_type,
        dynamic_hidden_states,
        dynamic_n_layer,
        dynamic_width_layers_relative_input,
        dynamic_activation,
        dynamic_last_layer_bias,
        observation_layer_type,
        observation_n_layer,
        observation_width_layers_relative_input,
        observation_activation,
    )
end

function build_data_driven_model(dd::DataDrivenParams, ::Type{SteadyStateNODE}, device_name::String) 
    initializer_input_dim = 3   #Vd is zero based on reference frame
    initializer_hidden_dim = dd.initializer_width_layers_relative_input + initializer_input_dim
    initializer_output_dim = dd.dynamic_hidden_states + 2 #initializer includes two references
    dynamic_input_dim = initializer_output_dim + 2 #input is dimension of 2
    dynamic_hidden_dim = dynamic_input_dim + dd.dynamic_width_layers_relative_input
    dynamic_output_dim = dd.dynamic_hidden_states 

    model_initializer =_instantiate_model_initializer(dd, initializer_input_dim, initializer_hidden_dim, initializer_output_dim)    
    model_dynamic = _instantiate_model_dynamic(dd, dynamic_input_dim, dynamic_hidden_dim, dynamic_output_dim)  

    surr = PowerSimulationsDynamicsSurrogates.SteadyStateNODE(
        name = device_name,
        initializer_structure = model_initializer,
        node_structure = model_dynamic,
        base_power = 100.0,
        ext = Dict{String, Any}(),
    )
    return surr 
end 
#at the start of training, you'll have something like: add_parameters_from_database(PSIDComponent) which will add any fixed parameters that only depend on the database.


function _instantiate_model_initializer(m, input_dim, hidden_dim, output_dim)
    type = m.initializer_layer_type
    n_layer = m.initializer_n_layer
    activation = m.initializer_activation
    if type == "dense"
        vector_layers = []
        if n_layer == 0
            push!(vector_layers, (input_dim, output_dim, true, "identity"))
        else
            push!(vector_layers, (input_dim, hidden_dim, true, activation))
            for i in 1:(n_layer - 1)
                push!(vector_layers, (hidden_dim, hidden_dim, true, activation))
            end
            push!(vector_layers, (hidden_dim, output_dim, true, "identity"))     
        end
    end 
    return vector_layers
end

function _instantiate_model_dynamic(m, input_dim, hidden_dim, output_dim)
    type = m.dynamic_layer_type
    n_layer = m.dynamic_n_layer
    activation = m.dynamic_activation
    vector_layers = []
    if type == "dense"
        if n_layer == 0
            push!(
                vector_layers,
                (input_dim, output_dim, m.dynamic_last_layer_bias, "identity"),
            )
        else
            push!(vector_layers, (input_dim, hidden_dim, true, activation))
            for i in 1:(n_layer - 1)
                push!(vector_layers, (hidden_dim, hidden_dim, true, activation))
            end
            push!(
                vector_layers,
                (hidden_dim, output_dim, m.dynamic_last_layer_bias, "identity"),
            )
        end
    end
    return vector_layers
end
