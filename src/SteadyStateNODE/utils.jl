
function get_SteadyStateNODE_states(dim_r::Int64)
    states = Symbol[]
    for i in 1:dim_r
        push!(states,Symbol(string("r",i)))
    end 
    return states, dim_r 
end

function activation_map(activation)
    d = Dict(
        "tanh" => tanh,
    )
    return d[activation]
end

function initializer_psids(SSNODE, x)
    layers = get_initializer_structure(SSNODE)
    p =  get_initializer_parameters(SSNODE)
    x_scale =  get_x_scale(SSNODE)
    x_bias = get_x_bias(SSNODE)
    x = x .* x_scale .+ x_bias 
    p_index = 0 
    for l in layers 
        input_dim = l[1]
        output_dim = l[2]
        bias = l[3]
        activation = l[4]
        x = reshape(p[p_index+1:p_index+input_dim*output_dim], (input_dim,output_dim)) * x
        p_index += input_dim*output_dim
        if bias
            x += p[p_index+1:p_index+output_dim]
            p_index += output_dim
        end 
        x = activation_map(activation).(x)
    end 
    return x
end 