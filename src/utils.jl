function hardtanh(x)
    max(-1, min(1, x))
end

function relu(x)
    max(0, x)
end

function activation_map(activation)
    d = Dict("tanh" => tanh, "hardtanh" => hardtanh, "relu" => relu, "identity" => (x) -> x)
    return d[activation]
end

function _push_layer_weights_and_biases!(Ws, bs, layers, p, p_index_start)
    p_index = p_index_start
    for l in layers
        input_dim = l[1]
        output_dim = l[2]
        bias = l[3]
        W = reshape(
            p[(p_index + 1):(p_index + input_dim * output_dim)],
            (output_dim, input_dim),
        )
        p_index += input_dim * output_dim
        push!(Ws, W)
        if bias
            b = p[(p_index + 1):(p_index + output_dim)]
            p_index += output_dim
            push!(bs, b)
        end
    end
    return p_index
end

function min_max_normalization(x, xmin, xmax, u, l)  #https://www.baeldung.com/cs/normalizing-inputs-artificial-neural-network
    x_prime = (x .- xmin) ./ (xmax .- xmin) .* (u .- l) .+ l
    return x_prime
end

function min_max_normalization_inverse(x_prime, xmin, xmax, u, l)
    x = (x_prime .- l) .* (xmax .- xmin) ./ (u .- l) .+ xmin
    return x
end
