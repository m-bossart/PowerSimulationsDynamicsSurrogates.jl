function min_max_normalization(x, xmin, xmax, u, l)  #https://www.baeldung.com/cs/normalizing-inputs-artificial-neural-network
    x_prime = (x .- xmin) ./ (xmax .- xmin) .* (u .- l) .+ l
    return x_prime
end

function min_max_normalization_inverse(x_prime, xmin, xmax, u, l)
    x = (x_prime .- l) .* (xmax .- xmin) ./ (u .- l) .+ xmin
    return x
end

function _input_scale(scaler::MinMaxScaler, x, indices)
    if get_scale_input(scaler)
        xmin = get_input_min(scaler)[indices]
        xmax = get_input_max(scaler)[indices]
        l, u = get_input_lims(scaler)
        x_prime = similar(x)
        for (i, row) in enumerate(eachrow(x))
            x_prime[i, :] = min_max_normalization(row, xmin, xmax, u, l)
        end
        return x_prime
    else
        return x
    end
end

function _target_scale_inverse(scaler::MinMaxScaler, x)
    if get_scale_target(scaler)
        xmin = get_target_min(scaler)
        xmax = get_target_max(scaler)
        l, u = get_target_lims(scaler)
        return min_max_normalization_inverse(x, xmin, xmax, u, l)
    else
        return x
    end
end

function standard_normalization(x, xmean, xstd)
    x_prime = (x .- xmean) ./ xstd
    return x_prime
end

function standard_normalization_inverse(x, xmean, xstd)
    x = (x_prime .* xstd) .+ xmean
    return x
end

function _input_scale(scaler::StandardScaler, x, indices)
    if get_scale_input(scaler)
        xmean = get_input_mean(scaler)[indices]
        xstd = get_input_std(scaler)[indices]
        x_prime = similar(x)
        for (i, row) in enumerate(eachrow(x))
            x_prime[i, :] = standard_normalization(row, xmean, xstd)
        end
        return x_prime
    else
        return x
    end
end

function _target_scale_inverse(scaler::StandardScaler, x)
    if get_scale_target(scaler)
        xmean = get_target_mean(scaler)
        xstd = get_target_std(scaler)
        return standard_normalization_inverse(x, xmean, xstd)
    else
        return x
    end
end
