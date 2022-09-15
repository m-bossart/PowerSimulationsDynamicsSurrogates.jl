
function get_SteadyStateNODE_states(dim_r::Int64)
    states = Symbol[]
    for i in 1:dim_r
        push!(states, Symbol(string("r", i)))
    end
    return states, dim_r
end

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

PSID.get_inner_vars_count(::SteadyStateNODE) = 0
function PSID._get_frequency_state(d::PSID.DynamicWrapper{SteadyStateNODE})
    return 0
end

function PSID.device_mass_matrix_entries!(
    mass_matrix::AbstractArray,
    dynamic_device::PSID.DynamicWrapper{SteadyStateNODE},
)
    global_index = PSID.get_global_index(dynamic_device)
    mass_matrix_ssnode_entries!(mass_matrix, dynamic_device, global_index)
    return
end

function mass_matrix_ssnode_entries!(
    mass_matrix,
    pvs::PSID.DynamicWrapper{SteadyStateNODE},
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    @debug "Using default mass matrix entries $pvs"
end

function PSID.device!(
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    ::AbstractArray{T},
    ::AbstractArray{T},
    dynamic_device::PSID.DynamicWrapper{SteadyStateNODE},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES}
    θ = dynamic_device.ext["θ0"]
    vd, vq =  PSID.ri_dq(θ) * [voltage_r, voltage_i]
    v_scaled = _input_scale(dynamic_device, [vd, vq])
    refs = dynamic_device.ext["refs"]
    output_ode .= _forward_pass_node(dynamic_device, device_states, v_scaled, refs)
    id_scale =  _forward_pass_observer(dynamic_device, device_states)[1]
    iq_scale = _forward_pass_observer(dynamic_device, device_states)[2]
    id, iq = _target_scale_inverse(dynamic_device, [id_scale, iq_scale])
    ir, ii = PSID.dq_ri(θ) * [id, iq]
    current_r[1] += ir 
    current_i[1] += ii 
    return
end

function PSID.initialize_dynamic_device!(
    dynamic_device::PSID.DynamicWrapper{SteadyStateNODE},
    source::PSY.Source,
    ::AbstractVector,
)
    n_states = PSY.get_n_states(dynamic_device)
    device_states = zeros(n_states)

    #PowerFlow Data
    P0 = PSY.get_active_power(source)
    Q0 = PSY.get_reactive_power(source)
    S0 = P0 + Q0 * 1im
    Vm = PSY.get_magnitude(PSY.get_bus(source))
    θ = PSY.get_angle(PSY.get_bus(source))
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)
    I_R = real(I)
    I_I = imag(I)

    Vd, Vq = PSID.ri_dq(θ) * [V_R, V_I] #Vd should be zero by definition 
    Id, Iq = PSID.ri_dq(θ) * [I_R, I_I]
    _, Vq_scale = _input_scale(dynamic_device, [Vd, Vq])
    Id_scale, Iq_scale = _input_scale(dynamic_device, [Id, Iq])
    function f!(out, x)
        input_scaled = _input_scale(dynamic_device, [Vd, Vq])
        out[1:n_states] .= _forward_pass_node(
            dynamic_device,
            x[1:n_states],
            input_scaled,
            x[(n_states + 1):(n_states + 2)],
        )
        out[n_states + 1] = _forward_pass_observer(dynamic_device, x[1:n_states])[1] - Id_scale
        out[n_states + 2] = _forward_pass_observer(dynamic_device, x[1:n_states])[2] - Iq_scale
    end
    x0 = _forward_pass_initializer(dynamic_device, Vq_scale,  [Id_scale, Iq_scale])
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in SteadyStateNODE failed")
    end
    device_states = sol.zero[1:n_states]
    refs = sol.zero[(n_states + 1):(n_states + 2)]
    dynamic_device.ext["θ0"] = θ
    dynamic_device.ext["refs"] = refs
    dynamic_device.ext["initializer_error"] = x0 .- sol.zero
    return device_states
end

function PSID.DynamicWrapper(
    device::PSY.Source,
    dynamic_device::SteadyStateNODE,
    bus_ix::Int,
    ix_range,
    ode_range,
    inner_var_range,
    sys_base_power,
    sys_base_freq,
)
    device_states = PSY.get_states(dynamic_device)
    ext = _allocate_weights_and_biases(dynamic_device)

    return PSID.DynamicWrapper(
        dynamic_device,
        sys_base_power,
        sys_base_freq,
        PSY.Source,
        PSID.BUS_MAP[PSY.get_bustype(PSY.get_bus(device))],
        Base.Ref(1.0),
        Base.Ref(0.0),
        Base.Ref(0.0),
        Base.Ref(0.0),
        Base.Ref(0.0),
        collect(inner_var_range),
        collect(ix_range),
        collect(ode_range),
        bus_ix,
        Base.ImmutableDict(Dict(device_states .=> ix_range)...),
        Base.ImmutableDict{Int, Vector{Int}}(),
        Base.ImmutableDict{Int, Vector{Int}}(),
        ext,
    )
end

function _allocate_weights_and_biases(dynamic_device::SteadyStateNODE)
    ext = Dict{String, Any}()

    #Allocate weights and biases in initializer
    Ws = []
    bs = []
    layers = get_initializer_structure(dynamic_device)
    p = get_initializer_parameters(dynamic_device)
    p_index_start = 0
    _ = _push_layer_weights_and_biases!(Ws, bs, layers, p, p_index_start)
    ext["initializer"] = Dict{String, Vector{AbstractArray}}("W" => Ws, "b" => bs)

    #Allocate weights and biases in node
    Ws = []
    bs = []
    layers = get_node_structure(dynamic_device)
    p = get_node_parameters(dynamic_device)
    p_index_start = 0
    _ = _push_layer_weights_and_biases!(Ws, bs, layers, p, p_index_start)
    ext["node"] = Dict{String, Vector{AbstractArray}}("W" => Ws, "b" => bs)

    #Allocate weights and biases in observer
    Ws = []
    bs = []
    layers = get_observer_structure(dynamic_device)
    p = get_observer_parameters(dynamic_device)
    p_index_start = 0
    _ = _push_layer_weights_and_biases!(Ws, bs, layers, p, p_index_start)
    ext["observer"] = Dict{String, Vector{AbstractArray}}("W" => Ws, "b" => bs)

    return ext
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

get_W_initializer(wrapper::PSID.DynamicWrapper{SteadyStateNODE}) =
    wrapper.ext["initializer"]["W"]
get_b_initializer(wrapper::PSID.DynamicWrapper{SteadyStateNODE}) =
    wrapper.ext["initializer"]["b"]

get_W_node(wrapper::PSID.DynamicWrapper{SteadyStateNODE}) = wrapper.ext["node"]["W"]
get_b_node(wrapper::PSID.DynamicWrapper{SteadyStateNODE}) = wrapper.ext["node"]["b"]

get_W_observer(wrapper::PSID.DynamicWrapper{SteadyStateNODE}) = wrapper.ext["observer"]["W"]
get_b_observer(wrapper::PSID.DynamicWrapper{SteadyStateNODE}) = wrapper.ext["observer"]["b"]

function _input_scale(wrapper::PSID.DynamicWrapper{SteadyStateNODE}, x)
    xmin = get_input_min(wrapper.device)
    xmax = get_input_max(wrapper.device)
    l, u = get_input_lims(wrapper.device)
    return min_max_normalization(x, xmin, xmax, u, l)
end

function _target_scale(wrapper::PSID.DynamicWrapper{SteadyStateNODE}, x)
    xmin = get_target_min(wrapper.device)
    xmax = get_target_max(wrapper.device)
    l, u = get_target_lims(wrapper.device)
    return min_max_normalization(x, xmin, xmax, u, l)
end

function _target_scale_inverse(wrapper::PSID.DynamicWrapper{SteadyStateNODE}, x)
    xmin = get_target_min(wrapper.device)
    xmax = get_target_max(wrapper.device)
    l, u = get_target_lims(wrapper.device)
    return min_max_normalization_inverse(x, xmin, xmax, u, l)
end

function min_max_normalization(x, xmin, xmax, u, l)  #https://www.baeldung.com/cs/normalizing-inputs-artificial-neural-network
    x_prime = (x .- xmin)./(xmax .- xmin) .* (u .- l) .+ l 
    return x_prime 
end  

function min_max_normalization_inverse(x_prime, xmin, xmax, u, l)
    x = (x_prime .- l) .* (xmax .- xmin) ./ (u .- l) .+ xmin
    return x 
end 

function _forward_pass_initializer(wrapper::PSID.DynamicWrapper{SteadyStateNODE}, input, target)
    x = vcat(input, target)
    W_index = 1
    b_index = 1
    W = get_W_initializer(wrapper)
    b = get_b_initializer(wrapper)
    for layer in get_initializer_structure(wrapper.device)
        x = W[W_index] * x
        W_index += 1
        if layer[3] == true
            x += b[b_index]
            b_index += 1
        end
        x = activation_map(layer[4]).(x)
    end
    return x
end

function _forward_pass_node(wrapper::PSID.DynamicWrapper{SteadyStateNODE}, r, ex, refs)
    input = vcat(r, ex, refs)
    W_index = 1
    b_index = 1
    W = get_W_node(wrapper)
    b = get_b_node(wrapper)
    for layer in get_node_structure(wrapper.device)
        input = W[W_index] * input
        W_index += 1
        if layer[3] == true
            input += b[b_index]
            b_index += 1
        end
        input = activation_map(layer[4]).(input)
    end
    return input
end

function _forward_pass_observer(wrapper::PSID.DynamicWrapper{SteadyStateNODE}, r)
    W_index = 1
    b_index = 1
    W = get_W_observer(wrapper)
    b = get_b_observer(wrapper)
    for layer in get_observer_structure(wrapper.device)
        r = W[W_index] * r
        W_index += 1
        if layer[3] == true
            r += b[b_index]
            b_index += 1
        end
        r = activation_map(layer[4]).(r)
    end
    return r
end
