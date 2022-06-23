
function get_SteadyStateNODE_states(dim_r::Int64)
    states = Symbol[]
    for i in 1:dim_r
        push!(states, Symbol(string("r", i)))
    end
    return states, dim_r
end

function activation_map(activation)
    d = Dict("tanh" => tanh)
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
    v_scaled = _exogenous_scale(dynamic_device, [voltage_r, voltage_i])
    refs = dynamic_device.ext["refs"]
    output_ode .= _forward_pass_node(dynamic_device, device_states, v_scaled, refs)
    current_r[1] += _forward_pass_observer(dynamic_device, device_states)[1]
    current_i[1] += _forward_pass_observer(dynamic_device, device_states)[2]
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
    function f!(out, x)
        exogenous_scaled = _exogenous_scale(dynamic_device, [V_R, V_I])
        out[1:n_states] .= _forward_pass_node(
            dynamic_device,
            x[1:n_states],
            exogenous_scaled,
            x[(n_states + 1):(n_states + 2)],
        )
        out[n_states + 1] = _forward_pass_observer(dynamic_device, x[1:n_states])[1] - I_R
        out[n_states + 2] = _forward_pass_observer(dynamic_device, x[1:n_states])[2] - I_I
    end
    x_scaled = _x_scale(dynamic_device, [P0, Q0, Vm, θ])
    x0 = _forward_pass_initializer(dynamic_device, x_scaled)
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in SteadyStateNODE failed")
    end
    device_states = sol.zero[1:n_states]
    refs = sol.zero[(n_states + 1):(n_states + 2)]
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

function _exogenous_scale(wrapper::PSID.DynamicWrapper{SteadyStateNODE}, ex)
    ex .= ex .* get_exogenous_scale(wrapper.device)
    ex .= get_exogenous_bias(wrapper.device) .+ ex
    return ex
end

function _x_scale(wrapper::PSID.DynamicWrapper{SteadyStateNODE}, x)
    x .= x .* get_x_scale(wrapper.device)
    x .= get_x_bias(wrapper.device) .+ x
    return x
end

function _forward_pass_initializer(wrapper::PSID.DynamicWrapper{SteadyStateNODE}, x)
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
