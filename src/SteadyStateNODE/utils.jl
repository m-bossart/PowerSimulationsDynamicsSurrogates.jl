
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
    output_ode .= _forward_pass_node(dynamic_device, device_states, [voltage_r, voltage_i])
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
    Vm = PSY.get_magnitude(PSY.get_bus(source))
    θ = PSY.get_angle(PSY.get_bus(source))
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    function f!(out, x) #TODO - explore the parameters of NLsolve, check in PSID
        out = _forward_pass_node(dynamic_device, x, [V_R, V_I])
    end
    x0 = _forward_pass_initializer(dynamic_device, [P0, Q0, Vm, θ])
    display(x0)
    display(_forward_pass_node(dynamic_device, x0, [V_R, V_I]))
    sol = NLsolve.nlsolve(f!, x0)
    display(sol)
    if !NLsolve.converged(sol)
        @warn("Initialization in SteadyStateNODE failed")
    end
    device_states = sol.zero
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

    return PSID.DynamicWrapper{typeof(dynamic_device)}(
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
    layers_exogenous = get_node_structure_exogenous(dynamic_device)
    layers_states = get_node_structure_states(dynamic_device)
    layers_common = get_node_structure_common(dynamic_device)
    p = get_node_parameters(dynamic_device)
    p_index_start = 0
    p_index_intermediate =
        _push_layer_weights_and_biases!(Ws, bs, layers_exogenous, p, p_index_start)
    p_index_intermediate =
        _push_layer_weights_and_biases!(Ws, bs, layers_states, p, p_index_intermediate)
    _ = _push_layer_weights_and_biases!(Ws, bs, layers_common, p, p_index_intermediate)
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

function _forward_pass_node(wrapper::PSID.DynamicWrapper{SteadyStateNODE}, r, ex)
    r_start = r
    W_index = 1
    b_index = 1
    W = get_W_node(wrapper)
    b = get_b_node(wrapper)
    for layer in get_node_structure_exogenous(wrapper.device)
        ex = W[W_index] * ex
        W_index += 1
        if layer[3] == true
            ex += b[b_index]
            b_index += 1
        end
        ex = activation_map(layer[4]).(ex)
    end
    for layer in get_node_structure_states(wrapper.device)
        r = W[W_index] * r
        W_index += 1
        if layer[3] == true
            r += b[b_index]
            b_index += 1
        end
        r = activation_map(layer[4]).(r)
    end
    r = ex .+ r
    for layer in get_node_structure_common(wrapper.device)
        r = W[W_index] * r
        W_index += 1
        if layer[3] == true
            r += b[b_index]
            b_index += 1
        end
        r = activation_map(layer[4]).(r)
    end
    return r .- r_start
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
