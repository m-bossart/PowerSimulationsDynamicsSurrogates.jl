"""
    mutable struct SteadyStateNODEObs <: DynamicInjection
        name::String
        initializer_structure::Vector{Tuple{Int64, Int64, Bool, String}}
        initializer_parameters::Vector{Float64}
        node_structure::Vector{Tuple{Int64, Int64, Bool, String}}
        node_parameters::Vector{Float64}
        observer_structure::Vector{Tuple{Int64, Int64, Bool, String}}
        observer_parameters::Vector{Float64}
        input_min::Vector{Float64}
        input_max::Vector{Float64}
        input_lims::Tuple{Float64, Float64}
        target_min::Vector{Float64}
        target_max::Vector{Float64}
        target_lims::Tuple{Float64, Float64}
        base_power::Float64
        states::Vector{Symbol}
        n_states::Int
        ext::Dict{String, Any}
        internal::InfrastructureSystemsInternal
    end

Experimental surrogate

# Arguments
- `name::String`
- `initializer_structure::Vector{Tuple{Int64, Int64, Bool, String}}`: layers of the initializer
- `initializer_parameters::Vector{Float64}`: parameters of the initializer
- `node_structure::Vector{Tuple{Int64, Int64, Bool, String}}`: layers of the node 
- `node_parameters::Vector{Float64}`: parameters of the node
- `observer_structure::Vector{Tuple{Int64, Int64, Bool, String}}`: layers of the initializer
- `observer_parameters::Vector{Float64}`: parameters of the node
- `input_min::Vector{Float64}`: minimum values of inputs
- `input_max::Vector{Float64}`: maximum values of inputs
- `input_lims::Tuple{Float64, Float64}`: limits for inputs
- `target_min::Vector{Float64}`: minimum values of targets
- `target_max::Vector{Float64}`: maximum values of targets
- `target_lims::Tuplse{Float64, Float64`: limits for targets
- `base_power::Float64`: Base power
- `states::Vector{Symbol}`: The states of GenericDER depend on the Flags
- `n_states::Int`: The states of GenericDER depend on the Flags
- `ext::Dict{String, Any}`
- `internal::InfrastructureSystemsInternal`: power system internal reference, do not modify
"""
mutable struct SteadyStateNODEObs <: PSY.DynamicInjection
    name::String
    "layers of the initializer"
    initializer_structure::Vector{Tuple{Int64, Int64, Bool, String}}
    "parameters of the initializer"
    initializer_parameters::Vector{Float64}
    "layers of the initializer"
    node_structure::Vector{Tuple{Int64, Int64, Bool, String}}
    "layers of the initializer"
    node_parameters::Vector{Float64}
    "layers of the initializer"
    observer_structure::Vector{Tuple{Int64, Int64, Bool, String}}
    "parameters of the node"
    observer_parameters::Vector{Float64}
    "minimum values of inputs"
    input_min::Vector{Float64}
    "maximum values of inputs"
    input_max::Vector{Float64}
    "limits for inputs"
    input_lims::Tuple{Float64, Float64}
    "minimum values of targets"
    target_min::Vector{Float64}
    "maximum values of targets"
    target_max::Vector{Float64}
    "limits for targets"
    target_lims::Tuple{Float64, Float64}
    "Base power"
    base_power::Float64
    "The states of GenericDER depend on the Flags"
    states::Vector{Symbol}
    "The states of GenericDER depend on the Flags"
    n_states::Int
    ext::Dict{String, Any}
    "power system internal reference, do not modify"
    internal::IS.InfrastructureSystemsInternal
end

function SteadyStateNODEObs(
    name,
    initializer_structure = [(0, 0, true, "init")], #input dim, output dim, bias?, activation 
    initializer_parameters = [0.0],
    node_structure = [(0, 0, true, "init")],
    node_parameters = [0.0],
    observer_structure = [(0, 0, true, "init")],
    observer_parameters = [0.0],
    input_min = [0.0],
    input_max = [0.0],
    input_lims = (-1, 1),
    target_min = [0.0],
    target_max = [0.0],
    target_lims = (-1, 1),
    base_power = 100.0,
    ext = Dict{String, Any}(),
)
    SteadyStateNODEObs(
        name,
        initializer_structure,
        initializer_parameters,
        node_structure,
        node_parameters,
        observer_structure,
        observer_parameters,
        input_min,
        input_max,
        input_lims,
        target_min,
        target_max,
        target_lims,
        base_power,
        ext,
        get_SteadyStateNODEObs_states(initializer_structure[end][2] - 2)[1],
        get_SteadyStateNODEObs_states(initializer_structure[end][2] - 2)[2],
        IS.InfrastructureSystemsInternal(),
    )
end

function SteadyStateNODEObs(;
    name,
    initializer_structure = [(0, 0, true, "init")],
    initializer_parameters = [0.0],
    node_structure = [(0, 0, true, "init")],
    node_parameters = [0.0],
    observer_structure = [(0, 0, true, "init")],
    observer_parameters = [0.0],
    input_min = [0.0],
    input_max = [0.0],
    input_lims = (-1, 1),
    target_min = [0.0],
    target_max = [0.0],
    target_lims = (-1, 1),
    base_power = 100.0,
    states = get_SteadyStateNODEObs_states(initializer_structure[end][2] - 2)[1],
    n_states = get_SteadyStateNODEObs_states(initializer_structure[end][2] - 2)[2],
    ext = Dict{String, Any}(),
    internal = IS.InfrastructureSystemsInternal(),
)
    SteadyStateNODEObs(
        name,
        initializer_structure,
        initializer_parameters,
        node_structure,
        node_parameters,
        observer_structure,
        observer_parameters,
        input_min,
        input_max,
        input_lims,
        target_min,
        target_max,
        target_lims,
        base_power,
        states,
        n_states,
        ext,
        internal,
    )
end

# Constructor for demo purposes; non-functional.
function SteadyStateNODEObs(::Nothing)
    SteadyStateNODEObs(;
        name = "init",
        initializer_structure = [(0, 0, true, "init")],
        initializer_parameters = [0],
        node_structure = [(0, 0, true, "init")],
        node_parameters = [0],
        observer_structure = [(0, 0, true, "init")],
        observer_parameters = Any[0],
        input_min = [0.0],
        input_max = [0.0],
        input_lims = (-1, 1),
        target_min = [0.0],
        target_max = [0.0],
        target_lims = (-1, 1),
        base_power = 0,
        ext = Dict{String, Any}(),
    )
end

#If function exists in PSY, overload it here. 
"""Get [`SteadyStateNODEObs`](@ref) `name`."""
PSY.get_name(value::SteadyStateNODEObs) = value.name
"""Get [`SteadyStateNODEObs`](@ref) `initializer_structure`."""
get_initializer_structure(value::SteadyStateNODEObs) = value.initializer_structure
"""Get [`SteadyStateNODEObs`](@ref) `initializer_parameters`."""
get_initializer_parameters(value::SteadyStateNODEObs) = value.initializer_parameters
"""Get [`SteadyStateNODEObs`](@ref) `node_structure`."""
get_node_structure(value::SteadyStateNODEObs) = value.node_structure
"""Get [`SteadyStateNODEObs`](@ref) `node_structure_states`."""
get_node_parameters(value::SteadyStateNODEObs) = value.node_parameters
"""Get [`SteadyStateNODEObs`](@ref) `observer_structure`."""
get_observer_structure(value::SteadyStateNODEObs) = value.observer_structure
"""Get [`SteadyStateNODEObs`](@ref) `observer_parameters`."""
get_observer_parameters(value::SteadyStateNODEObs) = value.observer_parameters
"""Get [`SteadyStateNODEObs`](@ref) `input_min`."""
get_input_min(value::SteadyStateNODEObs) = value.input_min
"""Get [`SteadyStateNODEObs`](@ref) `input_max`."""
get_input_max(value::SteadyStateNODEObs) = value.input_max
"""Get [`SteadyStateNODEObs`](@ref) `input_lims`."""
get_input_lims(value::SteadyStateNODEObs) = value.input_lims
"""Get [`SteadyStateNODEObs`](@ref) `target_min`."""
get_target_min(value::SteadyStateNODEObs) = value.target_min
"""Get [`SteadyStateNODEObs`](@ref) `target_max`."""
get_target_max(value::SteadyStateNODEObs) = value.target_max
"""Get [`SteadyStateNODEObs`](@ref) `target_lims`."""
get_target_lims(value::SteadyStateNODEObs) = value.target_lims
"""Get [`SteadyStateNODEObs`](@ref) `base_power`."""
PSY.get_base_power(value::SteadyStateNODEObs) = value.base_power
"""Get [`SteadyStateNODEObs`](@ref) `states`."""
PSY.get_states(value::SteadyStateNODEObs) = value.states
"""Get [`SteadyStateNODEObs`](@ref) `n_states`."""
PSY.get_n_states(value::SteadyStateNODEObs) = value.n_states
"""Get [`SteadyStateNODEObs`](@ref) `ext`."""
PSY.get_ext(value::SteadyStateNODEObs) = value.ext
"""Get [`SteadyStateNODEObs`](@ref) `internal`."""
PSY.get_internal(value::SteadyStateNODEObs) = value.internal

"""Set [`SteadyStateNODEObs`](@ref) `base_power`."""
PSY.set_base_power!(value::SteadyStateNODEObs, val) = value.base_power = val
"""Set [`SteadyStateNODEObs`](@ref) `ext`."""
PSY.set_ext!(value::SteadyStateNODEObs, val) = value.ext = val
"""Set [`SteadyStateNODEObs`](@ref) `initializer_parameters`."""
set_initializer_parameters!(value::SteadyStateNODEObs, val) =
    value.initializer_parameters = val
"""Set [`SteadyStateNODEObs`](@ref) `node_parameters`."""
set_node_parameters!(value::SteadyStateNODEObs, val) = value.node_parameters = val
"""Set [`SteadyStateNODEObs`](@ref) `observer_parameters`."""
set_observer_parameters!(value::SteadyStateNODEObs, val) = value.observer_parameters = val

function get_SteadyStateNODEObs_states(dim_r::Int64)
    states = Symbol[]
    for i in 1:dim_r
        push!(states, Symbol(string("r", i)))
    end
    return states, dim_r
end

PSID.get_inner_vars_count(::SteadyStateNODEObs) = 0
function PSID._get_frequency_state(d::PSID.DynamicWrapper{SteadyStateNODEObs})
    return 0
end

function PSID.device_mass_matrix_entries!(
    mass_matrix::AbstractArray,
    dynamic_device::PSID.DynamicWrapper{SteadyStateNODEObs},
)
    global_index = PSID.get_global_index(dynamic_device)
    mass_matrix_ssnode_entries!(mass_matrix, dynamic_device, global_index)
    return
end

function mass_matrix_ssnode_entries!(
    mass_matrix,
    pvs::PSID.DynamicWrapper{SteadyStateNODEObs},
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
    dynamic_device::PSID.DynamicWrapper{SteadyStateNODEObs},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES}
    θ = dynamic_device.ext["θ0"]
    vd, vq = PSID.ri_dq(θ) * [voltage_r, voltage_i]
    v_scaled = _input_scale(dynamic_device, [vd, vq])
    refs = dynamic_device.ext["refs"]
    output_ode .= _forward_pass_node(dynamic_device, device_states, v_scaled, refs)
    id_scale = _forward_pass_observer(dynamic_device, device_states)[1]
    iq_scale = _forward_pass_observer(dynamic_device, device_states)[2]
    id, iq = _target_scale_inverse(dynamic_device, [id_scale, iq_scale])
    ir, ii = PSID.dq_ri(θ) * [id, iq]
    current_r[1] += ir
    current_i[1] += ii
    return
end

function PSID.initialize_dynamic_device!(
    dynamic_device::PSID.DynamicWrapper{SteadyStateNODEObs},
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
    Id_scale, Iq_scale = _target_scale(dynamic_device, [Id, Iq])    #Was input_scale! bug fixed
    function f!(out, x)
        input_scaled = _input_scale(dynamic_device, [Vd, Vq])
        out[1:n_states] .= _forward_pass_node(
            dynamic_device,
            x[1:n_states],
            input_scaled,
            x[(n_states + 1):(n_states + 2)],
        )
        x1, x2 = _target_scale_inverse(
            dynamic_device,
            [
                _forward_pass_observer(dynamic_device, x[1:n_states])[1],
                _forward_pass_observer(dynamic_device, x[1:n_states])[2],
            ],
        )
        out[n_states + 1] = x1 - Id
        out[n_states + 2] = x2 - Iq
    end
    x0 = _forward_pass_initializer(dynamic_device, Vq_scale, [Id_scale, Iq_scale])
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in SteadyStateNODEObs failed")
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
    dynamic_device::SteadyStateNODEObs,
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

function _allocate_weights_and_biases(dynamic_device::SteadyStateNODEObs)
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

get_W_initializer(wrapper::PSID.DynamicWrapper{SteadyStateNODEObs}) =
    wrapper.ext["initializer"]["W"]
get_b_initializer(wrapper::PSID.DynamicWrapper{SteadyStateNODEObs}) =
    wrapper.ext["initializer"]["b"]

get_W_node(wrapper::PSID.DynamicWrapper{SteadyStateNODEObs}) = wrapper.ext["node"]["W"]
get_b_node(wrapper::PSID.DynamicWrapper{SteadyStateNODEObs}) = wrapper.ext["node"]["b"]

get_W_observer(wrapper::PSID.DynamicWrapper{SteadyStateNODEObs}) =
    wrapper.ext["observer"]["W"]
get_b_observer(wrapper::PSID.DynamicWrapper{SteadyStateNODEObs}) =
    wrapper.ext["observer"]["b"]

function _input_scale(wrapper::PSID.DynamicWrapper{SteadyStateNODEObs}, x)
    xmin = get_input_min(wrapper.device)
    xmax = get_input_max(wrapper.device)
    l, u = get_input_lims(wrapper.device)
    return min_max_normalization(x, xmin, xmax, u, l)
end

function _target_scale(wrapper::PSID.DynamicWrapper{SteadyStateNODEObs}, x)
    xmin = get_target_min(wrapper.device)
    xmax = get_target_max(wrapper.device)
    l, u = get_target_lims(wrapper.device)
    return min_max_normalization(x, xmin, xmax, u, l)
end

function _target_scale_inverse(wrapper::PSID.DynamicWrapper{SteadyStateNODEObs}, x)
    xmin = get_target_min(wrapper.device)
    xmax = get_target_max(wrapper.device)
    l, u = get_target_lims(wrapper.device)
    return min_max_normalization_inverse(x, xmin, xmax, u, l)
end

function _forward_pass_initializer(
    wrapper::PSID.DynamicWrapper{SteadyStateNODEObs},
    input,
    target,
)
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

function _forward_pass_node(wrapper::PSID.DynamicWrapper{SteadyStateNODEObs}, r, ex, refs)
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

function _forward_pass_observer(wrapper::PSID.DynamicWrapper{SteadyStateNODEObs}, r)
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

#TODO - this function doesn't work yet. Need the dynamic wrapper to computer the output current.
#Issue is similar to this: https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/issues/269
#= function PSID.compute_output_current(
    res::PSID.SimulationResults,
    dynamic_device::SteadyStateNODEObs,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    dt::Union{Nothing, Float64},
)
    name = PSY.get_name(dynamic_device)
    n_states = PSY.get_n_states(dynamic_device)
    ts, _ = PSID.post_proc_state_series(res, (name, :r1), dt)
    states = []
    for ix in 1:n_states
        _, r = PSID.post_proc_state_series(res, (name, Symbol("r$(ix)")), dt)
        if ix == 1
            states = r'
        else
            hcat(states, r')
        end
    end
    I_R = _forward_pass_observer.(states)[1, :]
    I_I = _forward_pass_observer.(states)[2, :]

    return ts, I_R, I_I
end =#
