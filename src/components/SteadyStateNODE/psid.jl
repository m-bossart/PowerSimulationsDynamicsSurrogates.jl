PSID.is_valid(::SteadyStateNODE) = nothing

PSID.get_inner_vars_count(::SteadyStateNODE) = 0
function PSID._get_frequency_state(d::PSID.DynamicWrapper{SteadyStateNODE})
    return 0
end

function PSID._get_refs(x::PSID.DynamicWrapper{SteadyStateNODE})
    return (s1 = 0.0, s2 = 0.0, θ0 = 0.0)
end

function PSID._get_refs_metadata(x::PSID.DynamicWrapper{SteadyStateNODE})
    return (
        s1 = PSID.ParamsMetadata(PSID.DEVICE_SETPOINT, false, true),
        s2 = PSID.ParamsMetadata(PSID.DEVICE_SETPOINT, false, true),
        θ0 = PSID.ParamsMetadata(PSID.DEVICE_SETPOINT, false, true),
    )
end

PSID.get_params(x::SteadyStateNODE) =
    (; θ = (node = PSY.get_ext(x)["ps_node"], init = PSY.get_ext(x)["ps_init"]))
PSID.get_params_metadata(::SteadyStateNODE) = (;
    θ = (
        node = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
        init = PSID.ParamsMetadata(PSID.DEVICE_PARAM, false, true),
    ),
)

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
    device_states::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{T},
    p::AbstractArray{<:PSID.ACCEPTED_REAL_TYPES},
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    ::AbstractArray{T},
    ::AbstractArray{T},
    dynamic_device::PSID.DynamicWrapper{SteadyStateNODE},
    h,
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES}
    ext_wrapper = PSID.get_ext(dynamic_device)
    device = PSID.get_device(dynamic_device)
    ext_device = PSY.get_ext(device)

    model_node = ext_device["model_node"]
    ps_node = p[:params][:θ][:node]
    st_node = ext_device["st_node"]

    θ = p[:refs][:θ0]
    refs = [p[:refs][:s1], p[:refs][:s2]]
    vd, vq = PSID.ri_dq(θ) * [voltage_r, voltage_i]
    node_input = (device_states, [vd, vq], refs)
    y, _ = model_node(node_input, ps_node, st_node)
    output_ode .= y
    ir, ii = PSID.dq_ri(θ) * [device_states[1], device_states[2]]
    current_r[1] += ir
    current_i[1] += ii
    return
end

function PSID.initialize_dynamic_device!(
    dynamic_device::PSID.DynamicWrapper{SteadyStateNODE},
    source::PSY.Source,
    ::AbstractVector,
    p::AbstractVector,
    device_states::AbstractVector,
)
    ext_wrapper = PSID.get_ext(dynamic_device)
    device = PSID.get_device(dynamic_device)
    ext_device = PSY.get_ext(device)

    model_init = ext_device["model_init"]
    ps_init = p[:params][:θ][:init]
    st_init = ext_device["st_init"]
    model_node = ext_device["model_node"]
    ps_node = p[:params][:θ][:node]
    st_node = ext_device["st_node"]

    n_states = PSY.get_n_states(dynamic_device)

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
    function f!(out, x, ps_node)
        node_input = (x[1:n_states], [Vd, Vq], x[(n_states + 1):(n_states + 2)])
        y, _ = model_node(node_input, ps_node, st_node)
        out[1:n_states] .= y
        out[n_states + 1] = x[1] - Id
        out[n_states + 2] = x[2] - Iq
    end
    input_init = ([Vq], [Id, Iq])
    x0, _ = model_init(input_init, ps_init, st_init)
    prob = NonlinearSolve.NonlinearProblem{true}(f!, x0, ps_node)
    sol = NonlinearSolve.solve(
        prob,
        NonlinearSolve.TrustRegion();
        reltol = PSID.STRICT_NLSOLVE_F_TOLERANCE,
        abstol = PSID.STRICT_NLSOLVE_F_TOLERANCE,
    )
    if !SciMLBase.successful_retcode(sol)
        @warn("Initialization of SteadyStateNODE failed")
    else
        sol_x0 = sol.u
        #@warn "nlsolve result in SteadyStateNODE $sol_x0"
        device_states .= sol_x0[1:n_states]
        refs = sol_x0[(n_states + 1):(n_states + 2)]

        @view(p[:refs])[:s1] = refs[1]
        @view(p[:refs])[:s2] = refs[2]
        @view(p[:refs])[:θ0] = θ
        ext_wrapper["initializer_error"] = x0 .- sol_x0
    end
    return
end

function PSID.DynamicWrapper(
    static_device::PSY.Source,
    dynamic_device::SteadyStateNODE,
    bus_ix::Int,
    ix_range,
    ode_range,
    inner_var_range,
    sys_base_power,
    sys_base_freq,
)
    device_states = PSY.get_states(dynamic_device)

    return PSID.DynamicWrapper(
        dynamic_device,
        sys_base_power,
        sys_base_freq,
        static_device,
        PSID.BUS_MAP[PSY.get_bustype(PSY.get_bus(static_device))],
        1.0,
        collect(inner_var_range),
        collect(ix_range),
        collect(ode_range),
        bus_ix,
        Base.ImmutableDict(Dict(device_states .=> ix_range)...),
        Base.ImmutableDict{Int, Vector{Int}}(),
        Base.ImmutableDict{Int, Vector{Int}}(),
        Dict{String, Any}(),
    )
end

#TODO - this function doesn't work yet. Need the dynamic wrapper to computer the output current.
#Issue is similar to this: https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/issues/269
#= function PSID.compute_output_current(
    res::PSID.SimulationResults,
    dynamic_device::SteadyStateNODE,
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
