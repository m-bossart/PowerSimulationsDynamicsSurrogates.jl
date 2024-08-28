#This is a workaround to enable reverse mode training based on the following restrictions:
# --> ReverseDiffAdjoint is not compatible with Enzyme
# --> Zygote is not compatible with mutation (PSID initialization routines mutate x0, p)
#Uses Zygote.Buffer to initialize only the ML models in a system without hitting Zygote mutation restriction.

function _make_buffer(x)
    x_buff = Zygote.Buffer(x)
    for i in eachindex(x)
        x_buff[i] = x[i]
    end
    return x_buff
end

function PSID._non_mutating_initialization_of_ml_surrogates(
    x0,
    p,
    sim_inputs::PSID.SimulationInputs,
)
    x0_buffer = _make_buffer(x0)
    p_buffer = _make_buffer(p)

    for dynamic_device in PSID.get_dynamic_injectors(sim_inputs)
        if (typeof(PSID.get_device(dynamic_device)) == TerminalDataSurrogate) ||
           (typeof(PSID.get_device(dynamic_device)) == SteadyStateNODE)
            static = PSID.get_static_device(dynamic_device)
            p_device = p_buffer[PSID._get_wrapper_name(dynamic_device)]
            x0_device = x0_buffer[PSID.get_ix_range(dynamic_device)]
            p_buffer_device = _make_buffer(p_device)
            x0_buffer_device = _make_buffer(x0_device)
            PSID.initialize_dynamic_device!(
                dynamic_device,
                static,
                zeros(1),
                p_buffer_device,
                x0_buffer_device,
            )  #does not use inner vars 
            p_buffer[PSID._get_wrapper_name(dynamic_device)] = copy(p_buffer_device)
            x0_buffer[PSID.get_ix_range(dynamic_device)] = copy(x0_buffer_device)
        end
    end
    x0_modified = copy(x0_buffer)
    p_modified = copy(p_buffer)
    return x0_modified, p_modified
end

function PSID.initialize_dynamic_device!(
    dynamic_device::PSID.DynamicWrapper{SteadyStateNODE},
    source::PSY.Source,
    ::AbstractVector,
    p::Zygote.Buffer,
    device_states::Zygote.Buffer,
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

    Vd, Vq = [
        sin(θ) -cos(θ)
        cos(θ) sin(θ)
    ] * [V_R, V_I] #Vd should be zero by definition 
    Id, Iq = [
        sin(θ) -cos(θ)
        cos(θ) sin(θ)
    ] * [I_R, I_I]
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
        for i in 1:n_states
            device_states[i] = sol_x0[i]
        end
        refs = sol_x0[(n_states + 1):(n_states + 2)]
        p_refs = _make_buffer(p[:refs])
        p_refs[1] = refs[1] #s1 
        p_refs[2] = refs[2] #s2
        p_refs[3] = θ       #θ0
        p[:refs] = copy(p_refs)
        ext_wrapper["initializer_error"] = x0 .- sol_x0
    end
    return
end

function PSID.initialize_dynamic_device!(
    dynamic_device::PSID.DynamicWrapper{TerminalDataSurrogate},
    source::PSY.Source,
    ::AbstractVector,
    p::Zygote.Buffer,
    device_states::Zygote.Buffer,
)
    ext_wrapper = PSID.get_ext(dynamic_device)
    device = PSID.get_device(dynamic_device)
    ext_device = PSY.get_ext(device)

    P0 = PSY.get_active_power(source)
    Q0 = PSY.get_reactive_power(source)
    Vm = PSY.get_magnitude(PSY.get_bus(source))
    θ = PSY.get_angle(PSY.get_bus(source))
    S0 = P0 + Q0 * 1im

    VR0 = Vm * cos(θ)
    VI0 = Vm * sin(θ)
    V = VR0 + VI0 * 1im
    I = conj(S0 / V)
    IR0 = real(I)
    II0 = imag(I)

    ext_wrapper["v0"] = [VR0, VI0]
    ext_wrapper["i0"] = [IR0, II0]
    ext_wrapper["voltage_violation_warning"] = false

    model = ext_device["model"]
    ps = p[:params][:θ]
    st = ext_device["st"]
    window_size = get_window_size(device)
    #@error typeof(device_states[1])
    v_ss = Zygote.Buffer(Array{Float64}(undef, (2, window_size, 1)))
    for i in eachindex(v_ss)
        if isodd(i)
            v_ss[i] = VR0
        else
            v_ss[i] = VI0
        end
    end
    i_ss = Zygote.Buffer(Array{Float64}(undef, (2, window_size, 1)))
    for i in eachindex(i_ss)
        if isodd(i)
            i_ss[i] = IR0
        else
            i_ss[i] = II0
        end
    end
    v0 = Zygote.Buffer(Array{Float64}(undef, (2, 1)))
    v0[1] = VR0
    v0[2] = VI0
    i0 = Zygote.Buffer(Array{Float64}(undef, (2, 1)))
    i0[1] = IR0
    i0[2] = II0
    #@assert size(copy(v_ss))[1] == 2
    ss_input = (copy(v0), copy(i0), copy(v_ss), copy(i_ss))
    y, st = model(ss_input, ps, st)
    p_refs = _make_buffer(p[:refs])
    if get_steadystate_offset_correction(device)
        device_states[1] = IR0
        device_states[2] = II0
        device_states[3] = VR0
        device_states[4] = VI0
        p_refs[1] = y[1] - IR0
        p_refs[2] = y[2] - II0
    else
        p_refs[1] = 0.0
        p_refs[2] = 0.0
        device_states[1] = y[1]
        device_states[2] = y[2]
        device_states[3] = VR0
        device_states[4] = VI0
    end
    p[:refs] = copy(p_refs)
    return
end
