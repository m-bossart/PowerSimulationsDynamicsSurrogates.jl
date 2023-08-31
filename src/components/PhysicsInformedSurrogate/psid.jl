#= 
What are entries in ext? when are they created? 
    1 . "nn" -> weights and biases of the model in matrix form for forward pass calculation
             -> created when creating the  static wrapper
    2.  "references" -> vector of references to be passed as part of the model input
             -> created during initialization 
    3.  "past_values" -> cache of past values to be passed as part of the model input
             -> created during initialization; modified every timestep via callback
 =#

function PSID.StaticWrapper(device::PhysicsInformedSurrogate, bus_ix::Int)
    bus = PSY.get_bus(device)
    ext = Dict{String, Any}()
    model_architecture = get_model_architecture(device)
    model_parameters = get_model_parameters(device)
    _allocate_model_parameters!(ext, model_architecture, model_parameters)
    return PSID.StaticWrapper{PhysicsInformedSurrogate, PSID.BUS_MAP[PSY.get_bustype(bus)]}(
        device,
        Base.Ref(1.0),
        Base.Ref(PSY.get_internal_voltage(device)),
        Base.Ref(PSY.get_internal_angle(device)),
        Base.Ref(PSY.get_active_power(device)),
        Base.Ref(PSY.get_reactive_power(device)),
        bus_ix,
        ext,
    )
end

function PSID.initialize_static_device!(
    device::PSID.StaticWrapper{PhysicsInformedSurrogate, T},
) where {T <: PSID.BusCategory}
    ext = PSID.get_ext(device)

    P0 = PSY.get_active_power(device)
    Q0 = PSY.get_reactive_power(device)
    Vm = PSY.get_magnitude(PSY.get_bus(device))
    θ = PSY.get_angle(PSY.get_bus(device))
    S0 = P0 + Q0 * 1im

    x0_dict, refs = init_underlying_device(
        get_underlying_dynamic_model(PSID.get_device(device)),
        P0,
        Q0,
        Vm,
        θ,
    )
    ext["references"] = [refs["P_ref"], refs["Q_ref"], refs["V_ref"]]

    VR0 = Vm * cos(θ)
    VI0 = Vm * sin(θ)

    past_values = Dict(
        "t" => repeat([0.0], device.device.n_past_timesteps),
        "vr" => repeat([VR0], device.device.n_past_timesteps),
        "vi" => repeat([VI0], device.device.n_past_timesteps),
        "ω_sys" => repeat([1.0], device.device.n_past_timesteps),
    )
    for (key, value) in x0_dict
        past_values[string(key)] = repeat([value], device.device.n_past_timesteps)
    end
    ext["past_values"] = past_values
    PSID.set_V_ref(device, Vm)
    PSID.set_θ_ref(device, θ)
    return
end

function PSID.device!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    global_vars::AbstractArray{T},
    ::AbstractArray{T},
    device::PSID.StaticWrapper{PhysicsInformedSurrogate, U},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES, U <: PSID.BusCategory}
    _ = mdl_solution_prediction_surrogate!(     #retrun x (states) for use in callback for caching values 
        voltage_r,
        voltage_i,
        current_r,
        current_i,
        global_vars,
        device,
        t,
    )
    return
end

function mdl_solution_prediction_surrogate!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    global_vars, #::AbstractArray{T},
    static_device::PSID.StaticWrapper{PhysicsInformedSurrogate},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES}
    ext = PSID.get_ext(static_device)
    ω_sys = global_vars[PSID.GLOBAL_VAR_SYS_FREQ_INDEX]
    underlying_dynamic_model = get_underlying_dynamic_model(static_device.device)
    n_states = PSY.get_n_states(underlying_dynamic_model)
    θ_ref_frame = ext["past_values"][_get_ref_frame_state(
        static_device.device.underlying_dynamic_model,
    )][1]

    static_input, dynamic_input =
        _build_inputs(static_device.device, ext, voltage_r, voltage_i, t, θ_ref_frame)
    @assert size(dynamic_input) == (static_device.device.n_past_timesteps, 3 + n_states)  #vd, vq, Δt, states...                
    @assert size(static_input) == (1, 3)                                       #Pref, Qref, Vref                    
    data_scaler = get_data_scaler(PSID.get_device(static_device))
    model_architecture = get_model_architecture(PSID.get_device(static_device))
    weights_dict = static_device.ext["nn"]

    static_input_scaled = _input_scale(data_scaler, static_input, 1:3)
    dynamic_input_scaled = _input_scale(data_scaler, dynamic_input, 4:(n_states + 6))

    model_output_scaled = _call_model(
        weights_dict,
        model_architecture,
        static_input_scaled,
        dynamic_input_scaled,
    )
    model_output = _target_scale_inverse(data_scaler, model_output_scaled)                    #TODO - the output is the states... need to handle the state order properly

    ir, ii = _calculate_current_from_states(underlying_dynamic_model, model_output)   #TODO- implement! 

    current_r[1] += ir
    current_i[1] += ii
    return model_output
end

function _calculate_current_from_states(underlying_dynamic_model, x0)
    ir = 0.0
    ii = 0.0
    return ir, ii
end

function _build_inputs(
    static_device::PhysicsInformedSurrogate,
    ext::Dict{String, Any},
    voltage_r,
    voltage_i,
    t,
    ref_frame_angle,
)
    n_states = PSY.get_n_states(get_underlying_dynamic_model(static_device))
    states = PSY.get_states(get_underlying_dynamic_model(static_device))
    dynamic_input = Array{Any}(undef, static_device.n_past_timesteps, n_states + 3)  #3 -> vr, vi, t
    static_input = Array{Any}(undef, 1, 3) # 3-> P_ref, Q_ref, V_ref 

    for i in 1:(static_device.n_past_timesteps)
        if i == 1
            dynamic_input[i, :] = vcat(
                PSID.ri_dq(ref_frame_angle) * [voltage_r; voltage_i],
                zeros(n_states),              # TODO - should these be set to zero before or after scaling?
                [t - ext["past_values"]["t"][1]],
            )
        else
            dynamic_input[i, :] = vcat(
                PSID.ri_dq(ref_frame_angle) *
                [ext["past_values"]["vr"][i]; ext["past_values"]["vi"][i]],
                [ext["past_values"][string(x)][i] for x in states], # zeros(n_states), #TODO -get past values of states, need strategy to consider order of states
                [ext["past_values"]["t"][i] - ext["past_values"]["t"][i - 1]],
            )
        end
    end
    static_input[1, :] = ext["references"]
    return static_input, dynamic_input
end

mutable struct PhysicsInformedSurrogateCacheValues <: PSID.Perturbation
    device::PhysicsInformedSurrogate
end

function PSID._add_callback!(
    tstops::Vector{Float64},
    callback_vector::Vector{SciMLBase.DiscreteCallback},
    ix::Int,
    pert::PhysicsInformedSurrogateCacheValues,
    sim::PSID.Simulation,
    inputs::PSID.SimulationInputs,
)
    wrapped_device_ix = _find_device_index(inputs, pert.device)
    global_vars = PSID.get_global_vars_update_pointers(sim.inputs)  #not sure if this is correct
    f =
        (u, t, integrator) -> begin
            static_device = PSID.get_static_injectors(integrator.p)[wrapped_device_ix]
            n_buses = PSID.get_n_buses(sim.sys)
            bus_ix = PSID.get_bus_ix(static_device)
            voltage_r = integrator.u[bus_ix]
            voltage_i = integrator.u[bus_ix + n_buses]

            current_r = [0.0]
            current_i = [0.0]
            states = mdl_solution_prediction_surrogate!(
                voltage_r,
                voltage_i,
                current_r,
                current_i,
                global_vars,
                static_device,
                t,
            )
            past_values_dict = PSID.get_ext(static_device)["past_values"]
            past_values_dict["t"] = vcat(t, past_values_dict["t"][1:(end - 1)])
            past_values_dict["ω_sys"] = vcat(
                global_vars[PSID.GLOBAL_VAR_SYS_FREQ_INDEX],
                past_values_dict["ω_sys"][1:(end - 1)],
            )
            past_values_dict["vr"] =
                vcat(voltage_r, past_values_dict["vr"][1:(end - 1)])
            past_values_dict["vi"] =
                vcat(voltage_i, past_values_dict["vi"][1:(end - 1)])
            for (ix, k) in enumerate(
                PSY.get_states(get_underlying_dynamic_model(static_device.device)),
            )
                past_values_dict[string(k)] =
                    vcat(states[ix], past_values_dict[string(k)][1:(end - 1)])    #State ordering?? 
            end
            PSID.get_ext(static_device)["past_values"] = past_values_dict
        end
    callback_vector[ix] = DiffEqCallbacks.FunctionCallingCallback(f)
end

#NOTE:  will be hard to post-process compute the output current of the surrogate. For now, rely on measuring current on a connecting line. 
