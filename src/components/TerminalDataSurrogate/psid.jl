#= 
What are entries in ext? when are they created? 
    1 . "nn" -> weights and biases of the model in matrix form for forward pass calculation
             -> created when creating the  static wrapper
    2.  "references" -> vector of references to be passed as part of the model input
             -> created during initialization 
    3.  "past_values" -> cache of past values to be passed as part of the model input
             -> created during initialization; modified every timestep via callback
 =#

function PSID.StaticWrapper(device::TerminalDataSurrogate, bus_ix::Int)
    bus = PSY.get_bus(device)
    ext = Dict{String, Any}()
    model_architecture = get_model_architecture(device)
    model_parameters = get_model_parameters(device)
    _allocate_model_parameters!(ext, model_architecture, model_parameters)
    return PSID.StaticWrapper{TerminalDataSurrogate, PSID.BUS_MAP[PSY.get_bustype(bus)]}(
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
    device::PSID.StaticWrapper{TerminalDataSurrogate, T},
) where {T <: PSID.BusCategory}
    ext = PSID.get_ext(device)

    P0 = PSY.get_active_power(device)
    Q0 = PSY.get_reactive_power(device)
    Vm = PSY.get_magnitude(PSY.get_bus(device))
    θ = PSY.get_angle(PSY.get_bus(device))
    S0 = P0 + Q0 * 1im

    _, refs = init_underlying_device(
        get_underlying_dynamic_model(PSID.get_device(device)),
        P0,
        Q0,
        Vm,
        θ,
    )
    ext["references"] = [refs["P_ref"], refs["Q_ref"], refs["V_ref"]]

    VR0 = Vm * cos(θ)
    VI0 = Vm * sin(θ)
    set_θ_ref_frame!(device.device, atan(VI0, VR0))

    V = VR0 + VI0 * 1im
    I = conj(S0 / V)
    IR0 = real(I)
    II0 = imag(I)
    past_values = Dict(
        "t" => repeat([0.0], device.device.n_past_timesteps),
        "vr" => repeat([VR0], device.device.n_past_timesteps),
        "vi" => repeat([VI0], device.device.n_past_timesteps),
        "ir" => repeat([IR0], device.device.n_past_timesteps),
        "ii" => repeat([II0], device.device.n_past_timesteps),
        "ω_sys" => repeat([1.0], device.device.n_past_timesteps),
    )
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
    device::PSID.StaticWrapper{TerminalDataSurrogate, U},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES, U <: PSID.BusCategory}
    mdl_solution_prediction_surrogate!(
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
    static_device::PSID.StaticWrapper{TerminalDataSurrogate},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES}
    ext = PSID.get_ext(static_device)
    ω_sys = global_vars[PSID.GLOBAL_VAR_SYS_FREQ_INDEX]
    θ_ref_frame = get_θ_ref_frame(static_device.device)

    static_input, dynamic_input =
        _build_inputs(static_device.device, ext, voltage_r, voltage_i, t, θ_ref_frame)
    @assert size(dynamic_input) == (static_device.device.n_past_timesteps, 5)  #vd, vq, id, iq, Δt
    @assert size(static_input) == (1, 3)                                       #Pref, Qref, Vref
    data_scaler = get_data_scaler(PSID.get_device(static_device))
    model_architecture = get_model_architecture(PSID.get_device(static_device))
    parameter_arrays = static_device.ext["nn"]

    static_input_scaled = _input_scale(data_scaler, static_input, 1:3)  #dispatch on DataScaler - indicate indices to scale (e.g 1:3, 4:8) for static and dynamic 
    dynamic_input_scaled = _input_scale(data_scaler, dynamic_input, 4:8)

    model_output_scaled = _call_model(
        parameter_arrays,
        model_architecture,
        static_input_scaled,
        dynamic_input_scaled,
    )   #dispatch on DataDrivenModelArchitecture

    id, iq = _target_scale_inverse(data_scaler, model_output_scaled) #dispatch on DataScaler
    ir, ii = PSID.dq_ri(θ_ref_frame) * [id; iq]
    current_r[1] += ir
    current_i[1] += ii
    return
end

function _build_inputs(
    static_device::TerminalDataSurrogate,
    ext::Dict{String, Any},
    voltage_r,
    voltage_i,
    t,
    ref_frame_angle,
)
    dynamic_input = Array{Any}(undef, static_device.n_past_timesteps, 5)
    static_input = Array{Any}(undef, 1, 3)
    for i in 1:(static_device.n_past_timesteps)
        if i == 1
            dynamic_input[i, :] = vcat(
                [t - ext["past_values"]["t"][1]],
                PSID.ri_dq(ref_frame_angle) * [voltage_r; voltage_i] .- PSID.ri_dq(ref_frame_angle) * [ext["past_values"]["vr"][1]; ext["past_values"]["vi"][1]],
                PSID.ri_dq(ref_frame_angle) * [ext["past_values"]["ir"][1]; ext["past_values"]["ii"][1]] .- PSID.ri_dq(ref_frame_angle) * [ext["past_values"]["ir"][2]; ext["past_values"]["ii"][2]],
            )
        else    #voltage and time are offset from current by one timestep
            dynamic_input[i, :] = vcat(
                [ext["past_values"]["t"][i-1] - ext["past_values"]["t"][i]],
                PSID.ri_dq(ref_frame_angle) * [ext["past_values"]["vr"][i-1]; ext["past_values"]["vi"][i-1]] .- PSID.ri_dq(ref_frame_angle) * [ext["past_values"]["vr"][i]; ext["past_values"]["vi"][i]],
                PSID.ri_dq(ref_frame_angle) * [ext["past_values"]["ir"][1]; ext["past_values"]["ii"][1]] .- PSID.ri_dq(ref_frame_angle) * [ext["past_values"]["ir"][i+1]; ext["past_values"]["ii"][i+1]],
            )
        end
    end
    static_input[1, :] = ext["references"]
    return static_input, dynamic_input
end

mutable struct TerminalDataSurrogateCacheValues <: PSID.Perturbation
    device::TerminalDataSurrogate
end

function PSID._add_callback!(
    tstops::Vector{Float64},
    callback_vector::Vector{SciMLBase.DiscreteCallback},
    ix::Int,
    pert::TerminalDataSurrogateCacheValues,
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
            mdl_solution_prediction_surrogate!(
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
            past_values_dict["ir"] =
                vcat(current_r[1], past_values_dict["ir"][1:(end - 1)])
            past_values_dict["ii"] =
                vcat(current_i[1], past_values_dict["ii"][1:(end - 1)])
            PSID.get_ext(static_device)["past_values"] = past_values_dict
        end
    callback_vector[ix] = DiffEqCallbacks.FunctionCallingCallback(f)
end

#NOTE:  will be hard to post-process compute the output current of the surrogate. For now, rely on measuring current on a connecting line. 
