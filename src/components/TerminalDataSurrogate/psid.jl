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
    _allocate_weights_and_biases!(ext, model_architecture)
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

function _allocate_weights_and_biases!(
    ext::Dict{String, Any},
    model_architecture::FullyConnected,
)
    Ws = []
    bs = []
    layers = get_nn_structure(model_architecture)
    p = get_nn_parameters(model_architecture)
    p_index_start = 0
    _ = _push_layer_weights_and_biases!(Ws, bs, layers, p, p_index_start)
    ext["nn"] = Dict{String, Vector{AbstractArray}}("W" => Ws, "b" => bs)
    return
end

get_W_nn(wrapper::PSID.StaticWrapper{TerminalDataSurrogate}) = wrapper.ext["nn"]["W"]
get_b_nn(wrapper::PSID.StaticWrapper{TerminalDataSurrogate}) = wrapper.ext["nn"]["b"]

function PSID.initialize_static_device!(
    device::PSID.StaticWrapper{TerminalDataSurrogate, T},
) where {T <: PSID.BusCategory}
    ext = PSID.get_ext(device)
    #Get PowerFlow Data
    P0 = PSY.get_active_power(device)
    Q0 = PSY.get_reactive_power(device)
    Vm = PSY.get_magnitude(PSY.get_bus(device))
    θ = PSY.get_angle(PSY.get_bus(device))
    S0 = P0 + Q0 * 1im

    #TODO - we want to actually call PSID initialization methods which initialize the device.
    ext["references"] = []

    VR0 = Vm * cos(θ)
    VI0 = Vm * sin(θ)
    set_θ_ref_frame!(device.device, atan(VI0, VR0))

    V = VR0 + VI0 * 1im
    I = conj(S0 / V)
    IR0 = real(I)
    II0 = imag(I)
    ext["past_values"] = repeat([0.0, VR0, VI0, IR0, II0], device.device.n_past_timesteps)
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
    vd, vq = PSID.ri_dq(θ_ref_frame) * [voltage_r; voltage_i]
    #TODO - need a function that transforms the past values (go from vr,vi to vd,vq, calculate Δt from t, etc.)
    model_input = vcat(ext["references"], ext["past_values"], [vd, vq, ω_sys])

    model_input_scaled = _input_scale(static_device, model_input)
    model_output_scaled = _forward_pass_nn(static_device, model_input_scaled)
    id, iq = _target_scale_inverse(static_device, model_output_scaled)
    ir, ii = PSID.dq_ri(θ_ref_frame) * [id; iq]
    current_r[1] += ir
    current_i[1] += ii
    return
end

function _forward_pass_nn(wrapper::PSID.StaticWrapper{TerminalDataSurrogate}, x)
    W_index = 1
    b_index = 1
    W = get_W_nn(wrapper)
    b = get_b_nn(wrapper)
    for layer in get_nn_structure(wrapper.device.model_architecture)
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

function _calculate_input_features(
    static_device::PSID.StaticWrapper{TerminalDataSurrogate},
    past_values_dict::Dict{String, Vector{Float64}},
    voltage_r,
    voltage_i,
    t,
)
    if static_device.device.nn_features == :direct
        return vcat(
            past_values_dict["Vr"],
            past_values_dict["Vi"],
            past_values_dict["Ir"],
            past_values_dict["Ii"],
            past_values_dict["t"],
            voltage_r,
            voltage_i,
            t,
        )
    end
end

function _calculate_current_from_output(
    static_device::PSID.StaticWrapper{TerminalDataSurrogate},
    y,
    voltage_r,
    voltage_i,
    t,
)
    if static_device.device.nn_features == :direct
        return y
    end
end
function _input_scale(wrapper::PSID.StaticWrapper{TerminalDataSurrogate}, x)
    xmin = get_input_min(wrapper.device.model_architecture)
    xmax = get_input_max(wrapper.device.model_architecture)
    l, u = get_input_lims(wrapper.device.model_architecture)
    return min_max_normalization(x, xmin, xmax, u, l)
end
function _target_scale_inverse(wrapper::PSID.StaticWrapper{TerminalDataSurrogate}, x)
    xmin = get_target_min(wrapper.device.model_architecture)
    xmax = get_target_max(wrapper.device.model_architecture)
    l, u = get_target_lims(wrapper.device.model_architecture)
    return min_max_normalization_inverse(x, xmin, xmax, u, l)
end

mutable struct SolutionSurrogateCacheValues <: PSID.Perturbation
    device::TerminalDataSurrogate
end

function _is_same_device(
    device1::PSID.StaticWrapper{T},
    device2::U,
) where {T <: PSY.StaticInjection, U <: PSY.StaticInjection}
    if typeof(device1.device) != typeof(device2)                #Had to change this line from PSID because my static injector is parameterized based on other types. 
        return false
    end
    if PSY.get_name(device1) == PSY.get_name(device2)
        return true
    elseif PSY.get_name(device1) != PSY.get_name(device2)
        return false
    else
        error("comparison failed for $device1 and $device2")
    end
end

function _find_device_index(inputs::PSID.SimulationInputs, device::PSY.StaticInjection)
    wrapped_devices = PSID.get_static_injectors(inputs)
    wrapped_device_ixs = findall(x -> _is_same_device(x, device), wrapped_devices)
    if isempty(wrapped_device_ixs)
        error(
            "Device $(typeof(device))-$(PSY.get_name(device)) not found in the simulation inputs",
        )
    end
    return wrapped_device_ixs[1]
end

function PSID._add_callback!(
    tstops::Vector{Float64},
    callback_vector::Vector{SciMLBase.DiscreteCallback},
    ix::Int,
    pert::SolutionSurrogateCacheValues,
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

            current_r = [1.0]
            current_i = [1.0]
            mdl_solution_prediction_surrogate!(
                voltage_r,
                voltage_i,
                current_r,
                current_i,
                global_vars,
                static_device,
                t,
            )

            new_values = [
                voltage_r,
                voltage_i,
                current_r[1],
                current_i[1],
                global_vars[PSID.GLOBAL_VAR_SYS_FREQ_INDEX],
            ]
            past_values_vector = PSID.get_ext(static_device)["past_values"]
            past_values_vector = circshift(past_values_vector, 5)
            past_values_vector[1:5] = new_values
            PSID.get_ext(static_device)["past_values"] = past_values_vector
        end
    callback_vector[ix] = DiffEqCallbacks.FunctionCallingCallback(f)
end

#TODO - will be hard to post-process compute the output current coming out of one of these surrogates.
#Might need to rely on measuring the current someplace else or calculating based on multiple measured currents. 
