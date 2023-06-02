function PSID.StaticWrapper(device::T, bus_ix::Int) where {T <: LearnedSolutionSurrogate}
    bus = PSY.get_bus(device)
    ext = Dict{String, Any}()
    _allocate_weights_and_biases!(ext, device)  
    return PSID.StaticWrapper{T, PSID.BUS_MAP[PSY.get_bustype(bus)]}(
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

function _allocate_weights_and_biases!(ext::Dict{String, Any}, static_device::SolutionPredictionSurrogate)
    Ws = []
    bs = []
    layers = get_nn_structure(static_device)
    p = get_nn_parameters(static_device)
    p_index_start = 0
    _ = _push_layer_weights_and_biases!(Ws, bs, layers, p, p_index_start)
    ext["nn"] = Dict{String, Vector{AbstractArray}}("W" => Ws, "b" => bs)
    return 
end

get_W_nn(wrapper::PSID.StaticWrapper{SolutionPredictionSurrogate}) =
    wrapper.ext["nn"]["W"]
get_b_nn(wrapper::PSID.StaticWrapper{SolutionPredictionSurrogate}) =
    wrapper.ext["nn"]["b"]

function PSID.initialize_static_device!(
    device::PSID.StaticWrapper{SolutionPredictionSurrogate, T},
) where {T <: PSID.BusCategory}
    l = get_length_cache(device.device)
    #Get PowerFlow Data
    P0 = PSY.get_active_power(device)
    Q0 = PSY.get_reactive_power(device)
    Vm = PSY.get_magnitude(PSY.get_bus(device))
    θ = PSY.get_angle(PSY.get_bus(device))
    S0 = P0 + Q0 * 1im

    VR0 = Vm * cos(θ)
    VI0 = Vm * sin(θ)
    V = VR0 + VI0 * 1im
    I = conj(S0 / V)
    IR0 = real(I)
    II0 = imag(I)
    ext = PSID.get_ext(device)
    ext["past_values"] =  Dict("t" => zeros(l), "Vr" => ones(l) .* VR0,   "Vi" => ones(l) .* VI0, "Ir" => ones(l) .* IR0,  "Ii" => ones(l) .* II0)
    PSID.set_V_ref(device, Vm)
    PSID.set_θ_ref(device, θ)
    return
end

function PSID.device!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    ::AbstractArray{T},
    ::AbstractArray{T},
    device::PSID.StaticWrapper{SolutionPredictionSurrogate, U},
    t,
) where {T <: PSID.ACCEPTED_REAL_TYPES, U <: PSID.BusCategory}
    mdl_solution_prediction_surrogate!(voltage_r, voltage_i, current_r, current_i, device, t)  #Device (the static wrapper) has past values
    return
end

function mdl_solution_prediction_surrogate!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    static_device::PSID.StaticWrapper{SolutionPredictionSurrogate},
    t
) where {T <: PSID.ACCEPTED_REAL_TYPES}
    past_values_dict = PSID.get_ext(static_device)["past_values"]
    x_prime = _calculate_input_features(static_device, past_values_dict, voltage_r, voltage_i, t)     
    x = _input_scale(static_device, x_prime)
    y = _forward_pass_nn(static_device, x) 
    y_prime = _target_scale_inverse(static_device, y)
    ir, ii = _calculate_current_from_output(static_device, y_prime, voltage_r, voltage_i, t)
    current_r[1] += ir
    current_i[1] += ii
    return
end

function _forward_pass_nn(
    wrapper::PSID.StaticWrapper{SolutionPredictionSurrogate},
    x,
)
    W_index = 1
    b_index = 1
    W = get_W_nn(wrapper)
    b = get_b_nn(wrapper)
    for layer in get_nn_structure(wrapper.device)
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


function _calculate_input_features(static_device::PSID.StaticWrapper{SolutionPredictionSurrogate}, past_values_dict::Dict{String, Vector{Float64}}, voltage_r, voltage_i, t)
    if static_device.device.nn_features == :direct 
        return vcat(past_values_dict["Vr"], past_values_dict["Vi"], past_values_dict["Ir"], past_values_dict["Ii"], past_values_dict["t"], voltage_r, voltage_i, t)
    end 
end 

function _calculate_current_from_output(static_device::PSID.StaticWrapper{SolutionPredictionSurrogate}, y, voltage_r, voltage_i, t)
    if static_device.device.nn_features == :direct 
        return y
    end 
end 
function _input_scale(wrapper::PSID.StaticWrapper{SolutionPredictionSurrogate}, x)
    xmin = get_input_min(wrapper.device)
    xmax = get_input_max(wrapper.device)
    l, u = get_input_lims(wrapper.device)
    return min_max_normalization(x, xmin, xmax, u, l)
end
function _target_scale_inverse(wrapper::PSID.StaticWrapper{SolutionPredictionSurrogate}, x)
    xmin = get_target_min(wrapper.device)
    xmax = get_target_max(wrapper.device)
    l, u = get_target_lims(wrapper.device)
    return min_max_normalization_inverse(x, xmin, xmax, u, l)
end

mutable struct SolutionSurrogateCacheValues <: PSID.Perturbation
    device::SolutionPredictionSurrogate
end

function PSID._add_callback!(tstops::Vector{Float64}, callback_vector::Vector{SciMLBase.DiscreteCallback}, ix::Int, pert::SolutionSurrogateCacheValues, sim::PSID.Simulation, inputs::PSID.SimulationInputs) 
    wrapped_device_ix = PSID._find_device_index(inputs, pert.device)
    f = (u, t, integrator) -> begin 
    static_device = PSID.get_static_injectors(integrator.p)[wrapped_device_ix]
    n_buses = PSID.get_n_buses(sim.sys)
    bus_ix = PSID.get_bus_ix(static_device)
    voltage_r = integrator.u[bus_ix]
    voltage_i = integrator.u[bus_ix+n_buses]

    #Calculate the device current before changing "past_values"
    past_values_dict = PSID.get_ext(static_device)["past_values"]
    x_prime = _calculate_input_features(static_device, past_values_dict, voltage_r, voltage_i, t)     
    x = _input_scale(static_device, x_prime)
    y = _forward_pass_nn(static_device, x) 
    y_prime = _target_scale_inverse(static_device, y)
    ir, ii = _calculate_current_from_output(static_device, y_prime, voltage_r, voltage_i, t)
    t_old = past_values_dict["t"]
    Vr_old = past_values_dict["Vr"]
    Vi_old =  past_values_dict["Vi"]
    Ir_old = past_values_dict["Ir"]
    Ii_old =past_values_dict["Ii"]
    t_new = circshift(t_old, 1)
    t_new[1] = t 
    Vr_new = circshift(Vr_old, 1)
    Vr_new[1] = voltage_r 
    Vi_new = circshift(Vi_old, 1)
    Vi_new[1] = voltage_i 
    Ir_new = circshift(Ir_old, 1)
    Ir_new[1] = ir 
    Ii_new = circshift(Ii_old, 1)
    Ii_new[1] = ii 
    past_values_dict["t"] = t_new
    past_values_dict["Vr"] = Vr_new
    past_values_dict["Vi"] = Vi_new
    past_values_dict["Ir"] = Ir_new
    past_values_dict["Ii"] = Ii_new
    end 
    callback_vector[ix] = DiffEqCallbacks.FunctionCallingCallback(f)
    #TODO - modify this so the function actually caches the values! 
end 



#wrapped_device_ix = _find_zip_load_ix(inputs, pert.device)
#= wrapped_device_ix = _find_device_index(inputs, pert.device)
return (integrator) -> begin
    wrapped_device = get_dynamic_injectors(integrator.p)[wrapped_device_ix]
    ix_range = get_ix_range(wrapped_device)
    @debug "Changing connection status $(PSY.get_name(wrapped_device)), setting states $ix_range to 0.0"
    if integrator.du !== nothing
        @debug "setting du $ix_range to 0.0"
        integrator.du[ix_range] .= 0.0
    end
    integrator.u[ix_range] .= 0.0
    set_connection_status(wrapped_device, 0)
end =#


#TODO - will be hard to post-process compute the output current coming out of one of these surrogates.
#Might need to rely on measuring the current someplace else or calculating based on multiple measured currents. 