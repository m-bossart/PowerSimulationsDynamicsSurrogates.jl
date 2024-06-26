function hardtanh(x)
    max(-1, min(1, x))
end

sigmoid(x) = 1 / (1 + exp(-x))

function relu(x)
    max(0, x)
end

function activation_map(activation)
    d = Dict("tanh" => tanh, "hardtanh" => hardtanh, "relu" => relu, "identity" => (x) -> x)
    return d[activation]
end

function _get_ref_frame_state(device::D) where {D <: PSY.DynamicInjection}
    @error "_get_ref_frame_state not implemented for type $(typeof(device))"
end

function _get_ref_frame_state(
    device::PSY.DynamicInverter{
        PSY.AverageConverter,
        PSY.OuterControl{PSY.ActivePowerDroop, PSY.ReactivePowerDroop},
        PSY.VoltageModeControl,
        PSY.FixedDCSource,
        PSY.KauraPLL,
        PSY.LCLFilter,
    },
)
    return "θ_pll"
end

#Gets the initial condition for the underlying device by building a toy 1-bus system with source + device
function init_underlying_device(
    device::D,
    P0,
    Q0,
    Vm0,
    θ0,
) where {D <: PSY.DynamicInjection}
    sys = PSY.System(100.0)
    b = PSY.ACBus(
        1,
        "Bus1",
        PSY.ACBusTypes.REF,
        θ0,
        Vm0,
        (min = 0.0, max = 2.0),
        1.0,
        nothing,
    )
    PSY.add_component!(sys, b; skip_validation = true)
    gen = PSY.ThermalStandard(;
        name = PSY.get_name(device),
        available = true,
        status = true,
        bus = b,
        active_power = P0,
        reactive_power = Q0,
        rating = 1.0,
        active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = nothing,
        ramp_limits = nothing,
        operation_cost = PSY.ThreePartCost(nothing),
        base_power = PSY.get_base_power(device),
        time_limits = nothing,
        must_run = false,
        prime_mover_type = PSY.PrimeMovers.OT,
        fuel = PSY.ThermalFuels.OTHER,
        services = PSY.Device[],
    )
    source = PSY.Source(;
        name = "gen_perturb",
        available = true,
        bus = b,
        active_power = -P0,
        reactive_power = -Q0,
        R_th = 0.001,
        X_th = 0.01,
        internal_voltage = Vm0,
        internal_angle = θ0,
    )
    PSY.add_component!(sys, gen)
    PSY.add_component!(sys, source)
    PSY.add_component!(sys, device, gen)
    rng_state = copy(Random.default_rng())
    sim = PSID.Simulation!(
        PSID.MassMatrixModel,
        sys,
        pwd(),
        (0.0, 1.0);
        disable_timer_outputs = true,
    )
    copy!(Random.default_rng(), rng_state)
    refs = PSID.get_setpoints(sim)[PSY.get_name(device)]
    x0 = PSID.read_initial_conditions(sim)[PSY.get_name(device)]
    return x0, refs
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

function min_max_normalization(x, xmin, xmax, u, l)  #https://www.baeldung.com/cs/normalizing-inputs-artificial-neural-network
    x_prime = (x .- xmin) ./ (xmax .- xmin) .* (u .- l) .+ l
    return x_prime
end

function min_max_normalization_inverse(x_prime, xmin, xmax, u, l)
    x = (x_prime .- l) .* (xmax .- xmin) ./ (u .- l) .+ xmin
    return x
end

function to_json_with_surrogates(sys, full_path)
    dir = dirname(full_path)
    mkpath(joinpath(dir, "surrogate_models"))
    for g in PSY.get_components(DataDrivenSurrogate, sys)
        if typeof(g) == SteadyStateNODE
            model_node = PSY.get_ext(g)["model_node"]
            p_node = PSY.get_ext(g)["ps_node"]
            st_node = PSY.get_ext(g)["st_node"]
            model_init = PSY.get_ext(g)["model_init"]
            p_init = PSY.get_ext(g)["ps_init"]
            st_init = PSY.get_ext(g)["st_init"]

            mkpath(joinpath(dir, "surrogate_models", PSY.get_name(g)))
            BSON.@save joinpath(dir, "surrogate_models", PSY.get_name(g), "node") model_node p_node st_node
            BSON.@save joinpath(dir, "surrogate_models", PSY.get_name(g), "init") model_init p_init st_init
            PSY.set_ext!(
                g,
                Dict{String, Any}(
                    "node_path" =>
                        joinpath(dir, "surrogate_models", PSY.get_name(g), "node"),
                    "init_path" =>
                        joinpath(dir, "surrogate_models", PSY.get_name(g), "init"),
                    "model_node" => nothing,
                    "ps_node" => nothing,
                    "st_node" => nothing,
                    "model_init" => nothing,
                    "ps_init" => nothing,
                    "st_init" => nothing,
                ),
            )
        else
            model = PSY.get_ext(g)["model"]
            p = PSY.get_ext(g)["ps"]
            st = PSY.get_ext(g)["st"]
            dir = dirname(full_path)

            mkpath(joinpath(dir, "surrogate_models", PSY.get_name(g)))
            BSON.@save joinpath(dir, "surrogate_models", PSY.get_name(g), "nn") model p st
            PSY.set_ext!(
                g,
                Dict{String, Any}(
                    "model_path" =>
                        joinpath(dir, "surrogate_models", PSY.get_name(g), "nn"),
                    "model" => nothing,
                    "ps" => nothing,
                    "st" => nothing,
                ),
            )
        end
    end
    PSY.to_json(sys, full_path; force = true)
end

function deserialize_with_surrogates(full_path)
    sys = PSY.System(full_path)
    for g in PSY.get_components(DataDrivenSurrogate, sys)
        if typeof(g) == SteadyStateNODE
            BSON.@load PSY.get_ext(g)["node_path"] model_node p_node st_node
            BSON.@load PSY.get_ext(g)["init_path"] model_init p_init st_init
            PSY.set_ext!(
                g,
                Dict{String, Any}(
                    "node_path" => nothing,
                    "init_path" => nothing,
                    "model_node" => model_node,
                    "ps_node" => p_node,
                    "st_node" => st_node,
                    "model_init" => model_init,
                    "ps_init" => p_init,
                    "st_init" => st_init,
                ),
            )
        else
            BSON.@load PSY.get_ext(g)["model_path"] model p st
            PSY.set_ext!(
                g,
                Dict{String, Any}(
                    "model_path" => nothing,
                    "model" => model,
                    "ps" => p,
                    "st" => st,
                ),
            )
        end
    end
    return sys
end
