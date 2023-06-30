mutable struct AllStatesData <: SurrogateDataset
    type::String
    tsteps::AbstractArray
    tstops::AbstractArray
    stable::Bool
    solve_time::Float64
    device_states_data::Dict{String, Dict{Symbol, AbstractArray}}
    device_ref_data::Dict{String, Dict{Symbol, Float64}}
end

function AllStatesData(;
    type = "AllStatesData",
    tsteps = [],
    tstops = [],
    stable = false,
    solve_time = 0.0,
    device_states_data = Dict{String, Dict{Symbol, AbstractArray}}(),
    device_ref_data = Dict{String, Dict{Symbol, Float64}}(),
)
    return AllStatesData(
        type,
        tsteps,
        tstops,
        stable,
        solve_time,
        device_states_data,
        device_ref_data,
    )
end

function fill_surrogate_data!(
    data::AllStatesData,
    device_details,
    data_collection_params,
    sim_full,
    results,
    save_indices,
)
    sys = sim_full.sys
    states_data_dict = Dict{String, Dict{Symbol, AbstractArray}}()
    references_data_dict = Dict{String, Dict{Symbol, Float64}}()
    for (device_name, orientation_details) in device_details
        all_components_with_name =
            PSY.get_components_by_name(PSY.Component, sys, device_name)
        exclude_dynamic_injectors =
            filter!(x -> !(typeof(x) <: PSY.DynamicInjection), all_components_with_name)
        @assert length(exclude_dynamic_injectors) == 1
        device = exclude_dynamic_injectors[1]
        _fill_states_data!(
            states_data_dict,
            results,
            device,
            save_indices,
            orientation_details,
            data_collection_params,
            sys,
        )
        _fill_references_data!(
            references_data_dict,
            results,
            device,
            save_indices,
            orientation_details,
            data_collection_params,
            sys,
        )
    end
    data.tstops = unique(results.solution.t)
    data.tsteps = unique(results.solution.t)[save_indices]
    data.stable = true
    data.solve_time = results.time_log[:timed_solve_time]
    data.device_states_data = states_data_dict
    data.device_ref_data = references_data_dict
end

function _fill_states_data!(
    states_data_dict,
    results,
    device::S,
    save_indices,
    orientation_details,
    data_collection_params,
    sys,
) where {S <: PSY.StaticInjection}
    state_trajectories_dict = Dict{Symbol, AbstractArray}()
    dynamic_device_name = PSY.get_name(device)
    dynamic_device = PSY.get_dynamic_injector(device)
    @assert dynamic_device !== nothing
    for s in PSY.get_states(dynamic_device)
        state_trajectories_dict[s] =
            PSID.get_state_series(results, (dynamic_device_name, s))[2][save_indices]
    end
    states_data_dict[dynamic_device_name] = state_trajectories_dict
end

function _fill_references_data!(
    refs_data_dict,
    results,
    device::S,
    save_indices,
    orientation_details,
    data_collection_params,
    sys,
) where {S <: PSY.StaticInjection}
    references = Dict{Symbol, Float64}()
    dynamic_device_name = PSY.get_name(device)
    dynamic_device = PSY.get_dynamic_injector(device)
    @assert dynamic_device !== nothing
    references[:P_ref] = PSY.get_P_ref(dynamic_device)
    references[:ω_ref] = 1.0
    references[:V_ref] = PSY.get_V_ref(dynamic_device)
    references[:Q_ref] = PSY.get_reactive_power(device)
    if typeof(dynamic_device) <: PSY.DynamicGenerator
        references[:eq_p] = PSY.get_eq_p(PSY.get_machine(dynamic_device))
    end
    refs_data_dict[dynamic_device_name] = references
end

function EmptyTrainDataSet(T::Type{AllStatesData})
    return AllStatesData()
end

function generate_empty_plot(T::Type{AllStatesData})
    return PlotlyJS.plot()
end

function add_data_trace!(p, data::AllStatesData; name = "", color = "")
    for (device_name, device_data_dict) in data.device_states_data
        for (s, value) in device_data_dict
            PlotlyJS.add_trace!(
                p,
                PlotlyJS.scatter(;
                    x = data.tsteps,
                    y = value,
                    marker = PlotlyJS.attr(color = color),
                    name = string(s, name),
                ),
                row = 1,
                col = 1,
            )
        end
    end
    return p
end
