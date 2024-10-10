mutable struct BusData <: SurrogateDataset
    type::String
    tsteps::AbstractArray
    tstops::AbstractArray
    stable::Bool
    solve_time::Float64
    bus_data::Dict{Int64, Dict{Symbol, AbstractArray}}
end

function BusData(;
    type = "BusData",
    tsteps = [],
    tstops = [],
    stable = false,
    solve_time = 0.0,
    bus_data = Dict{Int64, Dict{Symbol, AbstractArray}}(),
)
    return BusData(type, tsteps, tstops, stable, solve_time, bus_data)
end

function fill_surrogate_data!(
    data::BusData,
    device_details,
    data_collection_params,
    sim_full,
)
    if sim_full.status == PSID.SIMULATION_FINALIZED
        results = PSID.read_results(sim_full)
        if length(data_collection_params.tsave) == 0
            save_indices = 1:length(unique(results.solution.t))
        else
            save_indices = indexin(data_collection_params.tsave, unique(results.solution.t))
        end
        sys = sim_full.sys
        bus_data_dict = Dict{Int64, Dict{Symbol, AbstractArray}}()
        @show device_details
        for (bus_name, _) in device_details
            bus_number = PSY.get_number(PSY.get_component(PSY.ACBus, sys, bus_name))
            _fill_bus_data!(
                bus_data_dict,
                bus_number,
                results,
                save_indices,
                data_collection_params,
                sys,
            )
        end
        data.tstops = unique(results.solution.t)
        data.tsteps = unique(results.solution.t)[save_indices]
        data.stable = true
        data.solve_time = results.time_log[:timed_solve_time]
        data.bus_data = bus_data_dict
    else
        @error "Simulation was not stable, not recording any data for BusData"
    end
end

function _fill_bus_data!(
    bus_data_dict,
    bus_number,
    results,
    save_indices,
    data_collection_params,
    sys,
)
    voltage_trajectories_dict = Dict{Symbol, AbstractArray}()
    vm = PSID.get_voltage_magnitude_series(results, bus_number)[2][save_indices]
    vθ = PSID.get_voltage_angle_series(results, bus_number)[2][save_indices]
    vr = vm .* cos.(vθ)
    vi = vm .* sin.(vθ)
    voltage_trajectories_dict[:vm] = vm
    voltage_trajectories_dict[:vθ] = vθ
    voltage_trajectories_dict[:vr] = vr
    voltage_trajectories_dict[:vi] = vi
    bus_data_dict[bus_number] = voltage_trajectories_dict
end

function EmptyTrainDataSet(T::Type{BusData})
    return BusData()
end

function generate_empty_plot(T::Type{BusData})
    p = PlotlyJS.make_subplots(
        rows = 2,
        cols = 2,
        specs = [
            PlotlyJS.Spec() PlotlyJS.Spec()
            PlotlyJS.Spec() PlotlyJS.Spec()
        ],
        subplot_titles = ["vm" "vθ" "vr" "vi"],
        vertical_spacing = 0.1,
    )
    return p
end

function add_data_trace!(p, data::BusData; name = "", color = "")
    for (bus_number, bus_data) in data.bus_data
        PlotlyJS.add_trace!(
            p,
            PlotlyJS.scatter(;
                x = data.tsteps,
                y = bus_data[:vm],
                marker = PlotlyJS.attr(color = color),
                name = string(bus_number, name),
            ),
            row = 1,
            col = 1,
        )
        PlotlyJS.add_trace!(
            p,
            PlotlyJS.scatter(;
                x = data.tsteps,
                y = bus_data[:vθ],
                marker = PlotlyJS.attr(color = color),
                name = string(bus_number, name),
            ),
            row = 1,
            col = 2,
        )
        PlotlyJS.add_trace!(
            p,
            PlotlyJS.scatter(;
                x = data.tsteps,
                y = bus_data[:vr],
                marker = PlotlyJS.attr(color = color),
                name = string(bus_number, name),
            ),
            row = 2,
            col = 1,
        )
        PlotlyJS.add_trace!(
            p,
            PlotlyJS.scatter(;
                x = data.tsteps,
                y = bus_data[:vi],
                marker = PlotlyJS.attr(color = color),
                name = string(bus_number, name),
            ),
            row = 2,
            col = 2,
        )
    end
    return p
end
