#= mutable struct FullSolutionDataParams <: SurrogateDatasetParams
    type::String
end

function FullSolutionDataParams(;type = "FullSolutionDataParams")
    return FullSolutionDataParams(type)
end
 =#
mutable struct FullSolutionData <: SurrogateDataset
    type::String
    tsteps::AbstractArray
    tstops::AbstractArray
    stable::Bool
    solve_time::Float64
    u::Vector{AbstractArray}
end

function FullSolutionData(;
    type = "FullSolutionData",
    tsteps = [],
    tstops = [],
    stable = false,
    solve_time = 0.0,
    u = Vector{AbstractArray}(),
)
    return FullSolutionData(type, tsteps, tstops, stable, solve_time, u)
end

function fill_surrogate_data!(
    data::FullSolutionData,
    device_details,
    data_collection_params,
    sim_full,
    results,
    save_indices,
)
    data.tstops = unique(results.solution.t)
    data.tsteps = unique(results.solution.t)[save_indices]
    data.u = results.solution.u[save_indices]
    data.stable = true
    data.solve_time = results.time_log[:timed_solve_time]
end

function EmptyTrainDataSet(T::Type{FullSolutionData})
    return FullSolutionData()
end

function generate_empty_plot(T::Type{FullSolutionData})
    return PlotlyJS.plot()
end

function add_data_trace!(p, data::FullSolutionData; name = "")
    n_values = length(data.u[1])
    for n in 1:n_values
        y = [u[n] for u in data.u]
        PlotlyJS.add_trace!(
            p,
            PlotlyJS.scatter(; x = data.tsteps, y = y, title = "All values together"),
        )
    end
end
