mutable struct FullSolutionData <: SurrogateDataset
    type::String
    tsteps::AbstractArray
    tstops::AbstractArray
    stable::Bool
    solve_time::Float64
    psid_result::Union{Nothing, PSID.SimulationResults}
end

function FullSolutionData(;
    type = "FullSolutionData",
    tsteps = [],
    tstops = [],
    stable = false,
    solve_time = 0.0,
    psid_result = nothing, # Vector{AbstractArray}(),
)
    return FullSolutionData(type, tsteps, tstops, stable, solve_time, psid_result)
end

function fill_surrogate_data!(
    data::FullSolutionData,
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
        data.tstops = unique(results.solution.t)
        data.tsteps = unique(results.solution.t)[save_indices]
        data.psid_result = results #  results.solution.u[save_indices]
        data.stable = true
        data.solve_time = results.time_log[:timed_solve_time]
    else
        @error "Simulation was not stable, not recording any data for FullSolutionData"
    end
end

function EmptyTrainDataSet(T::Type{FullSolutionData})
    return FullSolutionData()
end

function generate_empty_plot(T::Type{FullSolutionData})
    return PlotlyJS.plot()
end

function add_data_trace!(p, data::FullSolutionData; name = "", color = "") end
