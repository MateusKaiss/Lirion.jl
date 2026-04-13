struct ConsumerTimeseriesData
    has_timeseries::Bool
    time::Vector{Int}
    bus::Vector{Int}
    demand::Vector{Float64}
    q_demand::Vector{Float64}
end

function build_consumer_timeseries_data(case_raw::CaseRaw)::ConsumerTimeseriesData
    ts = case_raw.consumer_timeseries
    if ts === nothing || isempty(ts.time)
        return ConsumerTimeseriesData(false, Int[], Int[], Float64[], Float64[])
    end
    ts = ts::ConsumerTimeseriesRaw{Float64}
    return ConsumerTimeseriesData(true, ts.time, ts.bus, ts.demand, ts.q_demand)
end
