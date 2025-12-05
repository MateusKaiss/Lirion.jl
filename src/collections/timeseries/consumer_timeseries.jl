struct ConsumerTimeseriesData
    has_timeseries::Bool
    time::Vector{Int}
    bus::Vector{Int}
    demand::Vector{Float64}
    q_demand::Vector{Float64}
end

function build_consumer_timeseries_data(data::Dict{String, Any})::ConsumerTimeseriesData
    if !haskey(data, "consumer_timeseries") ||
       !haskey(data["consumer_timeseries"], "time") ||
       isempty(data["consumer_timeseries"]["time"])
        return ConsumerTimeseriesData(
            false,
            Int[],
            Int[],
            Float64[],
            Float64[],
        )
    end

    ts = data["consumer_timeseries"]

    return ConsumerTimeseriesData(
        true,
        ts["time"],
        ts["bus"],
        get(ts, "demand", zeros(Float64, length(ts["bus"]))),
        get(ts, "q_demand", zeros(Float64, length(ts["bus"]))),
    )
end
