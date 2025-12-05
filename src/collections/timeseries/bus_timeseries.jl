struct BusTimeseriesData
    has_timeseries::Bool
    time::Vector{Int}
    bus::Vector{Int}
    pd::Vector{Float64}
    qd::Vector{Float64}
    gs::Vector{Float64}
    bs::Vector{Float64}
end

function build_bus_timeseries_data(data::Dict{String, Any})::BusTimeseriesData
    if !haskey(data, "bus_timeseries") ||
       !haskey(data["bus_timeseries"], "time") ||
       isempty(data["bus_timeseries"]["time"])
        return BusTimeseriesData(
            false,
            Int[],
            Int[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
        )
    end

    ts = data["bus_timeseries"]

    return BusTimeseriesData(
        true,
        ts["time"],
        ts["bus"],
        get(ts, "pd", zeros(Float64, length(ts["bus"]))),
        get(ts, "qd", zeros(Float64, length(ts["bus"]))),
        get(ts, "gs", zeros(Float64, length(ts["bus"]))),
        get(ts, "bs", zeros(Float64, length(ts["bus"]))),
    )
end
