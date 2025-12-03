struct ProsumerTimeseriesData
    has_timeseries::Bool
    time::Vector{Int}
    bus::Vector{Int}
    p_pv::Vector{Float64}
    p_load::Vector{Float64}
    q_generation::Vector{Float64}
    q_demand::Vector{Float64}
    p_curtail::Vector{Float64}
end

function build_prosumer_timeseries_data(data::Dict{String,Any})::ProsumerTimeseriesData
    if !haskey(data, "prosumer_timeseries") ||
       !haskey(data["prosumer_timeseries"], "time") ||
       isempty(data["prosumer_timeseries"]["time"])

        return ProsumerTimeseriesData(
            false,
            Int[],
            Int[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
        )
    end

    ts = data["prosumer_timeseries"]

    return ProsumerTimeseriesData(
        true,
        ts["time"],
        ts["bus"],
        get(ts, "p_pv",      zeros(Float64, length(ts["bus"]))),
        get(ts, "p_load",    zeros(Float64, length(ts["bus"]))),
        get(ts, "q_generation", zeros(Float64, length(ts["bus"]))),
        get(ts, "q_demand",     zeros(Float64, length(ts["bus"]))),
        get(ts, "p_curtail",    zeros(Float64, length(ts["bus"]))),
    )
end