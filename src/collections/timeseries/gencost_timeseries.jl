struct GencostTimeseriesData
    has_timeseries::Bool
    time::Vector{Int}
    gen_id::Vector{Int}
    model::Vector{Int}
    startup::Vector{Float64}
    shutdown::Vector{Float64}
    ncost::Vector{Int}
    cost_coefficients::Vector{Vector{Float64}}
    tou::Vector{Float64}
    feed_in_tariff::Vector{Float64}
end

function build_gencost_timeseries_data(data::Dict{String, Any})::GencostTimeseriesData
    if !haskey(data, "gencost_timeseries") ||
       !haskey(data["gencost_timeseries"], "time") ||
       isempty(data["gencost_timeseries"]["time"])
        return GencostTimeseriesData(
            false,
            Int[],
            Int[],
            Int[],
            Float64[],
            Float64[],
            Int[],
            Vector{Vector{Float64}}(),
            Float64[],
            Float64[],
        )
    end

    ts = data["gencost_timeseries"]

    return GencostTimeseriesData(
        true,
        ts["time"],
        ts["gen_id"],
        ts["model"],
        ts["startup"],
        ts["shutdown"],
        ts["ncost"],
        ts["cost_coefficients"],
        get(ts, "tou", ones(Float64, length(ts["time"]))),
        get(ts, "feed_in_tariff", zeros(Float64, length(ts["time"]))),
    )
end
