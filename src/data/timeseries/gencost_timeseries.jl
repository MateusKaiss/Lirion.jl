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

function build_gencost_timeseries_data(raw::CaseRaw)::GencostTimeseriesData
    ts = raw.gencost_timeseries
    if ts === nothing || isempty(ts.time)
        return GencostTimeseriesData(
            false, Int[], Int[], Int[], Float64[], Float64[], Int[],
            Vector{Vector{Float64}}(), Float64[], Float64[],
        )
    end
    ts = ts::GencostTimeseriesRaw{Float64}
    return GencostTimeseriesData(
        true, ts.time, ts.gen_id, ts.model, ts.startup, ts.shutdown,
        ts.ncost, ts.cost_coefficients, ts.tou, ts.feed_in_tariff,
    )
end
