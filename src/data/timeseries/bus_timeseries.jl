# ============================================================
# Bus time-series data
# ============================================================

struct BusTimeseriesData
    has_timeseries::Bool
    time::Vector{Int}
    bus::Vector{Int}
    pd::Vector{Float64}
    qd::Vector{Float64}
    gs::Vector{Float64}
    bs::Vector{Float64}
end

function build_bus_timeseries_data(case_raw::CaseRaw)::BusTimeseriesData
    ts = case_raw.bus_timeseries
    if ts === nothing || isempty(ts.time)
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
    ts = ts::BusTimeseriesRaw{Float64}
    return BusTimeseriesData(true, ts.time, ts.bus, ts.pd, ts.qd, ts.gs, ts.bs)
end
