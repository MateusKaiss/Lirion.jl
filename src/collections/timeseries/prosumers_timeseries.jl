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

function build_prosumer_timeseries_data(raw::CaseRaw)::ProsumerTimeseriesData
    ts = raw.prosumer_timeseries
    if ts === nothing || isempty(ts.time)
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
    ts = ts::ProsumerTimeseriesRaw{Float64}
    return ProsumerTimeseriesData(
        true,
        ts.time,
        ts.bus,
        ts.p_pv,
        ts.p_load,
        fill(0.0, length(ts.bus)),  # q_generation (não existe no raw atual)
        fill(0.0, length(ts.bus)),  # q_demand
        ts.p_curtail,
    )
end
