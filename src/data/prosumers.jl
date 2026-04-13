struct ProsumerData
    is_timeseries::Bool
    number_of_prosumers::Int
    time::Vector{Int}
    bus_index::Vector{Int}
    generation::Vector{Float64}
    demand::Vector{Float64}
    q_generation::Vector{Float64}
    q_demand::Vector{Float64}
    curtailment::Vector{Float64}
    cost_injection::Vector{Float64}
    inverter_capacity::Vector{Float64}
    inverter_capacity_sq::Vector{Float64}
end

function get_prosumer_parameters(case_raw::CaseRaw)
    T = case_raw.T
    has_ts =
        case_raw.prosumer_timeseries !== nothing &&
        !isempty(case_raw.prosumer_timeseries.time)

    if has_ts && T > 1
        ts = case_raw.prosumer_timeseries::ProsumerTimeseriesRaw{Float64}
        nb = length(unique(ts.bus))
        s = case_raw.prosumer
        return (
            is_timeseries = true,
            number_of_prosumers = nb,
            time = ts.time,
            bus_index = ts.bus,
            q_generation = fill(0.0, length(ts.bus)),
            generation = ts.p_pv,
            q_demand = fill(0.0, length(ts.bus)),
            demand = ts.p_load,
            curtailment = ts.p_curtail,
            cost_injection = s === nothing ? fill(0.0, nb) : getfield(s, :cost_injection),
            inverter_capacity = s === nothing ? fill(Inf, nb) :
                                getfield(s, :inverter_capacity),
        )
    end

    s = case_raw.prosumer
    nb = s === nothing ? 0 : length(s.bus)
    return (
        is_timeseries = false,
        number_of_prosumers = nb,
        time = ones(Int, nb),
        bus_index = s === nothing ? Int[] : s.bus,
        generation = s === nothing ? Float64[] : s.generation,
        demand = s === nothing ? Float64[] : s.demand,
        q_generation = s === nothing ? Float64[] : s.q_generation,
        q_demand = s === nothing ? Float64[] : s.q_demand,
        curtailment = fill(0.0, nb),
        cost_injection = s === nothing ? Float64[] : s.cost_injection,
        inverter_capacity = s === nothing ? Float64[] : s.inverter_capacity,
    )
end

function build_prosumer_data(case_raw::CaseRaw)::ProsumerData
    params = get_prosumer_parameters(case_raw)
    inverter_cap_sq = params.inverter_capacity .^ 2
    return ProsumerData(
        params.is_timeseries,
        params.number_of_prosumers,
        params.time,
        params.bus_index,
        params.generation,
        params.demand,
        params.q_generation,
        params.q_demand,
        params.curtailment,
        params.cost_injection,
        params.inverter_capacity,
        inverter_cap_sq,
    )
end
