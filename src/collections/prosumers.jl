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
end

function get_prosumer_parameters(raw::CaseRaw)
    T = raw.T
    has_ts = raw.prosumer_timeseries !== nothing && !isempty(raw.prosumer_timeseries.time)

    if has_ts && T > 1
        ts = raw.prosumer_timeseries::ProsumerTimeseriesRaw{Float64}
        nb = length(unique(ts.bus))
        s = raw.prosumer
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

    s = raw.prosumer
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

function build_prosumer_data(raw::CaseRaw)::ProsumerData
    p = get_prosumer_parameters(raw)
    return ProsumerData(
        p.is_timeseries,
        p.number_of_prosumers,
        p.time,
        p.bus_index,
        p.generation,
        p.demand,
        p.q_generation,
        p.q_demand,
        p.curtailment,
        p.cost_injection,
        p.inverter_capacity,
    )
end
