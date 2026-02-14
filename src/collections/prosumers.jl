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

function get_prosumer_parameters(data)
    has_ts = haskey(data, "prosumer_timeseries")
    T = get(data, "T", 1)
    s = get(data, "prosumer", Dict())

    if has_ts
        ts = data["prosumer_timeseries"]
        if T > 1 && length(ts["time"]) > 0
            nb = length(unique(ts["bus"]))
            return (
                is_timeseries = true,
                number_of_prosumers = nb,
                time = ts["time"],
                bus_index = ts["bus"],
                q_generation = get(ts, "q_generation", fill(0.0, length(ts["bus"]))),
                generation = get(ts, "p_pv", fill(0.0, length(ts["bus"]))),
                q_demand = get(ts, "q_demand", fill(0.0, length(ts["bus"]))),
                demand = get(ts, "p_load", fill(0.0, length(ts["bus"]))),
                curtailment = get(ts, "p_curtail", fill(0.0, length(ts["bus"]))),
                cost_injection = get(s, "cost_injection", fill(0.0, nb)),
                inverter_capacity = get(s, "inverter_capacity", fill(Inf, nb)),
            )
        end
    end

    nb = get(s, "bus", Int[]) |> length
    return (
        is_timeseries = false,
        number_of_prosumers = nb,
        time = ones(Int, nb),
        bus_index = get(s, "bus", Int[]),
        generation = get(s, "generation", fill(0.0, nb)),
        demand = get(s, "demand", fill(0.0, nb)),
        q_generation = get(s, "q_generation", fill(0.0, nb)),
        q_demand = get(s, "q_demand", fill(0.0, nb)),
        curtailment = fill(0.0, nb),
        cost_injection = get(s, "cost_injection", fill(0.0, nb)),
        inverter_capacity = get(s, "inverter_capacity", fill(Inf, nb)),
    )
end

function build_prosumer_data(data::Dict{String, Any})::ProsumerData
    params = get_prosumer_parameters(data)

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
    )
end
