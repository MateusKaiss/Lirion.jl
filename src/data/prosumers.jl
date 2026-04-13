struct ProsumerData
    is_timeseries::Bool
    number_of_prosumers::Int
    bus_index::Vector{Int}
    generation::Matrix{Float64}
    demand::Matrix{Float64}
    q_generation::Matrix{Float64}
    q_demand::Matrix{Float64}
    curtailment::Matrix{Float64}
    cost_injection::Vector{Float64}
    inverter_capacity::Vector{Float64}
    inverter_capacity_sq::Vector{Float64}
end

function _pivot_ts_to_matrix(
    time::Vector{Int},
    bus::Vector{Int},
    values::Vector{Float64},
    bus_to_local::Dict{Int, Int},
    T::Int,
    np::Int,
)
    out = zeros(Float64, np, T)
    @inbounds for i in eachindex(time)
        t = time[i]
        (t < 1 || t > T) && continue
        p = get(bus_to_local, bus[i], 0)
        p == 0 && continue
        out[p, t] = values[i]
    end
    return out
end

function build_prosumer_data(case_raw::CaseRaw)::ProsumerData
    T = max(case_raw.T, 1)
    s = case_raw.prosumer
    has_ts =
        case_raw.prosumer_timeseries !== nothing &&
        !isempty(case_raw.prosumer_timeseries.time)

    if s === nothing
        empty_mat = zeros(Float64, 0, T)
        return ProsumerData(
            false, 0, Int[],
            empty_mat, empty_mat, empty_mat, empty_mat, empty_mat,
            Float64[], Float64[], Float64[],
        )
    end

    np = length(s.bus)
    bus_to_local = Dict{Int, Int}()
    sizehint!(bus_to_local, np)
    @inbounds for p in 1:np
        bus_to_local[s.bus[p]] = p
    end

    inv_cap = copy(s.inverter_capacity)
    inv_cap_sq = inv_cap .^ 2
    cost_inj = copy(s.cost_injection)

    if has_ts && T > 1
        ts = case_raw.prosumer_timeseries::ProsumerTimeseriesRaw{Float64}
        generation =
            _pivot_ts_to_matrix(ts.time, ts.bus, ts.p_pv, bus_to_local, T, np)
        demand =
            _pivot_ts_to_matrix(ts.time, ts.bus, ts.p_load, bus_to_local, T, np)
        curtail =
            _pivot_ts_to_matrix(ts.time, ts.bus, ts.p_curtail, bus_to_local, T, np)
        # Q is not in the prosumer timeseries table — broadcast static values.
        q_generation = repeat(reshape(s.q_generation, np, 1), 1, T)
        q_demand = repeat(reshape(s.q_demand, np, 1), 1, T)

        return ProsumerData(
            true, np, copy(s.bus),
            generation, demand, q_generation, q_demand, curtail,
            cost_inj, inv_cap, inv_cap_sq,
        )
    end

    # Static mode — broadcast scalars to (np, T).
    generation = repeat(reshape(s.generation, np, 1), 1, T)
    demand = repeat(reshape(s.demand, np, 1), 1, T)
    q_generation = repeat(reshape(s.q_generation, np, 1), 1, T)
    q_demand = repeat(reshape(s.q_demand, np, 1), 1, T)
    curtail = zeros(Float64, np, T)

    return ProsumerData(
        false, np, copy(s.bus),
        generation, demand, q_generation, q_demand, curtail,
        cost_inj, inv_cap, inv_cap_sq,
    )
end
