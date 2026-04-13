function get_model_times(case::CaseRaw)
    time_vecs = Vector{Vector{Int}}()

    if case.bus_timeseries !== nothing && !isempty(case.bus_timeseries.time)
        push!(time_vecs, case.bus_timeseries.time)
    end
    if case.prosumer_timeseries !== nothing && !isempty(case.prosumer_timeseries.time)
        push!(time_vecs, case.prosumer_timeseries.time)
    end
    if case.gencost_timeseries !== nothing && !isempty(case.gencost_timeseries.time)
        push!(time_vecs, case.gencost_timeseries.time)
    end
    if case.consumer_timeseries !== nothing && !isempty(case.consumer_timeseries.time)
        push!(time_vecs, case.consumer_timeseries.time)
    end

    isempty(time_vecs) && return Int[]
    all_times = unique(vcat(time_vecs...))
    sort!(all_times)
    return all_times
end

function check_time_alignment(case::CaseRaw)
    time_sets = Dict{String, Set{Int}}()

    if case.bus_timeseries !== nothing && !isempty(case.bus_timeseries.time)
        time_sets["bus_timeseries"] = Set(case.bus_timeseries.time)
    end
    if case.prosumer_timeseries !== nothing && !isempty(case.prosumer_timeseries.time)
        time_sets["prosumer_timeseries"] = Set(case.prosumer_timeseries.time)
    end
    if case.gencost_timeseries !== nothing && !isempty(case.gencost_timeseries.time)
        time_sets["gencost_timeseries"] = Set(case.gencost_timeseries.time)
    end
    if case.consumer_timeseries !== nothing && !isempty(case.consumer_timeseries.time)
        time_sets["consumer_timeseries"] = Set(case.consumer_timeseries.time)
    end

    if length(time_sets) > 1
        ref_key, ref_times = first(time_sets)
        for (key, ts) in time_sets
            if ts != ref_times
                @warn "Table $key has different time points from $ref_key. Time vectors will be unified using full union. ($key: $(sort(collect(ts))) vs $ref_key: $(sort(collect(ref_times))))"
            end
        end
    end
end

function ensure_bus_timeseries!(case::CaseRaw, T::Int)
    if case.bus_timeseries !== nothing && !isempty(case.bus_timeseries.time)
        return
    end
    if T == 1 || case.bus === nothing
        return
    end

    @info "Creating bus_timeseries based on static table bus."

    b = case.bus
    nb = length(b.bus_i)

    case.bus_timeseries = BusTimeseriesRaw{Float64}(
        ones(Int, nb),
        copy(b.bus_i),
        copy(Float64.(b.pd)),
        copy(Float64.(b.qd)),
        copy(Float64.(b.gs)),
        copy(Float64.(b.bs)),
    )
    return
end

function ensure_gencost_timeseries!(case::CaseRaw, T::Int)
    if case.gencost_timeseries !== nothing && !isempty(case.gencost_timeseries.time)
        return
    end
    if T == 1 || case.gencost === nothing
        return
    end

    @info "Creating gencost_timeseries based on fixed gencost table."

    g = case.gencost
    ng = length(g.cost_coefficients)

    time = Int[]
    gen_id = Int[]
    model = Int[]
    startup = Float64[]
    shutdown = Float64[]
    ncost = Int[]
    cost_coefficients = Vector{Vector{Float64}}()
    tou = Float64[]
    feed_in_tariff = Float64[]

    for i in 1:ng
        push!(time, 1)
        push!(gen_id, i)
        push!(model, g.model[i])
        push!(startup, Float64(g.startup[i]))
        push!(shutdown, Float64(g.shutdown[i]))
        push!(ncost, g.ncost[i])
        push!(cost_coefficients, copy(Float64.(g.cost_coefficients[i])))
        push!(tou, 1.0)
        push!(feed_in_tariff, 0.0)
    end

    case.gencost_timeseries = GencostTimeseriesRaw{Float64}(
        time, gen_id, model, startup, shutdown, ncost, cost_coefficients, tou,
        feed_in_tariff,
    )
    return
end

function ensure_prosumer_timeseries!(case::CaseRaw, T::Int)
    if case.prosumer_timeseries !== nothing && !isempty(case.prosumer_timeseries.time)
        return
    end
    if T == 1 || case.prosumer === nothing
        return
    end

    @info "Creating prosumer_timeseries based on static table prosumer."

    s = case.prosumer
    nb = length(s.bus)

    case.prosumer_timeseries = ProsumerTimeseriesRaw{Float64}(
        ones(Int, nb),
        copy(s.bus),
        copy(Float64.(s.generation)),
        copy(Float64.(s.demand)),
        fill(0.0, nb),
    )
    return
end

function ensure_consumer_timeseries!(case::CaseRaw, T::Int)
    if case.consumer_timeseries !== nothing && !isempty(case.consumer_timeseries.time)
        return
    end
    if T == 1 || case.consumer === nothing
        return
    end

    @info "Creating consumer_timeseries based on static table consumer."

    c = case.consumer
    nb = length(c.bus)

    case.consumer_timeseries = ConsumerTimeseriesRaw{Float64}(
        ones(Int, nb),
        copy(c.bus),
        copy(Float64.(c.demand)),
        copy(Float64.(c.q_demand)),
    )
    return
end

function _extend_generic_timeseries!(
    times::Vector{Int},
    ids::Vector{Int},
    cols::Vector{AbstractVector},
    T::Int,
    label::String,
)
    isempty(times) && return

    last_time = maximum(times)
    last_time >= T && return

    @warn "$label has data only up to t=$last_time. Repeating values for t=$(last_time+1):$T automatically."

    unique_ids = unique(ids)
    for id in unique_ids
        idxs = [
            i for
            (i, (ii, tt)) in enumerate(zip(ids, times)) if ii == id && tt == last_time
        ]
        isempty(idxs) && continue
        idx = idxs[1]

        for t in (last_time+1):T
            push!(times, t)
            push!(ids, id)
            for col in cols
                push!(col, col[idx])
            end
        end
    end
    return
end

function extend_bus_timeseries!(case::CaseRaw, T::Int)
    case.bus_timeseries === nothing && return
    ts = case.bus_timeseries
    cols = AbstractVector[ts.pd, ts.qd, ts.gs, ts.bs]
    return _extend_generic_timeseries!(ts.time, ts.bus, cols, T, "bus_timeseries")
end

function extend_gencost_timeseries!(case::CaseRaw, T::Int)
    case.gencost_timeseries === nothing && return
    ts = case.gencost_timeseries
    cols = AbstractVector[
        ts.model,
        ts.startup,
        ts.shutdown,
        ts.ncost,
        ts.cost_coefficients,
        ts.tou,
        ts.feed_in_tariff,
    ]
    return _extend_generic_timeseries!(ts.time, ts.gen_id, cols, T, "gencost_timeseries")
end

function extend_prosumer_timeseries!(case::CaseRaw, T::Int)
    case.prosumer_timeseries === nothing && return
    ts = case.prosumer_timeseries
    cols = AbstractVector[ts.p_pv, ts.p_load, ts.p_curtail]
    return _extend_generic_timeseries!(ts.time, ts.bus, cols, T, "prosumer_timeseries")
end

function extend_consumer_timeseries!(case::CaseRaw, T::Int)
    case.consumer_timeseries === nothing && return
    ts = case.consumer_timeseries
    cols = AbstractVector[ts.demand, ts.q_demand]
    return _extend_generic_timeseries!(ts.time, ts.bus, cols, T, "consumer_timeseries")
end

function sort_timeseries_by_time_and_id!(case::CaseRaw)
    if case.bus_timeseries !== nothing
        ts = case.bus_timeseries
        order = sortperm([(ts.time[i], ts.bus[i]) for i in eachindex(ts.time)])
        case.bus_timeseries = slice_bus_ts(ts, order)
    end
    if case.gencost_timeseries !== nothing
        ts = case.gencost_timeseries
        order = sortperm([(ts.time[i], ts.gen_id[i]) for i in eachindex(ts.time)])
        case.gencost_timeseries = slice_gencost_ts(ts, order)
    end
    if case.prosumer_timeseries !== nothing
        ts = case.prosumer_timeseries
        order = sortperm([(ts.time[i], ts.bus[i]) for i in eachindex(ts.time)])
        case.prosumer_timeseries = slice_prosumer_ts(ts, order)
    end
    if case.consumer_timeseries !== nothing
        ts = case.consumer_timeseries
        order = sortperm([(ts.time[i], ts.bus[i]) for i in eachindex(ts.time)])
        case.consumer_timeseries = slice_consumer_ts(ts, order)
    end
    return
end
