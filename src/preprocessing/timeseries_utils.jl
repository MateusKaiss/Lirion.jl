function get_model_times(data)
    time_vecs = []

    if haskey(data, "bus_timeseries") && !isempty(data["bus_timeseries"])
        push!(time_vecs, data["bus_timeseries"]["time"])
    end
    if haskey(data, "prosumer_timeseries") && !isempty(data["prosumer_timeseries"])
        push!(time_vecs, data["prosumer_timeseries"]["time"])
    end
    if haskey(data, "gencost_timeseries") && !isempty(data["gencost_timeseries"])
        push!(time_vecs, data["gencost_timeseries"]["time"])
    end
    if haskey(data, "consumer_timeseries") && !isempty(data["consumer_timeseries"])
        push!(time_vecs, data["consumer_timeseries"]["time"])
    end

    if isempty(time_vecs)
        return Int[]
    end

    all_times = unique(vcat(time_vecs...))
    sort!(all_times)
    return all_times
end

function check_time_alignment(data)
    time_sets = Dict()

    for key in (
        "bus_timeseries",
        "prosumer_timeseries",
        "gencost_timeseries",
        "consumer_timeseries",
    )
        if haskey(data, key) && haskey(data[key], "time") && !isempty(data[key]["time"])
            time_sets[key] = Set(data[key]["time"])
        end
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

function ensure_bus_timeseries!(data::Dict{String, Any}, T::Int)
    ts = get(data, "bus_timeseries", Dict("time" => Int[]))
    if !isempty(ts["time"]) || T == 1 || !haskey(data, "bus")
        return
    end

    @info "Creating bus_timeseries based on static table bus."

    b = data["bus"]
    nb = length(b["bus_i"])

    return data["bus_timeseries"] = Dict(
        "time" => ones(Int, nb),
        "bus" => copy(b["bus_i"]),
        "pd" => copy(b["pd"]),
        "qd" => copy(b["qd"]),
        "gs" => copy(b["gs"]),
        "bs" => copy(b["bs"]),
    )
end

function ensure_gencost_timeseries!(data::Dict{String, Any}, T::Int)
    if haskey(data, "gencost_timeseries") &&
       haskey(data["gencost_timeseries"], "time") &&
       !isempty(data["gencost_timeseries"]["time"])
        return
    end

    if T == 1 || !haskey(data, "gencost")
        return
    end

    @info "Creating gencost_timeseries based on fixed gencost table."

    gencost = data["gencost"]
    ng = length(gencost["cost_coefficients"])
    ts = Dict{String, Vector}()

    ts["time"] = Int[]
    ts["gen_id"] = Int[]
    ts["model"] = Int[]
    ts["startup"] = Float64[]
    ts["shutdown"] = Float64[]
    ts["ncost"] = Int[]
    ts["cost_coefficients"] = Vector{Vector{Float64}}()
    ts["tou"] = Float64[]
    ts["feed_in_tariff"] = Float64[]

    for g in 1:ng
        push!(ts["time"], 1)
        push!(ts["gen_id"], g)
        push!(ts["model"], gencost["model"][g])
        push!(ts["startup"], gencost["startup"][g])
        push!(ts["shutdown"], gencost["shutdown"][g])
        push!(ts["ncost"], gencost["ncost"][g])
        push!(ts["cost_coefficients"], copy(gencost["cost_coefficients"][g]))
        push!(ts["tou"], 1.0)
        push!(ts["feed_in_tariff"], 0.0)
    end

    return data["gencost_timeseries"] = ts
end

function ensure_prosumer_timeseries!(data::Dict{String, Any}, T::Int)
    ts = get(data, "prosumer_timeseries", Dict("time" => Int[]))
    if !isempty(ts["time"]) || T == 1 || !haskey(data, "prosumer")
        return
    end

    @info "Creating prosumer_timeseries based on static table prosumer."

    s = data["prosumer"]
    nb = length(s["bus"])

    return data["prosumer_timeseries"] = Dict(
        "time" => ones(Int, nb),
        "bus" => copy(s["bus"]),
        "p_pv" => copy(s["generation"]),
        "p_load" => copy(s["demand"]),
        "p_curtail" => fill(0.0, nb),
    )
end

function ensure_consumer_timeseries!(data::Dict{String, Any}, T::Int)
    ts = get(data, "consumer_timeseries", Dict("time" => Int[]))
    if !isempty(ts["time"]) || T == 1 || !haskey(data, "consumer")
        return
    end

    @info "Creating consumer_timeseries based on static table consumer."

    c = data["consumer"]
    nb = length(c["bus"])

    data["consumer_timeseries"] = Dict(
        "time" => ones(Int, nb),
        "bus" => copy(c["bus"]),
        "demand" => copy(c["demand"]),
        "q_demand" => copy(c["q_demand"]),
    )

    return
end

function extend_bus_timeseries!(data::Dict{String, Any}, T::Int)
    if !haskey(data, "bus_timeseries")
        return
    end

    ts = data["bus_timeseries"]
    times = ts["time"]
    buses = ts["bus"]
    last_time = maximum(times)
    bus_ids = unique(buses)

    if last_time >= T
        return
    end

    @warn "bus_timeseries has data only up to t=$last_time. Repeating values ​​for t=$(last_time+1):$T automatically."

    for bus in bus_ids
        indices = [
            i for (i, (b, t)) in enumerate(zip(buses, times)) if b == bus && t == last_time
        ]
        if isempty(indices)
            continue
        end
        idx = indices[1]
        for t in (last_time+1):T
            for key in keys(ts)
                push!(ts[key], ts[key][idx])
            end
            ts["time"][end] = t
        end
    end
end

function extend_gencost_timeseries!(data::Dict{String, Any}, T::Int)
    if !haskey(data, "gencost_timeseries") || !haskey(data["gencost_timeseries"], "time")
        return
    end

    ts = data["gencost_timeseries"]
    times = ts["time"]
    gens = ts["gen_id"]
    last_time = maximum(times)
    gen_ids = unique(gens)

    if last_time >= T
        return
    end

    @warn "gencost_timeseries has data only up to t=$last_time. Repeating values for t=$(last_time+1):$T automatically."

    for gen in gen_ids
        indices =
            [i for (i, (g, t)) in enumerate(zip(gens, times)) if g == gen && t == last_time]
        if isempty(indices)
            continue
        end
        idx = indices[1]
        for t in (last_time+1):T
            for key in keys(ts)
                push!(ts[key], ts[key][idx])
            end
            ts["time"][end] = t
        end
    end
end

function extend_prosumer_timeseries!(data::Dict{String, Any}, T::Int)
    if !haskey(data, "prosumer_timeseries")
        return
    end

    ts = data["prosumer_timeseries"]
    times = ts["time"]
    prosumers = ts["bus"]

    if isempty(times)
        return
    end
    last_time = maximum(times)
    prosumer_ids = unique(prosumers)

    if last_time >= T
        return
    end

    @warn "prosumer_timeseries has data only up to t=$last_time. Repeating values ​​for t=$(last_time+1):$T automatically."

    for p in prosumer_ids
        indices = [
            i for (i, (pid, t)) in enumerate(zip(prosumers, times)) if
            pid == p && t == last_time
        ]
        if isempty(indices)
            continue
        end
        idx = indices[1]
        for t in (last_time+1):T
            for key in keys(ts)
                push!(ts[key], ts[key][idx])
            end
            ts["time"][end] = t
        end
    end
end

function extend_consumer_timeseries!(data::Dict{String, Any}, T::Int)
    if !haskey(data, "consumer_timeseries")
        return
    end

    ts = data["consumer_timeseries"]
    times = ts["time"]

    if isempty(times)
        return
    end
    consumers = ts["bus"]
    last_time = maximum(times)
    consumer_ids = unique(consumers)

    if last_time >= T
        return
    end

    @warn "consumer_timeseries has data only up to t=$last_time. Repeating values ​​for t=$(last_time+1):$T automatically."

    for c in consumer_ids
        indices = [
            i for (i, (cid, t)) in enumerate(zip(consumers, times)) if
            cid == c && t == last_time
        ]
        if isempty(indices)
            continue
        end
        idx = indices[1]
        for t in (last_time+1):T
            for key in keys(ts)
                push!(ts[key], ts[key][idx])
            end
            ts["time"][end] = t
        end
    end

    return
end

function sort_timeseries_by_time_and_id!(data::Dict{String, Any})
    id_fields = Dict(
        "bus_timeseries" => "bus",
        "gencost_timeseries" => "gen_id",
        "prosumer_timeseries" => "bus",
        "consumer_timeseries" => "bus",
    )

    for (key, id_field) in id_fields
        if haskey(data, key)
            ts = data[key]
            if haskey(ts, "time") && haskey(ts, id_field)
                order = sortperm([
                    (ts["time"][i], ts[id_field][i]) for i in eachindex(ts["time"])
                ])
                for k in keys(ts)
                    ts[k] = ts[k][order]
                end
            end
        end
    end
end
