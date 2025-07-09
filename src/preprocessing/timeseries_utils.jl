function check_time_alignment(data)
    time_sets = Dict()

    for key in ("bus_timeseries", "prosumer_timeseries", "gencost_timeseries", "consumer_timeseries")
        if haskey(data, key)
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

function get_model_times(data)
    time_vecs = []

    if haskey(data, "bus_timeseries")
        push!(time_vecs, data["bus_timeseries"]["time"])
    end
    if haskey(data, "prosumer_timeseries")
        push!(time_vecs, data["prosumer_timeseries"]["time"])
    end
    if haskey(data, "gencost_timeseries")
        push!(time_vecs, data["gencost_timeseries"]["time"])
    end
    if haskey(data, "consumer_timeseries")
        push!(time_vecs, data["consumer_timeseries"]["time"])
    end

    if isempty(time_vecs)
        return Int[]
    end

    all_times = unique(vcat(time_vecs...))
    sort!(all_times)
    return all_times
end

function ensure_prosumer_timeseries!(data::Dict{String, Any}, T::Int)
    ts = get(data, "prosumer_timeseries", Dict("time" => Int[]))
    if !isempty(ts["time"]) || T == 1 || !haskey(data, "prosumer")
        return
    end

    @info "Criando prosumer_timeseries com base na tabela estática prosumer."

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

    @info "Criando consumer_timeseries com base na tabela estática consumer."

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

function ensure_gencost_timeseries!(data::Dict{String, Any}, T::Int)
    if haskey(data, "gencost_timeseries") && !isempty(data["gencost_timeseries"]["time"])
        return
    end
    if T == 1 || !haskey(data, "gencost")
        return
    end

    @info "Criando gencost_timeseries com base na tabela gencost fixa."

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

function ensure_bus_timeseries!(data::Dict{String, Any}, T::Int)
    ts = get(data, "bus_timeseries", Dict("time" => Int[]))
    if !isempty(ts["time"]) || T == 1 || !haskey(data, "bus")
        return
    end

    @info "Criando bus_timeseries com base na tabela estática bus."

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

    @warn "bus_timeseries tem dados apenas até t=$last_time. Repetindo valores para t=$(last_time+1):$T automaticamente."

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
                push!(ts[key], ts[key][idx])  # copia o valor da linha
            end
            ts["time"][end] = t  # atualiza o tempo
        end
    end
end

function extend_gencost_timeseries!(data::Dict{String, Any}, T::Int)
    # Caso 1: não existe gencost_timeseries
    if !haskey(data, "gencost_timeseries")
        @info "gencost_timeseries não encontrado. Criando com base em gencost fixo."

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

        for t in 1:T
            for g in 1:ng
                push!(ts["time"], t)
                push!(ts["gen_id"], g)
                push!(ts["model"], gencost["model"][g])
                push!(ts["startup"], gencost["startup"][g])
                push!(ts["shutdown"], gencost["shutdown"][g])
                push!(ts["ncost"], gencost["ncost"][g])
                push!(ts["cost_coefficients"], gencost["cost_coefficients"][g])
                push!(ts["tou"], 1.0)
                push!(ts["feed_in_tariff"], 0.0)
            end
        end

        data["gencost_timeseries"] = ts
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

    @warn "gencost_timeseries tem dados apenas até t=$last_time. Repetindo valores para t=$(last_time+1):$T automaticamente."

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
    last_time = maximum(times)
    prosumer_ids = unique(prosumers)

    if last_time >= T
        return
    end

    @warn "prosumer_timeseries tem dados apenas até t=$last_time. Repetindo valores para t=$(last_time+1):$T automaticamente."

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
    consumers = ts["bus"]
    last_time = maximum(times)
    consumer_ids = unique(consumers)

    if last_time >= T
        return
    end

    @warn "consumer_timeseries tem dados apenas até t=$last_time. Repetindo valores para t=$(last_time+1):$T automaticamente."

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
