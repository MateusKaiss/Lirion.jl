
function prepare_admm_partitions(data::Dict{String, Any})
    bus_data = data["bus"]
    zones = bus_data["zone"]
    bus_ids = bus_data["bus_i"]

    unique_zones = unique(zones)
    partitions = Dict{Int, Dict{String, Any}}()

    for zone in unique_zones
        zone_mask = findall(==(zone), zones)
        zone_bus_ids = bus_ids[zone_mask]

        subdata = Dict{String, Any}()
        subdata["baseMVA"] = data["baseMVA"]
        subdata["T"] = data["T"]
        subdata["source_type"] = data["source_type"]
        subdata["source_version"] = data["source_version"]
        subdata["name"] = "$(data["name"])_zone_$zone"
        subdata["bus"] =
            Dict(k => v[zone_mask] for (k, v) in data["bus"] if length(v) == length(zones))

        for key in ["gen", "prosumer", "consumer", "storage"]
            if haskey(data, key)
                table = data[key]
                bus_key =
                    haskey(table, "bus") ? "bus" :
                    haskey(table, "gen_bus") ? "gen_bus" :
                    haskey(table, "storage_bus") ? "storage_bus" :
                    nothing

                if bus_key !== nothing
                    keep = findall(
                        i -> table[bus_key][i] in zone_bus_ids,
                        1:length(table[bus_key]),
                    )
                    if !isempty(keep)
                        subdata[key] = Dict(
                            k => v[keep] for
                            (k, v) in table if length(v) == length(table[bus_key])
                        )
                    else
                        subdata[key] = Dict(k => Any[] for k in keys(table))
                    end
                else
                    subdata[key] = Dict(k => Any[] for k in keys(table))
                end
            end
        end

        if haskey(data, "branch")
            f_bus = data["branch"]["f_bus"]
            t_bus = data["branch"]["t_bus"]
            keep = findall(
                i -> f_bus[i] in zone_bus_ids && t_bus[i] in zone_bus_ids,
                eachindex(f_bus),
            )
            subdata["branch"] = Dict(
                k => v[keep] for (k, v) in data["branch"] if length(v) == length(f_bus)
            )
        end

        for key in ["bus_timeseries", "prosumer_timeseries", "consumer_timeseries"]
            if haskey(data, key)
                ts = data[key]
                id_field = "bus"
                keep = findall(i -> ts[id_field][i] in zone_bus_ids, 1:length(ts[id_field]))
                subdata[key] =
                    Dict(k => v[keep] for (k, v) in ts if length(v) == length(ts[id_field]))
            end
        end

        if haskey(data, "gen")
            gen = data["gen"]
            bus_field = haskey(gen, "bus") ? "bus" : "gen_bus"
            if haskey(gen, bus_field)
                keep = findall(
                    i -> gen[bus_field][i] in zone_bus_ids,
                    1:length(gen[bus_field]),
                )
                if !isempty(keep)
                    subdata["gen"] = Dict(
                        k => v[keep] for
                        (k, v) in gen if length(v) == length(gen[bus_field])
                    )
                else
                    subdata["gen"] = Dict(k => Any[] for k in keys(gen))
                end

                if haskey(data, "gencost")
                    subdata["gencost"] = Dict(
                        k => v[keep] for
                        (k, v) in data["gencost"] if length(v) == length(gen[bus_field])
                    )
                end

                if haskey(data, "gencost_timeseries")
                    subdata["gencost_timeseries"] = Dict(
                        k => v[keep] for
                        (k, v) in data["gencost_timeseries"] if
                        length(v) == length(gen[bus_field])
                    )
                end
            else
                subdata["gen"] = Dict(k => Any[] for k in keys(gen))
                subdata["gencost"] = Dict{String, Vector{Any}}()
                subdata["gencost_timeseries"] = Dict{String, Vector{Any}}()
            end
        else
            subdata["gen"] = Dict{String, Vector{Any}}()
            subdata["gencost"] = Dict{String, Vector{Any}}()
            subdata["gencost_timeseries"] = Dict{String, Vector{Any}}()
        end

        preprocess_case_data!(subdata, DistributedOPFModel)
        partitions[zone] = subdata
    end

    return partitions
end

function find_shared_buses(data::Dict{String, Any}, partitions::Dict{Int, Dict{String, Any}})
        bus_to_zone = Dict{Int, Int}()
    for (zone_id, pdata) in partitions
        for b in pdata["bus"]["bus_i"]
            bus_to_zone[b] = zone_id
        end
    end

    shared_buses = Set{Int}()
    f_bus = data["branch"]["f_bus"]
    t_bus = data["branch"]["t_bus"]

    for (f, t) in zip(f_bus, t_bus)
        if haskey(bus_to_zone, f) && haskey(bus_to_zone, t)
            if bus_to_zone[f] != bus_to_zone[t]
                push!(shared_buses, f)
                push!(shared_buses, t)
            end
        end
    end

    return collect(shared_buses)
end