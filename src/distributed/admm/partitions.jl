function prepare_admm_partitions(data::Dict{String, Any})
    bus = data["bus"]
    zones = bus["zone"]; bus_ids = bus["bus_i"]
    unique_zones = unique(zones)

    partitions = Dict{Int, Dict{String, Any}}()
    zone_map = Dict(bus_ids[i] => zones[i] for i in eachindex(bus_ids))

    for zone in unique_zones
        zone_mask = findall(==(zone), zones)
        zone_bus_ids = bus_ids[zone_mask]

        extra = Set{Int}()
        if haskey(data, "branch")
            f_bus = data["branch"]["f_bus"]; t_bus = data["branch"]["t_bus"]
            @inbounds for i in eachindex(f_bus)
                if (f_bus[i] in zone_bus_ids) ⊻ (t_bus[i] in zone_bus_ids)
                    push!(extra, (f_bus[i] in zone_bus_ids) ? t_bus[i] : f_bus[i])
                end
            end
        end
        extra_bus_ids = collect(extra)
        all_bus_ids = union(zone_bus_ids, extra_bus_ids)
        bus_mask = findall(b -> b in all_bus_ids, bus_ids)

        sub = Dict{String, Any}()
        sub["baseMVA"] = data["baseMVA"]
        sub["T"] = get(data, "T", 1)
        sub["source_type"] = get(data, "source_type", "matpower")
        sub["source_version"] = get(data, "source_version", "")
        sub["name"] = "$(get(data, "name", "case"))_zone_$zone"
        sub["bus"] = Dict(k => v[bus_mask] for (k, v) in data["bus"] if length(v)==length(bus_ids))

        is_virtual = [b ∈ extra for b in sub["bus"]["bus_i"]]
        sub["bus"]["is_virtual"] = is_virtual
        for (i,b) in enumerate(sub["bus"]["bus_i"])
            if !(b in zone_bus_ids)
                sub["bus"]["pd"][i] = 0.0
                sub["bus"]["qd"][i] = 0.0
                sub["bus"]["bus_type"][i] = 1
            end
        end

        for key in ("gen","prosumer","consumer","storage")
            if haskey(data, key)
                tbl = data[key]
                bus_key = haskey(tbl,"bus") ? "bus" :
                          haskey(tbl,"gen_bus") ? "gen_bus" :
                          haskey(tbl,"storage_bus") ? "storage_bus" : nothing
                if bus_key !== nothing
                    keep = findall(i -> tbl[bus_key][i] in zone_bus_ids, eachindex(tbl[bus_key]))
                    sub[key] = !isempty(keep) ?
                        Dict(k => v[keep] for (k,v) in tbl if length(v)==length(tbl[bus_key])) :
                        Dict(k => Any[] for k in keys(tbl))
                else
                    sub[key] = Dict(k => Any[] for k in keys(tbl))
                end
            end
        end

        if haskey(data, "branch")
            f_bus = data["branch"]["f_bus"]; t_bus = data["branch"]["t_bus"]
            keep = findall(i -> (f_bus[i] in zone_bus_ids) || (t_bus[i] in zone_bus_ids),
                           eachindex(f_bus))

            sub["branch"] = Dict(k => v[keep] for (k,v) in data["branch"] if length(v)==length(f_bus))
            sub["branch_full"] = data["branch"]
            sub["branch_local_indices"] = keep
            edge_to_local = Dict{Int,Int}()
            @inbounds for (local_idx, global_idx) in enumerate(keep)
                edge_to_local[global_idx] = local_idx
            end
            sub["edge_to_local"] = edge_to_local
        end

        for key in ("bus_timeseries","prosumer_timeseries","consumer_timeseries")
            if haskey(data, key)
                ts = data[key]
                keep = findall(i -> ts["bus"][i] in zone_bus_ids, eachindex(ts["bus"]))
                sub[key] = Dict(k => v[keep] for (k,v) in ts if length(v)==length(ts["bus"]))
            end
        end

        if haskey(data,"gen")
            if haskey(data,"gencost")
                gen = sub["gen"]; bus_key = haskey(data["gen"],"bus") ? "bus" : "gen_bus"
                keep = findall(i -> data["gen"][bus_key][i] in zone_bus_ids, eachindex(data["gen"][bus_key]))
                sub["gencost"] = Dict(k => v[keep] for (k,v) in data["gencost"] if length(v)==length(data["gen"][bus_key]))
            end
            if haskey(data,"gencost_timeseries")
                gen = sub["gen"]; bus_key = haskey(data["gen"],"bus") ? "bus" : "gen_bus"
                keep = findall(i -> data["gen"][bus_key][i] in zone_bus_ids, eachindex(data["gen"][bus_key]))
                sub["gencost_timeseries"] = Dict(k => v[keep] for (k,v) in data["gencost_timeseries"] if length(v)==length(data["gen"][bus_key]))
            end
        else
            sub["gen"] = Dict{String, Vector{Any}}()
            sub["gencost"] = Dict{String, Vector{Any}}()
            sub["gencost_timeseries"] = Dict{String, Vector{Any}}()
        end

        add_case_index_maps!(sub)
        sub["bus_params"]             = get_bus_parameters(sub["bus"])
        sub["branch_params"]          = get_branch_parameters(sub["branch"], sub)
        sub["generator_params"]       = get_generator_parameters(sub["gen"])
        sub["power_balance_params"]   = get_power_balance_parameters(sub)
        sub["reactive_balance_params"]= get_reactive_power_parameters(sub)
        sub["objective_params"] = merge(
            get_objective_parameters(sub, sub["gen"], get(sub,"storage",Dict()),
                                     get(sub,"gencost",Dict()), get(sub,"prosumer",Dict()),
                                     get(sub,"gencost_timeseries",nothing)),
            (; sbase = sub["baseMVA"]),
        )

        sub["bus_set"] = Set(sub["bus"]["bus_i"])

        partitions[zone] = sub
    end
    bus_zone_appearance = Dict{Int, Set{Int}}()
    for (zone_id, sub) in partitions
        for b in sub["bus"]["bus_i"]
            push!(get!(bus_zone_appearance, b, Set{Int}()), zone_id)
        end
    end
    shared_buses = [b for (b, zs) in bus_zone_appearance if length(zs) > 1]
    for (_, sub) in partitions
        sub["shared_buses"] = [b for b in shared_buses if b in sub["bus_set"]]
    end
    return partitions, shared_buses
end

function _find_cut_edges(full::Dict{String, Any})
    zones = full["bus"]["zone"]
    bus_ids = full["bus"]["bus_i"]
    bus2zone = Dict(bus_ids[i] => zones[i] for i in eachindex(bus_ids))
    f_bus = full["branch"]["f_bus"]
    t_bus = full["branch"]["t_bus"] 
    cut = Int[]
    for i in eachindex(f_bus)
        if bus2zone[f_bus[i]] != bus2zone[t_bus[i]]
            push!(cut, i)
        end
    end
    return cut
end

function _zone_cut_edges(subdata::Dict{String, Any}, cut_edges::Vector{Int})
    local_keep = Set(subdata["branch_local_indices"])
    return [ℓ for ℓ in cut_edges if ℓ in local_keep]
end
