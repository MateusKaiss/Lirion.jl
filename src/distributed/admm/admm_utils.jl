function prepare_admm_partitions(data::Dict{String, Any})
    bus_data = data["bus"]
    zones = bus_data["zone"]
    bus_ids = bus_data["bus_i"]

    unique_zones = unique(zones)
    partitions = Dict{Int, Dict{String, Any}}()
    zone_map = Dict(bus_ids[i] => zones[i] for i in eachindex(bus_ids))
    interconnection_limits = Dict{Tuple{Int, Int}, Float64}()

    for zone in unique_zones
        zone_mask = findall(==(zone), zones)
        zone_bus_ids = bus_ids[zone_mask]

        # 🔧 Novo: adicionar também buses conectadas via ramos interzonais
        extra_bus_ids = Int[]
        if haskey(data, "branch")
            f_bus = data["branch"]["f_bus"]
            t_bus = data["branch"]["t_bus"]
            for i in eachindex(f_bus)
                if (f_bus[i] in zone_bus_ids) && !(t_bus[i] in zone_bus_ids)
                    push!(extra_bus_ids, t_bus[i])
                elseif (t_bus[i] in zone_bus_ids) && !(f_bus[i] in zone_bus_ids)
                    push!(extra_bus_ids, f_bus[i])
                end
            end
        end
        all_bus_ids = union(zone_bus_ids, extra_bus_ids)
        bus_mask = findall(b -> b in all_bus_ids, bus_ids)

        subdata = Dict{String, Any}()
        subdata["baseMVA"] = data["baseMVA"]
        subdata["T"] = data["T"]
        subdata["source_type"] = data["source_type"]
        subdata["source_version"] = data["source_version"]
        subdata["name"] = "$(data["name"])_zone_$zone"
        subdata["bus"] =
            Dict(k => v[bus_mask] for (k, v) in data["bus"] if length(v) == length(bus_ids))

        for key in ["gen", "prosumer", "consumer", "storage"]
            if haskey(data, key)
                table = data[key]
                bus_key =
                    haskey(table, "bus") ? "bus" :
                    haskey(table, "gen_bus") ? "gen_bus" :
                    haskey(table, "storage_bus") ? "storage_bus" : nothing
                if bus_key !== nothing
                    keep = findall(
                        i -> table[bus_key][i] in zone_bus_ids,
                        eachindex(table[bus_key]),
                    )
                    subdata[key] =
                        !isempty(keep) ?
                        Dict(
                            k => v[keep] for
                            (k, v) in table if length(v) == length(table[bus_key])
                        ) :
                        Dict(k => Any[] for k in keys(table))
                else
                    subdata[key] = Dict(k => Any[] for k in keys(table))
                end
            end
        end

        if haskey(data, "branch")
            f_bus = data["branch"]["f_bus"]
            t_bus = data["branch"]["t_bus"]
            rateA = data["branch"]["rate_a"]
            keep = findall(
                i -> f_bus[i] in zone_bus_ids || t_bus[i] in zone_bus_ids,
                eachindex(f_bus),
            )

            inter_out = Int[]
            inter_in = Int[]

            for i in keep
                zb_f = zone_map[f_bus[i]]
                zb_t = zone_map[t_bus[i]]
                if zb_f == zone && zb_t != zone
                    push!(inter_out, i)
                    interconnection_limits[(zb_f, zb_t)] =
                        get(interconnection_limits, (zb_f, zb_t), 0.0) + rateA[i]
                elseif zb_f != zone && zb_t == zone
                    push!(inter_in, i)
                    interconnection_limits[(zb_f, zb_t)] =
                        get(interconnection_limits, (zb_f, zb_t), 0.0) + rateA[i]
                end
            end

            subdata["inter_out_branches"] = inter_out
            subdata["inter_in_branches"] = inter_in
            bus_to_inter_out_branches = Dict{Int, Vector{Int}}()
            for l in inter_out
                f = data["branch"]["f_bus"][l]
                push!(get!(bus_to_inter_out_branches, f, Int[]), l)
            end
            bus_to_inter_in_branches = Dict{Int, Vector{Int}}()
            for l in inter_in
                t = data["branch"]["t_bus"][l]
                push!(get!(bus_to_inter_in_branches, t, Int[]), l)
            end

            subdata["bus_to_inter_out_branches"] = bus_to_inter_out_branches

            subdata["bus_to_inter_in_branches"] = bus_to_inter_in_branches

            subdata["branch"] = Dict(
                k => v[keep] for (k, v) in data["branch"] if length(v) == length(f_bus)
            )
        end

        for key in ["bus_timeseries", "prosumer_timeseries", "consumer_timeseries"]
            if haskey(data, key)
                ts = data[key]
                id_field = "bus"
                keep =
                    findall(i -> ts[id_field][i] in zone_bus_ids, eachindex(ts[id_field]))
                subdata[key] =
                    Dict(k => v[keep] for (k, v) in ts if length(v) == length(ts[id_field]))
            end
        end

        if haskey(data, "gen")
            gen = data["gen"]
            bus_field = haskey(gen, "bus") ? "bus" : "gen_bus"
            keep =
                findall(i -> gen[bus_field][i] in zone_bus_ids, eachindex(gen[bus_field]))
            subdata["gen"] =
                !isempty(keep) ?
                Dict(
                    k => v[keep] for (k, v) in gen if length(v) == length(gen[bus_field])
                ) :
                Dict(k => Any[] for k in keys(gen))

            if haskey(data, "gencost")
                subdata["gencost"] = Dict(
                    k => v[keep] for
                    (k, v) in data["gencost"] if length(v) == length(gen[bus_field])
                )
            end
            if haskey(data, "gencost_timeseries")
                subdata["gencost_timeseries"] = Dict(
                    k => v[keep] for (k, v) in data["gencost_timeseries"] if
                    length(v) == length(gen[bus_field])
                )
            end
        else
            subdata["gen"] = Dict{String, Vector{Any}}()
            subdata["gencost"] = Dict{String, Vector{Any}}()
            subdata["gencost_timeseries"] = Dict{String, Vector{Any}}()
        end

        add_case_index_maps!(subdata)

        subdata["bus_params"] = get_bus_parameters(subdata["bus"])
        subdata["branch_params"] = get_branch_parameters(subdata["branch"], subdata)
        subdata["generator_params"] = get_generator_parameters(subdata["gen"])
        subdata["power_balance_params"] = get_power_balance_parameters(subdata)
        subdata["reactive_balance_params"] = get_reactive_power_parameters(subdata)
        subdata["objective_params"] = merge(
            get_objective_parameters(
                subdata,
                subdata["gen"],
                get(subdata, "storage", Dict()),
                subdata["gencost"],
                get(subdata, "prosumer", Dict()),
                get(subdata, "gencost_timeseries", nothing),
            ),
            (; sbase = subdata["baseMVA"]),
        )

        partitions[zone] = subdata
    end

    return partitions, interconnection_limits
end

function find_shared_buses(
    data::Dict{String, Any},
    partitions::Dict{Int, Dict{String, Any}},
)
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

function compute_residuals(
    θ_avg::Dict{Int, Float64},
    θ_avg_prev::Dict{Int, Float64},
    partitions,
    shared_buses,
)
    primal = 0.0
    dual = 0.0

    for b in shared_buses
        for (_, pdata) in partitions
            if b in pdata["bus"]["bus_i"]
                model = pdata["model"]
                θ_val = value(model[:theta][b, 1])
                primal += (θ_val - θ_avg[b])^2
            end
        end
        if haskey(θ_avg_prev, b)
            dual += (θ_avg[b] - θ_avg_prev[b])^2
        end
    end

    return sqrt(primal), sqrt(dual)
end
