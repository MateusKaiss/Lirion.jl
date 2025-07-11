module ADMM

# Acesso ao tipo
import ..DistributedOPF: DistributedOPFModel

const PowerOptimization = Main.PowerOptimization

function run_admm(config::DistributedOPFModel, data)
    partitions = prepare_admm_partitions(data)

    # for k = 1:config.max_iter
    #     for zone_id in keys(partitions)
    #         # Ex: PowerOptimization.build_model(...)
    #     end
    # end

    return partitions
end

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
                if haskey(table, "bus")
                    keep = findall(
                        i -> table["bus"][i] in zone_bus_ids,
                        1:length(table["bus"]),
                    )
                    if !isempty(keep)
                        subdata[key] = Dict(
                            k => v[keep] for
                            (k, v) in table if length(v) == length(table["bus"])
                        )
                    else
                        subdata[key] = Dict(k => Any[] for k in keys(table))  # cria estrutura vazia
                    end
                else
                    subdata[key] = Dict(k => Any[] for k in keys(table))  # segurança extra
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

        # Tratamento especial para gen, gencost e gencost_timeseries
        if haskey(data, "gen")
            gen = data["gen"]
            if haskey(gen, "bus")
                keep = findall(i -> gen["bus"][i] in zone_bus_ids, 1:length(gen["bus"]))
                if !isempty(keep)
                    subdata["gen"] = Dict(
                        k => v[keep] for (k, v) in gen if length(v) == length(gen["bus"])
                    )
                else
                    subdata["gen"] = Dict(k => Any[] for k in keys(gen))
                end

                if haskey(data, "gencost")
                    subdata["gencost"] = Dict(
                        k => v[keep] for
                        (k, v) in data["gencost"] if length(v) == length(keep)
                    )
                end

                if haskey(data, "gencost_timeseries")
                    subdata["gencost_timeseries"] = Dict(
                        k => v[keep] for (k, v) in data["gencost_timeseries"] if
                        length(v) == length(keep)
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

        PowerOptimization.preprocess_case_data!(subdata)
        partitions[zone] = subdata
    end

    return partitions
end

end
