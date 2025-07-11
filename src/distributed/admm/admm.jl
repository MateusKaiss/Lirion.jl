module ADMM

# Acesso ao tipo
import ..DistributedOPF: DistributedOPFModel

function run_admm(config::DistributedOPFModel, data)

    # partitions = partition_system!(data)

    # for k = 1:config.max_iter
    #     for zone_id in keys(partitions)
    #         # Ex: PowerOptimization.build_model(...)
    #     end
    # end

    return partitions
end

# function partition_system!(data::Dict{String, Any})
#     bus_data = data["bus"]
#     zones = bus_data["zone"]
#     bus_ids = bus_data["bus_i"]

#     unique_zones = unique(zones)
#     partitions = Dict{Int, Dict{String, Any}}()

#     for zone in unique_zones
#         subdata = Dict{String, Any}()
#         zone_mask = findall(==(zone), zones)
#         zone_bus_ids = bus_ids[zone_mask]

#         sub_bus = Dict{String, Vector}()
#         for (k, v) in bus_data
#             if length(v) == length(zones)
#                 sub_bus[k] = v[zone_mask]
#             else
#                 sub_bus[k] = v 
#             end
#         end
#         subdata["bus"] = sub_bus

#         for key in ["gen", "prosumer", "consumer", "storage"]
#             if haskey(data, key)
#                 table = data[key]
#                 if !haskey(table, "bus")
#                     @warn "Tabela $key não possui coluna 'bus'. Ignorando particionamento desta tabela."
#                     continue
#                 end

#                 bus_col = table["bus"]
#                 keep = findall(x -> x in zone_bus_ids, bus_col)

#                 sub_table = Dict{String, Vector}()
#                 for (k, v) in table
#                     if length(v) == length(bus_col)
#                         sub_table[k] = v[keep]
#                     else
#                         sub_table[k] = v
#                     end
#                 end
#                 subdata[key] = sub_table
#             end
#         end

#         if haskey(data, "branch")
#             branch = data["branch"]
#             f_bus = branch["f_bus"]
#             t_bus = branch["t_bus"]
#             keep = findall(i -> (f_bus[i] in zone_bus_ids) && (t_bus[i] in zone_bus_ids), eachindex(f_bus))

#             sub_branch = Dict{String, Vector}()
#             for (k, v) in branch
#                 if length(v) == length(f_bus)
#                     sub_branch[k] = v[keep]
#                 else
#                     sub_branch[k] = v
#                 end
#             end
#             subdata["branch"] = sub_branch
#         end

#         partitions[zone] = subdata
#     end
#     return partitions
# end

end
