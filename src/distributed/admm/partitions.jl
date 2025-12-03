# function prepare_admm_partitions(data::Dict{String, Any})
#     bus = data["bus"]
#     bus_ids = bus["bus_i"]
#     zones = bus["zone"]
#     unique_zones = unique(zones)

#     partitions = Dict{Int, Dict{String, Any}}()
#     zone_map = Dict(bus_ids[i] => zones[i] for i in eachindex(bus_ids))

#     has_branch = haskey(data, "branch")
#     if has_branch
#         f_bus_full = data["branch"]["f_bus"]
#         t_bus_full = data["branch"]["t_bus"]
#         br_status =
#             haskey(data["branch"], "br_status") ?
#             data["branch"]["br_status"] :
#             ones(Int, length(f_bus_full))
#         @assert length(br_status) == length(f_bus_full)
#         active_mask_full = [br_status[i] == 1 for i in eachindex(f_bus_full)]
#     end

#     for zone in unique_zones
#         zone_mask = findall(==(zone), zones)
#         zone_bus_ids = @views bus_ids[zone_mask]
#         zone_bus_set = Set(zone_bus_ids)               # O(1) membership

#         extra = Set{Int}()
#         if has_branch
#             @inbounds for i in eachindex(f_bus_full)
#                 f_in = in(f_bus_full[i], zone_bus_set)
#                 t_in = in(t_bus_full[i], zone_bus_set)
#                 if f_in ⊻ t_in
#                     push!(extra, f_in ? t_bus_full[i] : f_bus_full[i])
#                 end
#             end
#         end
#         extra_bus_ids = collect(extra)
#         all_bus_set = union(zone_bus_set, extra)       # ambos Set
#         bus_mask = findall(b -> in(b, all_bus_set), bus_ids)

#         sub = Dict{String, Any}()
#         sub["baseMVA"] = data["baseMVA"]
#         sub["T"] = get(data, "T", 1)
#         sub["source_type"] = get(data, "source_type", "matpower")
#         sub["source_version"] = get(data, "source_version", "")
#         sub["name"] = "$(get(data,"name","case"))_zone_$zone"
#         sub["bus"] =
#             Dict(k => v[bus_mask] for (k, v) in data["bus"] if length(v) == length(bus_ids))

#         # flag virtual + poda de cargas em virtuais
#         is_virtual = [!in(b, zone_bus_set) for b in sub["bus"]["bus_i"]]
#         sub["bus"]["is_virtual"] = is_virtual
#         @inbounds for (i, b) in enumerate(sub["bus"]["bus_i"])
#             if is_virtual[i]
#                 sub["bus"]["pd"][i] = 0.0
#                 sub["bus"]["qd"][i] = 0.0
#                 sub["bus"]["bus_type"][i] = 1
#             end
#         end

#         # tabelas ancoradas em "bus"/"gen_bus"/"storage_bus": use zone_bus_set
#         for key in ("gen", "prosumer", "consumer", "storage")
#             if haskey(data, key)
#                 tbl = data[key]
#                 bus_key =
#                     haskey(tbl, "bus") ? "bus" :
#                     haskey(tbl, "gen_bus") ? "gen_bus" :
#                     haskey(tbl, "storage_bus") ? "storage_bus" : nothing
#                 if bus_key !== nothing
#                     keep = findall(
#                         i -> in(tbl[bus_key][i], zone_bus_set),
#                         eachindex(tbl[bus_key]),
#                     )
#                     sub[key] =
#                         !isempty(keep) ?
#                         Dict(
#                             k => v[keep] for
#                             (k, v) in tbl if length(v) == length(tbl[bus_key])
#                         ) :
#                         Dict(k => Any[] for k in keys(tbl))
#                 else
#                     sub[key] = Dict(k => Any[] for k in keys(tbl))
#                 end
#             end
#         end

#         if has_branch
#             keep_all = findall(
#                 i ->
#                     (active_mask_full[i]) &&
#                         (
#                             in(f_bus_full[i], zone_bus_set) ||
#                             in(t_bus_full[i], zone_bus_set)
#                         ),
#                 eachindex(f_bus_full),
#             )
#             sub["branch"] = Dict(
#                 k => v[keep_all] for
#                 (k, v) in data["branch"] if length(v) == length(f_bus_full)
#             )
#             sub["branch_full"] = data["branch"]
#             sub["branch_local_indices"] = keep_all                # globais dos ativos mantidos
#             # mapa global->local (só para os ativos mantidos)
#             edge_to_local = Dict{Int, Int}()
#             @inbounds for (loc, glob) in enumerate(keep_all)
#                 edge_to_local[glob] = loc
#             end
#             sub["edge_to_local"] = edge_to_local
#         end

#         # timeseries por bus
#         for key in ("bus_timeseries", "prosumer_timeseries", "consumer_timeseries")
#             if haskey(data, key)
#                 ts = data[key]
#                 keep = findall(i -> in(ts["bus"][i], zone_bus_set), eachindex(ts["bus"]))
#                 sub[key] =
#                     Dict(k => v[keep] for (k, v) in ts if length(v) == length(ts["bus"]))
#             end
#         end

#         # gencost / gencost_timeseries (use zone_bus_set)
#         if haskey(data, "gen")
#             if haskey(data, "gencost")
#                 gtbl = data["gen"]
#                 bkey = haskey(gtbl, "bus") ? "bus" : "gen_bus"
#                 keep = findall(i -> in(gtbl[bkey][i], zone_bus_set), eachindex(gtbl[bkey]))
#                 sub["gencost"] = Dict(
#                     k => v[keep] for
#                     (k, v) in data["gencost"] if length(v) == length(gtbl[bkey])
#                 )
#             end
#             if haskey(data, "gencost_timeseries")
#                 gtbl = data["gen"]
#                 bkey = haskey(gtbl, "bus") ? "bus" : "gen_bus"
#                 keep = findall(i -> in(gtbl[bkey][i], zone_bus_set), eachindex(gtbl[bkey]))
#                 sub["gencost_timeseries"] = Dict(
#                     k => v[keep] for (k, v) in data["gencost_timeseries"] if
#                     length(v) == length(gtbl[bkey])
#                 )
#             end
#         else
#             sub["gen"] = Dict{String, Vector{Any}}()
#             sub["gencost"] = Dict{String, Vector{Any}}()
#             sub["gencost_timeseries"] = Dict{String, Vector{Any}}()
#         end

#         add_case_index_maps!(sub)
#         sub["bus_params"] = get_bus_parameters(sub["bus"])
#         sub["branch_params"] = get_branch_parameters(sub["branch"], sub)
#         sub["generator_params"] = get_generator_parameters(sub["gen"])
#         sub["power_balance_params"] = get_power_balance_parameters(sub)
#         sub["reactive_balance_params"] = get_reactive_power_parameters(sub)
#         sub["objective_params"] = merge(
#             get_objective_parameters(
#                 sub,
#                 sub["gen"],
#                 get(sub, "storage", Dict()),
#                 get(sub, "gencost", Dict()),
#                 get(sub, "prosumer", Dict()),
#                 get(sub, "gencost_timeseries", nothing),
#             ),
#             (; sbase = sub["baseMVA"]),
#         )

#         sub["bus_set"] = Set(sub["bus"]["bus_i"])

#         partitions[zone] = sub
#     end

#     # mapeie compartilhadas uma vez
#     bus_zone_appearance = Dict{Int, Set{Int}}()
#     for (zone_id, sub) in partitions, b in sub["bus"]["bus_i"]
#         push!(get!(bus_zone_appearance, b, Set{Int}()), zone_id)
#     end
#     shared_buses = [b for (b, zs) in bus_zone_appearance if length(zs) > 1]
#     shared_set = Set(shared_buses)
#     for (_, sub) in partitions
#         local_b = sub["bus"]["bus_i"]
#         sub["shared_buses_local"] = [b for b in local_b if in(b, shared_set)]
#     end
#     return partitions, shared_buses
# end

function _set_model_silent!(m::JuMP.Model)
    _safe(k, v) = (
        try
            set_optimizer_attribute(m, k, v)
        catch
        end
    )

    try
        set_silent(m)
    catch
    end
    try
        set_optimizer_attribute(m, JuMP.MOI.Silent(), true)
    catch
    end

    solver = ""
    try
        backend = JuMP.backend(m)
        solver = lowercase(JuMP.MOI.get(backend, JuMP.MOI.SolverName()))
    catch
    end

    if occursin("ipopt", solver)
        _safe("print_level", 0)
        _safe("sb", "yes")
    elseif occursin("gurobi", solver)
        _safe("OutputFlag", 0)
    elseif occursin("cplex", solver)
        _safe("CPXPARAM_ScreenOutput", 0)
    elseif occursin("highs", solver)
        _safe("log_to_console", false)
        _safe("output_flag", false)
    elseif occursin("scs", solver)
        _safe("verbose", false)
    elseif occursin("clarabel", solver)
        _safe("verbose", false)
    elseif occursin("glpk", solver)
        _safe("msg_lev", 0)
    else
        _safe("verbose", false)
        _safe("output_flag", 0)
        _safe("print_level", 0)
    end
    return nothing
end

# function _find_cut_edges(full::Dict{String, Any})
#     zones = full["bus"]["zone"]
#     bus_ids = full["bus"]["bus_i"]
#     bus2zone = Dict(bus_ids[i] => zones[i] for i in eachindex(bus_ids))
#     f_bus = full["branch"]["f_bus"]
#     t_bus = full["branch"]["t_bus"]
#     br_status =
#         haskey(full["branch"], "br_status") ? full["branch"]["br_status"] :
#         ones(Int, length(f_bus))
#     cut = Int[]
#     for i in eachindex(f_bus)
#         if br_status[i] == 1 && bus2zone[f_bus[i]] != bus2zone[t_bus[i]]
#             push!(cut, i)
#         end
#     end
#     return cut
# end

function _zone_cut_edges(subdata::Dict{String, Any}, cut_edges::Vector{Int})
    local_keep = Set(get(subdata, "branch_local_indices", Int[]))
    return [ℓ for ℓ in cut_edges if ℓ in local_keep]
end

# ---------- Helpers ----------
@inline function _edge_local(sub::Dict{String, Any}, e::Int)
    edgemap = get(sub, "edge_to_local", nothing)
    return edgemap === nothing ? findfirst(==(e), sub["branch_local_indices"]) :
           get(edgemap, e, nothing)
end

@inline function _bus_set(sub::Dict{String, Any})
    return get(sub, "bus_set", Set(Int[b for b in sub["bus"]["bus_i"]]))
end
