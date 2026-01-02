function map_items_by_bus(bus_vec::Vector{T}) where {T}
    mapping = Dict{Int, Vector{Int}}()
    for (i, b) in enumerate(bus_vec)
        push!(get!(mapping, b, Int[]), i)
    end
    return mapping
end

function add_case_index_maps!(case::Dict{String, Any})
    if haskey(case, "bus")
        bus_ids = case["bus"]["bus_i"]
        case["bus_idx_map"] = Dict(b => i for (i, b) in enumerate(bus_ids))
    end

    if haskey(case, "gen")
        gen_bus = case["gen"]["gen_bus"]
        case["gen_idx_map"] = Dict(i => b for (i, b) in enumerate(gen_bus)) # posição (gerador) -> barra
        case["bus_to_gens"] = map_items_by_bus(gen_bus)                     # barras -> posições (gerador)
    end

    if haskey(case, "branch")
        f_bus = case["branch"]["f_bus"] # posição linha from -> barra
        t_bus = case["branch"]["t_bus"] # posição linha to -> barra

        case["bus_to_branches_from"] = map_items_by_bus(f_bus) # barras -> linhas
        case["bus_to_branches_to"] = map_items_by_bus(t_bus) # barras -> linhas

        bus_to_neighbors = Dict{Int, Vector{Int}}()
        for (f, t) in zip(f_bus, t_bus)
            push!(get!(bus_to_neighbors, f, Int[]), t)
            push!(get!(bus_to_neighbors, t, Int[]), f) # pega barra vizinha na qual ela está conectada
        end
        case["bus_to_neighbors"] = bus_to_neighbors
    end

    if haskey(case, "prosumer")
        case["bus_to_prosumers"] = map_items_by_bus(case["prosumer"]["bus"])  # posição prosumer -> barra
    end

    if haskey(case, "storage")
        case["bus_to_batteries"] = map_items_by_bus(case["storage"]["storage_bus"]) # posição bateria -> barra
    end

    return case
end

function convert_to_pu!(case::Dict{String, Any})
    Sbase = case["baseMVA"]

    # Barras
    for field in ["pd", "qd", "gs", "bs"]
        case["bus"][field] ./= Sbase
    end

    # Geradores
    for field in ["pg", "qg", "pmax", "pmin", "qmax", "qmin"]
        case["gen"][field] ./= Sbase
    end

    # Ramos (limites térmicos)
    for field in ["rate_a", "rate_b", "rate_c"]
        case["branch"][field] ./= Sbase
    end

    # Armazenamento
    if haskey(case, "storage")
        for field in [
            "ps",
            "qs",
            "energy",
            "energy_rating",
            "charge_rating",
            "discharge_rating",
            "thermal_rating",
            "qmin",
            "qmax",
            "q_loss",
            "p_loss",
        ]
            case["storage"][field] ./= Sbase
        end
    end

    return case
end

function convert_timeseries_to_pu!(case::Dict{String, Any})
    Sbase = case["baseMVA"]

    # Timeseries de barra
    if haskey(case, "bus_timeseries")
        for field in ["pd", "qd", "gs", "bs"]
            if haskey(case["bus_timeseries"], field)
                case["bus_timeseries"][field] ./= Sbase
            end
        end
    end

    # Timeseries de prosumer
    if haskey(case, "prosumer_timeseries")
        for field in ["p_load", "p_pv", "p_curtail"]
            if haskey(case["prosumer_timeseries"], field)
                case["prosumer_timeseries"][field] ./= Sbase
            end
        end
    end

    # Timeseries de consumidor
    if haskey(case, "consumer_timeseries")
        for field in ["demand", "q_demand"]
            if haskey(case["consumer_timeseries"], field)
                case["consumer_timeseries"][field] ./= Sbase
            end
        end
    end

    return case
end

# ---------------------------------------------------------------------
# NOTE / TODO (Timeseries review needed)
#
# This whole timeseries preprocessing pipeline still needs a full review.
# ---------------------------------------------------------------------

function preprocess_case_data!(
    data::Dict{String, Any},
    ::Type{T} = OPFModel,
) where {T <: OPFModel}
    if is_centralized_model(T)
        all_time_steps = get_model_times(data)
        time_horizon = isempty(all_time_steps) ? 1 : maximum(all_time_steps)
        data["T"] = time_horizon
    end

    time_horizon = data["T"]

    check_time_alignment(data)

    ensure_bus_timeseries!(data, time_horizon)
    ensure_gencost_timeseries!(data, time_horizon)
    ensure_prosumer_timeseries!(data, time_horizon)
    ensure_consumer_timeseries!(data, time_horizon)

    extend_bus_timeseries!(data, time_horizon)
    extend_gencost_timeseries!(data, time_horizon)
    extend_prosumer_timeseries!(data, time_horizon)
    extend_consumer_timeseries!(data, time_horizon)

    sort_timeseries_by_time_and_id!(data)

    add_case_index_maps!(data)

    convert_to_pu!(data)
    convert_timeseries_to_pu!(data)

    return nothing
end

getfieldsafe(data::Dict, field::String, default = missing) = get(data, field, default)

function get_bus_index_map(bus_ids::AbstractVector)
    return Dict(b => i for (i, b) in enumerate(bus_ids))
end

# ============================================================
# ADMM preprocessing
# ============================================================

function _find_interzonal_branches(full::Dict{String, Any})
    zones = full["bus"]["zone"]
    bus_ids = full["bus"]["bus_i"]
    bus2zone = Dict(bus_ids[i] => zones[i] for i in eachindex(bus_ids))

    f_bus = full["branch"]["f_bus"]
    t_bus = full["branch"]["t_bus"]
    br_status =
        haskey(full["branch"], "br_status") ? full["branch"]["br_status"] :
        ones(Int, length(f_bus))

    cut = Int[]
    for i in eachindex(f_bus)
        if br_status[i] == 1 && bus2zone[f_bus[i]] != bus2zone[t_bus[i]]
            push!(cut, i)
        end
    end
    return cut
end

@inline function branch_in_service_mask(branch::AbstractDict{String, V}) where {V}
    f = branch["f_bus"]
    br_status = haskey(branch, "br_status") ? branch["br_status"] : ones(Int, length(f))
    @assert length(br_status) == length(f)
    return br_status .== 1
end

@inline _slice_table_by_mask(
    tbl::AbstractDict{String, V},
    mask::Vector{Int},
    ref_len::Int,
) where {V} =
    Dict{String, Any}(k => v[mask] for (k, v) in tbl if length(v) == ref_len)

@inline _empty_like(tbl::AbstractDict{String, V}) where {V} =
    Dict{String, Any}(k => Any[] for k in keys(tbl))

@inline function _find_extra_neighbor_buses!(
    extra::Set{Int},
    zone_bus_set::Set{Int},
    f_bus_full::AbstractVector{<:Integer},
    t_bus_full::AbstractVector{<:Integer},
    branch_status_mask::AbstractVector{Bool},
)
    @inbounds for i in eachindex(f_bus_full)
        branch_status_mask[i] || continue
        f_in = in(f_bus_full[i], zone_bus_set)
        t_in = in(t_bus_full[i], zone_bus_set)
        if f_in ⊻ t_in
            push!(extra, f_in ? Int(t_bus_full[i]) : Int(f_bus_full[i]))
        end
    end
    return extra
end

@inline function _build_sub_meta(data::AbstractDict{String, V}, zone::Int) where {V}
    sub = Dict{String, Any}()
    sub["baseMVA"] = data["baseMVA"]
    sub["T"] = get(data, "T", 1)
    sub["source_type"] = get(data, "source_type", "matpower")
    sub["source_version"] = get(data, "source_version", "")
    sub["name"] = "$(get(data, "name", "case"))_zone_$zone"
    return sub
end

@inline function _mark_virtual_buses!(
    sub_bus::AbstractDict{String, V},
    zone_bus_set::Set{Int},
) where {V}
    bus_i_any = sub_bus["bus_i"]
    bus_i = Int.(bus_i_any)
    is_virtual = [!in(b, zone_bus_set) for b in bus_i]
    sub_bus["is_virtual"] = is_virtual

    @inbounds for i in eachindex(bus_i)
        if is_virtual[i]
            sub_bus["pd"][i] = 0.0
            sub_bus["qd"][i] = 0.0
            sub_bus["bus_type"][i] = 1
        end
    end
    return nothing
end

@inline function _infer_bus_key(tbl::AbstractDict{String, V}) where {V}
    return haskey(tbl, "bus") ? "bus" :
           haskey(tbl, "gen_bus") ? "gen_bus" :
           haskey(tbl, "storage_bus") ? "storage_bus" :
           nothing
end

@inline function _filter_table_by_zone_bus(
    tbl::AbstractDict{String, V},
    bus_key::String,
    zone_bus_set::Set{Int},
) where {V}
    keep = findall(i -> in(Int(tbl[bus_key][i]), zone_bus_set), eachindex(tbl[bus_key]))
    if isempty(keep)
        return _empty_like(tbl)
    end
    return Dict{String, Any}(
        k => v[keep] for (k, v) in tbl if length(v) == length(tbl[bus_key])
    )
end

@inline function _build_branch_sub!(
    sub::Dict{String, Any},
    data::AbstractDict{String, V},
    zone_bus_set::Set{Int},
    f_bus_full::AbstractVector{<:Integer},
    t_bus_full::AbstractVector{<:Integer},
    branch_status_mask::AbstractVector{Bool},
) where {V}
    keep_all = findall(
        i ->
            branch_status_mask[i] &&
                (in(f_bus_full[i], zone_bus_set) || in(t_bus_full[i], zone_bus_set)),
        eachindex(f_bus_full),
    )

    branch_tbl = data["branch"]
    sub["branch"] = _slice_table_by_mask(branch_tbl, keep_all, length(f_bus_full))
    sub["branch_full"] = branch_tbl
    sub["branch_local_indices"] = keep_all

    branch_global_to_local = Dict{Int, Int}()
    @inbounds for (loc, glob) in enumerate(keep_all)
        branch_global_to_local[glob] = loc
    end
    sub["branch_global_to_local"] = branch_global_to_local
    return keep_all
end

@inline function _filter_timeseries!(
    sub::Dict{String, Any},
    data::AbstractDict{String, V},
    zone_bus_set::Set{Int},
) where {V}
    for key in ("bus_timeseries", "prosumer_timeseries", "consumer_timeseries")
        haskey(data, key) || continue
        ts = data[key]
        keep = findall(i -> in(Int(ts["bus"][i]), zone_bus_set), eachindex(ts["bus"]))
        sub[key] = Dict{String, Any}(
            k => v[keep] for (k, v) in ts if length(v) == length(ts["bus"])
        )
    end
    return nothing
end

@inline function _filter_gencosts!(
    sub::Dict{String, Any},
    data::AbstractDict{String, V},
    zone_bus_set::Set{Int},
) where {V}
    if !haskey(data, "gen")
        sub["gen"] = Dict{String, Vector{Any}}()
        sub["gencost"] = Dict{String, Vector{Any}}()
        sub["gencost_timeseries"] = Dict{String, Vector{Any}}()
        return nothing
    end

    gtbl = data["gen"]
    bkey = haskey(gtbl, "bus") ? "bus" : "gen_bus"

    if haskey(data, "gencost")
        keep = findall(i -> in(Int(gtbl[bkey][i]), zone_bus_set), eachindex(gtbl[bkey]))
        sub["gencost"] = Dict{String, Any}(
            k => v[keep] for (k, v) in data["gencost"] if length(v) == length(gtbl[bkey])
        )
    end

    if haskey(data, "gencost_timeseries")
        keep = findall(i -> in(Int(gtbl[bkey][i]), zone_bus_set), eachindex(gtbl[bkey]))
        sub["gencost_timeseries"] = Dict{String, Any}(
            k => v[keep] for
            (k, v) in data["gencost_timeseries"] if length(v) == length(gtbl[bkey])
        )
    end

    return nothing
end

# ============================================================
# Prepare ADMM partitions
# ============================================================

function prepare_admm_partitions(data::Dict{String, Any})
    bus = data["bus"]

    bus_ids = Int.(bus["bus_i"])
    zones = Int.(bus["zone"])
    unique_zones = unique(zones)

    partitions = Dict{Int, Dict{String, Any}}()

    branch = data["branch"]
    f_bus_full = Int.(branch["f_bus"])
    t_bus_full = Int.(branch["t_bus"])
    branch_status_mask = branch_in_service_mask(branch)

    for zone in unique_zones
        zone_mask = findall(==(zone), zones)

        zone_bus_ids = @views bus_ids[zone_mask]
        zone_bus_set = Set{Int}(zone_bus_ids)

        extra = Set{Int}()
        _find_extra_neighbor_buses!(
            extra,
            zone_bus_set,
            f_bus_full,
            t_bus_full,
            branch_status_mask,
        )

        all_bus_set = union(zone_bus_set, extra)
        bus_mask = findall(b -> in(b, all_bus_set), bus_ids)

        sub = _build_sub_meta(data, zone)

        sub["bus"] = _slice_table_by_mask(data["bus"], bus_mask, length(bus_ids))
        _mark_virtual_buses!(sub["bus"], zone_bus_set)

        for key in ("gen", "prosumer", "consumer", "storage")
            haskey(data, key) || continue
            tbl = data[key]
            bus_key = _infer_bus_key(tbl)
            sub[key] =
                bus_key === nothing ? _empty_like(tbl) :
                _filter_table_by_zone_bus(tbl, bus_key, zone_bus_set)
        end

        _build_branch_sub!(
            sub,
            data,
            zone_bus_set,
            f_bus_full,
            t_bus_full,
            branch_status_mask,
        )

        _filter_timeseries!(sub, data, zone_bus_set)
        _filter_gencosts!(sub, data, zone_bus_set)

        add_case_index_maps!(sub)

        sub["bus_set"] = Set(Int.(sub["bus"]["bus_i"]))
        partitions[zone] = sub
    end

    bus_zone_appearance = Dict{Int, Set{Int}}()
    for (zone_id, sub) in partitions, b in Int.(sub["bus"]["bus_i"])
        push!(get!(bus_zone_appearance, b, Set{Int}()), zone_id)
    end
    shared_buses = [b for (b, zs) in bus_zone_appearance if length(zs) > 1]
    shared_set = Set(shared_buses)

    for (_, sub) in partitions
        local_b = Int.(sub["bus"]["bus_i"])
        sub["shared_buses_local"] = [b for b in local_b if in(b, shared_set)]
    end

    return partitions, shared_buses
end
