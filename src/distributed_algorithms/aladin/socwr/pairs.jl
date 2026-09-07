# ============================================================
# ALADIN-SOCWR — Consensus pair types and builders
# ============================================================

# ============================================================
# Pair types (consensus structures)
# ============================================================

# ------------------------------------------------------------
# Shared bus pair
# ------------------------------------------------------------

struct SharedBusPairSOCWR
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    v2_col_i::Int
    zone_j::Int
    v2_col_j::Int
end

# ------------------------------------------------------------
# Interzonal line pair (TieLine sectioning)
# ------------------------------------------------------------

struct InterzonalLinePairSOCWR
    consensus_idx::Int
    line_global::Int
    zone_i::Int
    line_local_i::Int
    zone_j::Int
    line_local_j::Int
end

# ------------------------------------------------------------
# Split bus pair (SplitBus sectioning)
# ------------------------------------------------------------

struct SplitBusPairSOCWR
    consensus_idx::Int
    bus_id::Int
    zones::Vector{Int}
    p_inj_cols::Vector{Int}
    q_inj_cols::Vector{Int}
end

# ------------------------------------------------------------
# Dummy line pair (DummyBus sectioning)
# ------------------------------------------------------------

struct DummyLinePairSOCWR
    consensus_idx::Int
    zone_i::Int
    local_i::Int
    cols_i::Vector{Int}
    pcoef_i::Vector{Float64}
    qcoef_i::Vector{Float64}
    pconst_i::Float64
    qconst_i::Float64
    zone_j::Int
    local_j::Int
    cols_j::Vector{Int}
    pcoef_j::Vector{Float64}
    qcoef_j::Vector{Float64}
    pconst_j::Float64
    qconst_j::Float64
end

# ============================================================
# Pair builders
# ============================================================

# ------------------------------------------------------------
# Shared bus pairs
# ------------------------------------------------------------

function _build_shared_bus_pairs(
    partitions::Partitions{PartitionSOCWR},
    zones::Vector{AladinZoneDataSOCWR},
)::Vector{SharedBusPairSOCWR}
    visits = Dict{Int, Vector{Tuple{Int, Int}}}()    # cidx → [(zone_pos, v2_col)]
    bus_id_by_cidx = Dict{Int, Int}()                # cidx → physical bus id

    for (zone_pos, part) in enumerate(partitions.parts)
        zone_data = zones[zone_pos]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            consensus_idx = part.shared_bus_consensus_idx[i]
            consensus_idx == 0 && continue
            v2_col = zone_data.shared_v2_cols[i]
            cidx_visits = get!(visits, consensus_idx, Tuple{Int, Int}[])
            push!(cidx_visits, (zone_pos, v2_col))
            bus_id_by_cidx[consensus_idx] = part.zone_shared_bus_ids[i]
        end
    end

    pairs = SharedBusPairSOCWR[]
    for consensus_idx in sort!(collect(keys(visits)))
        cidx_visits = visits[consensus_idx]
        if length(cidx_visits) < 2
            @warn "ALADIN-SOCWR pair formation: shared bus $(get(bus_id_by_cidx, consensus_idx, "?")) " *
                  "(consensus_idx $(consensus_idx)) is present in only " *
                  "$(length(cidx_visits)) zone(s) — coupling DROPPED. Check the partition."
            continue
        end
        sort!(cidx_visits; by = first)
        root_zone, root_col = cidx_visits[1]
        for k in 2:length(cidx_visits)
            peer_zone, peer_col = cidx_visits[k]
            push!(
                pairs,
                SharedBusPairSOCWR(
                    consensus_idx,
                    bus_id_by_cidx[consensus_idx],
                    root_zone, root_col,
                    peer_zone, peer_col,
                ),
            )
        end
    end
    return pairs
end

# ------------------------------------------------------------
# Interzonal line pairs
# ------------------------------------------------------------

function _build_interzonal_line_pairs(
    partitions::Partitions{PartitionSOCWR},
    interzonal::Vector{Int},
)::Vector{InterzonalLinePairSOCWR}
    pairs = InterzonalLinePairSOCWR[]
    for (consensus_idx, line_global) in enumerate(interzonal)
        line_visits = Tuple{Int, Int}[]
        for (zone_pos, part) in enumerate(partitions.parts)
            part.branch_idx_full_to_zone[line_global] == 0 && continue
            local_idx = findfirst(==(line_global), part.zone_edges)
            local_idx === nothing && continue
            push!(line_visits, (zone_pos, local_idx))
        end
        if length(line_visits) != 2
            # an interzonal line must live in exactly 2 zones.
            @warn "ALADIN-SOCWR pair formation: interzonal line $(line_global) " *
                  "(consensus_idx $(consensus_idx)) found in $(length(line_visits)) " *
                  "zone(s), expected 2 — coupling DROPPED. Check the partition."
            continue
        end
        zone_i, local_i = line_visits[1]
        zone_j, local_j = line_visits[2]
        push!(
            pairs,
            InterzonalLinePairSOCWR(
                consensus_idx, line_global,
                zone_i, local_i,
                zone_j, local_j,
            ),
        )
    end
    return pairs
end

# ------------------------------------------------------------
# Split bus pairs
# ------------------------------------------------------------

function _build_split_bus_pairs(
    zones::Vector{AladinZoneDataSOCWR},
    bus_groups::Vector{Tuple{Int, Vector{Int}}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)::Vector{SplitBusPairSOCWR}
    pairs = SplitBusPairSOCWR[]
    isempty(bus_groups) && return pairs
    sizehint!(pairs, length(bus_groups))

    for (consensus_idx, (bus_id, group_zones)) in enumerate(bus_groups)
        copy_zones = Int[]
        p_cols = Int[]
        q_cols = Int[]
        ok = true
        for zone in group_zones
            p_inj_var = get(p_inj_vars[zone], bus_id, nothing)
            q_inj_var = get(q_inj_vars[zone], bus_id, nothing)
            if p_inj_var === nothing || q_inj_var === nothing
                ok = false
                break
            end
            push!(copy_zones, zone)
            push!(p_cols, zones[zone].var_to_col[p_inj_var])
            push!(q_cols, zones[zone].var_to_col[q_inj_var])
        end
        if !(ok && length(copy_zones) >= 2)
            @warn "ALADIN-SOCWR pair formation: split bus $(bus_id) (group $(consensus_idx), " *
                  "zones $(group_zones)) has $(length(copy_zones)) usable injection " *
                  "copies$(ok ? "" : " (a zone lacks its p_inj/q_inj variable)") — " *
                  "coupling DROPPED. Check the partition."
            continue
        end
        push!(pairs, SplitBusPairSOCWR(consensus_idx, bus_id, copy_zones, p_cols, q_cols))
    end
    return pairs
end

# ------------------------------------------------------------
# Dummy line pairs
# ------------------------------------------------------------

function _build_dummy_line_pairs(
    partitions::Partitions{PartitionSOCWR},
    zones::Vector{AladinZoneDataSOCWR},
    sys_branch_pairs::Vector{NTuple{4, Int}},
)::Vector{DummyLinePairSOCWR}
    pairs = DummyLinePairSOCWR[]
    isempty(sys_branch_pairs) && return pairs
    sizehint!(pairs, length(sys_branch_pairs))

    for (consensus_idx, (line_global_a, zone_a, line_global_b, zone_b)) in
        enumerate(sys_branch_pairs)
        part_a = partitions.parts[zone_a]
        part_b = partitions.parts[zone_b]
        local_a = part_a.branch_idx_full_to_zone[line_global_a]
        local_b = part_b.branch_idx_full_to_zone[line_global_b]
        if local_a == 0 || local_b == 0
            @warn "ALADIN-SOCWR pair formation: dummy line pair $(consensus_idx) " *
                  "(lines $(line_global_a) in zone $(zone_a), $(line_global_b) in zone $(zone_b)) " *
                  "has a half missing from its zone (local ids $(local_a), $(local_b)) — " *
                  "coupling DROPPED. Check the partition."
            continue
        end

        zone_data_a = zones[zone_a]
        zone_data_b = zones[zone_b]
        model_a = zone_data_a.local_model
        model_b = zone_data_b.local_model

        prev_a = -model_a[:p_reverse_flow_branches][local_a, 1]
        qrev_a = -model_a[:q_reverse_flow_branches][local_a, 1]
        pflow_b = model_b[:p_flow_branches][local_b, 1]
        qflow_b = model_b[:q_flow_branches][local_b, 1]

        cols_i, pcoef_i, qcoef_i, pconst_i, qconst_i =
            _socwr_affine_pq_expansion(prev_a, qrev_a, zone_data_a.var_to_col)
        cols_j, pcoef_j, qcoef_j, pconst_j, qconst_j =
            _socwr_affine_pq_expansion(pflow_b, qflow_b, zone_data_b.var_to_col)

        push!(
            pairs,
            DummyLinePairSOCWR(
                consensus_idx,
                zone_a, local_a,
                cols_i, pcoef_i, qcoef_i, pconst_i, qconst_i,
                zone_b, local_b,
                cols_j, pcoef_j, qcoef_j, pconst_j, qconst_j,
            ),
        )
    end
    return pairs
end
