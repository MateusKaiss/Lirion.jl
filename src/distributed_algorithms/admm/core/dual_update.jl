# ============================================================
# Helpers
# ============================================================

# ------------------------------------------------------------
# Bus
# ------------------------------------------------------------

function _update_bus_duals_single!(
    duals::Matrix{Float64},
    partitions::Partitions,
    consensus::Vector{Float64},
    rho::Float64,
    snapshot_fn,
)
    for part in partitions
        x_local = snapshot_fn(part)
        p = part.position
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = part.shared_bus_consensus_idx[i]
            duals[j, p] += rho * (x_local[i] - consensus[j])
        end
    end
    return nothing
end

function _update_bus_duals_ac!(
    duals_theta::Matrix{Float64},
    duals_v_mag::Matrix{Float64},
    partitions::Partitions{PartitionAC},
    theta_consensus::Vector{Float64},
    v_mag_consensus::Vector{Float64},
    rho_theta::Float64,
    rho_v_mag::Float64,
)
    for part in partitions
        theta_loc = snapshot_theta_shared!(part)
        v_mag_loc = snapshot_v_mag_shared!(part)
        p = part.position
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = part.shared_bus_consensus_idx[i]
            duals_theta[j, p] +=
                rho_theta * _wrap_angle_diff(theta_loc[i] - theta_consensus[j])
            duals_v_mag[j, p] += rho_v_mag * (v_mag_loc[i] - v_mag_consensus[j])
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Branch
# ------------------------------------------------------------

function _update_branch_duals_single!(
    duals::AbstractMatrix{Float64},
    partitions::Partitions,
    consensus::AbstractVector{Float64},
    rho::Float64,
    interzonal::Vector{Int},
    extract_fn,
)
    for part in partitions
        p = part.position
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            val = extract_fn(part, l)
            duals[k, p] += rho * (val - consensus[k])
        end
    end
    return nothing
end

function _update_branch_duals_paired!(
    duals_a::AbstractMatrix{Float64},
    duals_b::AbstractMatrix{Float64},
    partitions::Partitions,
    consensus_a::AbstractVector{Float64},
    consensus_b::AbstractVector{Float64},
    rho_a::Float64,
    rho_b::Float64,
    interzonal::Vector{Int},
    extract_fn,
)
    for part in partitions
        p = part.position
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            a, b = extract_fn(part, l)
            duals_a[k, p] += rho_a * (a - consensus_a[k])
            duals_b[k, p] += rho_b * (b - consensus_b[k])
        end
    end
    return nothing
end

function _update_branch_duals_ac!(
    duals_p::AbstractMatrix{Float64},
    duals_q::AbstractMatrix{Float64},
    partitions::Partitions{PartitionAC},
    p_consensus::AbstractVector{Float64},
    q_consensus::AbstractVector{Float64},
    rho_p::Float64,
    rho_q::Float64,
    interzonal::Vector{Int},
)
    for part in partitions
        p = part.position
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            p_expr, q_expr = _pq_flow_from_exprs(part, l)
            p_val = JuMP.value(p_expr)
            q_val = JuMP.value(q_expr)
            isfinite(p_val) &&
                (duals_p[k, p] += rho_p * (p_val - p_consensus[k]))
            isfinite(q_val) &&
                (duals_q[k, p] += rho_q * (q_val - q_consensus[k]))
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Branch pairs
# ------------------------------------------------------------

function _update_branch_pairs_duals_single!(
    duals::AbstractMatrix{Float64},
    partitions::Partitions,
    consensus::AbstractVector{Float64},
    rho::Float64,
    pairs::Vector{NTuple{4, Int}},
    extract_fn,
)
    isempty(pairs) && return nothing
    @inbounds for (k, (ga, za, gb, zb)) in enumerate(pairs)
        pa = partitions.parts[za]
        la = pa.branch_idx_full_to_zone[ga]
        if la != 0
            vala = extract_fn(pa, la)
            duals[k, za] += rho * (vala - consensus[k])
        end
        pb = partitions.parts[zb]
        lb = pb.branch_idx_full_to_zone[gb]
        if lb != 0
            valb = extract_fn(pb, lb)
            duals[k, zb] += rho * (valb - consensus[k])
        end
    end
    return nothing
end

function _update_branch_pairs_duals_ac!(
    duals_p::AbstractMatrix{Float64},
    duals_q::AbstractMatrix{Float64},
    partitions::Partitions{PartitionAC},
    p_consensus::AbstractVector{Float64},
    q_consensus::AbstractVector{Float64},
    rho_p::Float64,
    rho_q::Float64,
    pairs::Vector{NTuple{4, Int}},
)
    isempty(pairs) && return nothing
    @inbounds for (k, (ga, za, gb, zb)) in enumerate(pairs)
        pa = partitions.parts[za]
        pb = partitions.parts[zb]
        la = pa.branch_idx_full_to_zone[ga]
        lb = pb.branch_idx_full_to_zone[gb]
        if la != 0
            p_expr, q_expr = _pq_flow_from_exprs(pa, la)
            p_val = JuMP.value(p_expr)
            q_val = JuMP.value(q_expr)
            isfinite(p_val) &&
                (duals_p[k, za] += rho_p * (p_val - p_consensus[k]))
            isfinite(q_val) &&
                (duals_q[k, za] += rho_q * (q_val - q_consensus[k]))
        end
        if lb != 0
            p_expr, q_expr = _pq_flow_from_exprs(pb, lb)
            p_val = JuMP.value(p_expr)
            q_val = JuMP.value(q_expr)
            isfinite(p_val) &&
                (duals_p[k, zb] += rho_p * (p_val - p_consensus[k]))
            isfinite(q_val) &&
                (duals_q[k, zb] += rho_q * (q_val - q_consensus[k]))
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Bus pairs
# ------------------------------------------------------------

function _update_bus_pairs_duals_signed!(
    duals::AbstractMatrix{Float64},
    partitions::Partitions,
    consensus::AbstractVector{Float64},
    rho::Float64,
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)
    isempty(bus_pairs) && return nothing
    @inbounds for (k, (bid, za, zb)) in enumerate(bus_pairs)
        z = consensus[k]
        va = _p_inj_value(partitions.parts[za], p_inj_vars[za], bid)
        vb = _p_inj_value(partitions.parts[zb], p_inj_vars[zb], bid)
        duals[k, za] += rho * (va - z)
        duals[k, zb] += rho * (vb + z)
    end
    return nothing
end

# ============================================================
# DC
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function dual_update!(ws::DCStateWorkspace)
    _update_bus_duals_single!(
        ws.vars.theta.duals,
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.theta.rho[],
        snapshot_theta_shared!,
    )
    return nothing
end

# ------------------------------------------------------------
# Power — TieLine
# ------------------------------------------------------------

function dual_update!(ws::DCPowerWorkspace{TieLine})
    _update_branch_duals_single!(
        ws.vars.p_flow.duals,
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.p_flow.rho[],
        ws.interzonal,
        _dc_flow_value,
    )
    return nothing
end

# ------------------------------------------------------------
# Power — DummyBus
# ------------------------------------------------------------

function dual_update!(ws::DCPowerWorkspace{DummyBus})
    _update_branch_pairs_duals_single!(
        ws.vars.p_dummy.duals,
        ws.partitions,
        ws.vars.p_dummy.consensus,
        ws.vars.p_dummy.rho[],
        ws.branch_pairs,
        _dc_flow_value,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function dual_update!(ws::DCHybridWorkspace{TieLine})
    _update_bus_duals_single!(
        ws.vars.theta.duals,
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.theta.rho[],
        snapshot_theta_shared!,
    )
    _update_branch_duals_single!(
        ws.vars.p_flow.duals,
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.p_flow.rho[],
        ws.interzonal,
        _dc_flow_value,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function dual_update!(ws::DCHybridWorkspace{DummyBus})
    _update_bus_duals_single!(
        ws.vars.theta.duals,
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.theta.rho[],
        snapshot_theta_shared!,
    )
    _update_branch_pairs_duals_single!(
        ws.vars.p_dummy.duals,
        ws.partitions,
        ws.vars.p_dummy.consensus,
        ws.vars.p_dummy.rho[],
        ws.branch_pairs,
        _dc_flow_value,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function dual_update!(ws::DCHybridWorkspace{SplitBus})
    _update_bus_duals_single!(
        ws.vars.theta.duals,
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.theta.rho[],
        snapshot_theta_shared!,
    )
    _update_inj_groups_duals!(
        ws.vars.p_inj.duals,
        ws.partitions,
        ws.vars.p_inj.consensus,
        ws.vars.p_inj.rho[],
        ws.inj_groups,
        ws.inj_copy_zone,
        ws.p_inj_vars,
    )
    return nothing
end

# ============================================================
# AC
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function dual_update!(ws::ACStateWorkspace)
    _update_bus_duals_ac!(
        ws.vars.theta.duals,
        ws.vars.v_mag.duals,
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
        ws.vars.theta.rho[],
        ws.vars.v_mag.rho[],
    )
    return nothing
end

# ------------------------------------------------------------
# Power — TieLine
# ------------------------------------------------------------

function dual_update!(ws::ACPowerWorkspace{TieLine})
    rho_p = ws.vars.p_flow.rho[]
    rho_q = ws.vars.q_flow.rho[]

    _update_branch_duals_ac!(
        ws.vars.p_flow.duals,
        ws.vars.q_flow.duals,
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        rho_p, rho_q,
        ws.interzonal,
    )
    return nothing
end

# ------------------------------------------------------------
# Power — DummyBus
# ------------------------------------------------------------

function dual_update!(ws::ACPowerWorkspace{DummyBus})
    rho_p = ws.vars.p_dummy.rho[]
    rho_q = ws.vars.q_dummy.rho[]

    _update_branch_pairs_duals_ac!(
        ws.vars.p_dummy.duals,
        ws.vars.q_dummy.duals,
        ws.partitions,
        ws.vars.p_dummy.consensus,
        ws.vars.q_dummy.consensus,
        rho_p, rho_q,
        ws.branch_pairs,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function dual_update!(ws::ACHybridWorkspace{TieLine})
    rho_p = ws.vars.p_flow.rho[]
    rho_q = ws.vars.q_flow.rho[]

    _update_bus_duals_ac!(
        ws.vars.theta.duals,
        ws.vars.v_mag.duals,
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
        ws.vars.theta.rho[],
        ws.vars.v_mag.rho[],
    )
    _update_branch_duals_ac!(
        ws.vars.p_flow.duals,
        ws.vars.q_flow.duals,
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        rho_p, rho_q,
        ws.interzonal,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function dual_update!(ws::ACHybridWorkspace{DummyBus})
    rho_p = ws.vars.p_dummy.rho[]
    rho_q = ws.vars.q_dummy.rho[]

    _update_bus_duals_ac!(
        ws.vars.theta.duals,
        ws.vars.v_mag.duals,
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
        ws.vars.theta.rho[],
        ws.vars.v_mag.rho[],
    )
    _update_branch_pairs_duals_ac!(
        ws.vars.p_dummy.duals,
        ws.vars.q_dummy.duals,
        ws.partitions,
        ws.vars.p_dummy.consensus,
        ws.vars.q_dummy.consensus,
        rho_p, rho_q,
        ws.branch_pairs,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function dual_update!(ws::ACHybridWorkspace{SplitBus})
    _update_bus_duals_ac!(
        ws.vars.theta.duals,
        ws.vars.v_mag.duals,
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
        ws.vars.theta.rho[],
        ws.vars.v_mag.rho[],
    )
    _update_inj_groups_duals!(
        ws.vars.p_inj.duals, ws.partitions,
        ws.vars.p_inj.consensus, ws.vars.p_inj.rho[],
        ws.inj_groups, ws.inj_copy_zone, ws.p_inj_vars,
    )
    _update_inj_groups_duals!(
        ws.vars.q_inj.duals, ws.partitions,
        ws.vars.q_inj.consensus, ws.vars.q_inj.rho[],
        ws.inj_groups, ws.inj_copy_zone, ws.q_inj_vars,
    )
    return nothing
end

# ============================================================
# SOCWR
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function dual_update!(ws::SOCWRStateWorkspace)
    _update_bus_duals_single!(
        ws.vars.v2.duals,
        ws.partitions,
        ws.vars.v2.consensus,
        ws.vars.v2.rho[],
        snapshot_v2_shared!,
    )
    _update_branch_duals_paired!(
        ws.vars.r.duals,
        ws.vars.t.duals,
        ws.partitions,
        ws.vars.r.consensus,
        ws.vars.t.consensus,
        ws.vars.r.rho[],
        ws.vars.t.rho[],
        ws.interzonal,
        _socwr_rt_values,
    )
    return nothing
end

# ------------------------------------------------------------
# Power — TieLine
# ------------------------------------------------------------

function dual_update!(ws::SOCWRPowerWorkspace{TieLine})
    rho_p = ws.vars.p_flow.rho[]
    rho_q = ws.vars.q_flow.rho[]

    _update_branch_duals_paired!(
        ws.vars.p_flow.duals,
        ws.vars.q_flow.duals,
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        rho_p, rho_q,
        ws.interzonal,
        _socwr_pq_values,
    )
    return nothing
end

# ------------------------------------------------------------
# Power — DummyBus
# ------------------------------------------------------------

function dual_update!(ws::SOCWRPowerWorkspace{DummyBus})
    rho_p = ws.vars.p_dummy.rho[]
    rho_q = ws.vars.q_dummy.rho[]

    _update_branch_pairs_duals_single!(
        ws.vars.p_dummy.duals,
        ws.partitions,
        ws.vars.p_dummy.consensus,
        rho_p,
        ws.branch_pairs,
        _socwr_p_value,
    )
    _update_branch_pairs_duals_single!(
        ws.vars.q_dummy.duals,
        ws.partitions,
        ws.vars.q_dummy.consensus,
        rho_q,
        ws.branch_pairs,
        _socwr_q_value,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function dual_update!(ws::SOCWRHybridWorkspace{TieLine})
    rho_p = ws.vars.p_flow.rho[]
    rho_q = ws.vars.q_flow.rho[]

    _update_bus_duals_single!(
        ws.vars.v2.duals,
        ws.partitions,
        ws.vars.v2.consensus,
        ws.vars.v2.rho[],
        snapshot_v2_shared!,
    )
    _update_branch_duals_paired!(
        ws.vars.p_flow.duals,
        ws.vars.q_flow.duals,
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        rho_p, rho_q,
        ws.interzonal,
        _socwr_pq_values,
    )
    _update_branch_duals_paired!(
        ws.vars.r.duals,
        ws.vars.t.duals,
        ws.partitions,
        ws.vars.r.consensus,
        ws.vars.t.consensus,
        ws.vars.r.rho[],
        ws.vars.t.rho[],
        ws.interzonal,
        _socwr_rt_values,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function dual_update!(ws::SOCWRHybridWorkspace{DummyBus})
    rho_p = ws.vars.p_dummy.rho[]
    rho_q = ws.vars.q_dummy.rho[]

    _update_bus_duals_single!(
        ws.vars.v2.duals,
        ws.partitions,
        ws.vars.v2.consensus,
        ws.vars.v2.rho[],
        snapshot_v2_shared!,
    )
    _update_branch_pairs_duals_single!(
        ws.vars.p_dummy.duals,
        ws.partitions,
        ws.vars.p_dummy.consensus,
        rho_p,
        ws.branch_pairs,
        _socwr_p_value,
    )
    _update_branch_pairs_duals_single!(
        ws.vars.q_dummy.duals,
        ws.partitions,
        ws.vars.q_dummy.consensus,
        rho_q,
        ws.branch_pairs,
        _socwr_q_value,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function dual_update!(ws::SOCWRHybridWorkspace{SplitBus})
    _update_bus_duals_single!(
        ws.vars.v2.duals,
        ws.partitions,
        ws.vars.v2.consensus,
        ws.vars.v2.rho[],
        snapshot_v2_shared!,
    )
    _update_inj_groups_duals!(
        ws.vars.p_inj.duals, ws.partitions,
        ws.vars.p_inj.consensus, ws.vars.p_inj.rho[],
        ws.inj_groups, ws.inj_copy_zone, ws.p_inj_vars,
    )
    _update_inj_groups_duals!(
        ws.vars.q_inj.duals, ws.partitions,
        ws.vars.q_inj.consensus, ws.vars.q_inj.rho[],
        ws.inj_groups, ws.inj_copy_zone, ws.q_inj_vars,
    )
    return nothing
end
