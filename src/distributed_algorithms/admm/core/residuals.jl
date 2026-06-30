# ============================================================
# Helpers
# ============================================================

@inline _wrap_angle_diff(Δ::Float64) = atan(sin(Δ), cos(Δ))

# ------------------------------------------------------------
# Bus
# ------------------------------------------------------------

function _primal_residual_bus_single(
    partitions::Partitions,
    consensus::Vector{Float64},
    snapshot_fn,
)
    acc = 0.0
    for part in partitions
        buf = snapshot_fn(part)
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            li = part.shared_bus_local_idx[i]
            cj = part.shared_bus_consensus_idx[i]
            acc += (buf[li] - consensus[cj])^2
        end
    end
    return sqrt(acc)
end

function _primal_residual_bus_ac_components(
    partitions::Partitions{PartitionAC},
    theta_consensus::Vector{Float64},
    v_mag_consensus::Vector{Float64},
)
    r_theta = 0.0
    r_v_mag = 0.0
    for part in partitions
        theta_buf = snapshot_theta!(part)
        v_mag_buf = snapshot_v_mag!(part)
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            li = part.shared_bus_local_idx[i]
            cj = part.shared_bus_consensus_idx[i]
            r_theta += _wrap_angle_diff(theta_buf[li] - theta_consensus[cj])^2
            r_v_mag += (v_mag_buf[li] - v_mag_consensus[cj])^2
        end
    end
    return sqrt(r_theta), sqrt(r_v_mag)
end

function _dual_residual_bus_single(
    partitions::Partitions,
    consensus::Vector{Float64},
    consensus_prev::Vector{Float64},
    rho::Float64,
)
    s = 0.0
    for part in partitions
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = part.shared_bus_consensus_idx[i]
            s += (consensus[j] - consensus_prev[j])^2
        end
    end
    return rho * sqrt(s)
end

function _dual_residual_bus_ac_components(
    partitions::Partitions{PartitionAC},
    theta_consensus::Vector{Float64},
    v_mag_consensus::Vector{Float64},
    theta_consensus_prev::Vector{Float64},
    v_mag_consensus_prev::Vector{Float64},
    rho_theta::Float64,
    rho_v_mag::Float64,
)
    s_theta = 0.0
    s_v_mag = 0.0
    for part in partitions
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = part.shared_bus_consensus_idx[i]
            s_theta += _wrap_angle_diff(theta_consensus[j] - theta_consensus_prev[j])^2
            s_v_mag += (v_mag_consensus[j] - v_mag_consensus_prev[j])^2
        end
    end
    return rho_theta * sqrt(s_theta), rho_v_mag * sqrt(s_v_mag)
end

# ------------------------------------------------------------
# Branch
# ------------------------------------------------------------

function _primal_residual_branch_single(
    partitions::Partitions,
    consensus::AbstractVector{Float64},
    interzonal::Vector{Int},
    extract_fn,
)
    acc = 0.0
    for part in partitions
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            acc += (extract_fn(part, l) - consensus[k])^2
        end
    end
    return sqrt(acc)
end

function _primal_residual_branch_paired_components(
    partitions::Partitions,
    consensus_a::AbstractVector{Float64},
    consensus_b::AbstractVector{Float64},
    interzonal::Vector{Int},
    extract_fn,
)
    acc_a = 0.0
    acc_b = 0.0
    for part in partitions
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            a, b = extract_fn(part, l)
            acc_a += (a - consensus_a[k])^2
            acc_b += (b - consensus_b[k])^2
        end
    end
    return sqrt(acc_a), sqrt(acc_b)
end

function _primal_residual_branch_ac_components(
    partitions::Partitions{PartitionAC},
    p_consensus::AbstractVector{Float64},
    q_consensus::AbstractVector{Float64},
    interzonal::Vector{Int},
)
    acc_p = 0.0
    acc_q = 0.0
    for part in partitions
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            p, q = _pq_flow_from_exprs(part, l; t = 1)
            p_val = JuMP.value(p)
            q_val = JuMP.value(q)
            isfinite(p_val) && (acc_p += (p_val - p_consensus[k])^2)
            isfinite(q_val) && (acc_q += (q_val - q_consensus[k])^2)
        end
    end
    return sqrt(acc_p), sqrt(acc_q)
end

function _dual_residual_branch_single(
    partitions::Partitions,
    consensus::AbstractVector{Float64},
    consensus_prev::AbstractVector{Float64},
    interzonal::Vector{Int},
    rho::Float64,
)
    s = 0.0
    for part in partitions
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            s += (consensus[k] - consensus_prev[k])^2
        end
    end
    return rho * sqrt(s)
end

function _dual_residual_branch_paired_components(
    partitions::Partitions,
    consensus_a::AbstractVector{Float64},
    consensus_b::AbstractVector{Float64},
    consensus_a_prev::AbstractVector{Float64},
    consensus_b_prev::AbstractVector{Float64},
    interzonal::Vector{Int},
    rho_a::Float64,
    rho_b::Float64,
)
    s_a = 0.0
    s_b = 0.0
    for part in partitions
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            diff_a = consensus_a[k] - consensus_a_prev[k]
            diff_b = consensus_b[k] - consensus_b_prev[k]
            s_a += diff_a^2
            s_b += diff_b^2
        end
    end
    return rho_a * sqrt(s_a), rho_b * sqrt(s_b)
end

# ------------------------------------------------------------
# Branch pairs
# ------------------------------------------------------------

function _primal_residual_branch_pairs(
    partitions::Partitions,
    consensus::AbstractVector{Float64},
    pairs::Vector{NTuple{4, Int}},
    extract_fn,
)
    isempty(pairs) && return 0.0
    acc = 0.0
    @inbounds for (k, (ga, za, gb, zb)) in enumerate(pairs)
        pa = partitions.parts[za]
        la = pa.branch_idx_full_to_zone[ga]
        if la != 0
            acc += (extract_fn(pa, la) - consensus[k])^2
        end
        pb = partitions.parts[zb]
        lb = pb.branch_idx_full_to_zone[gb]
        if lb != 0
            acc += (extract_fn(pb, lb) - consensus[k])^2
        end
    end
    return sqrt(acc)
end

function _primal_residual_branch_pairs_ac(
    partitions::Partitions{PartitionAC},
    p_consensus::AbstractVector{Float64},
    q_consensus::AbstractVector{Float64},
    pairs::Vector{NTuple{4, Int}},
)
    isempty(pairs) && return (0.0, 0.0)
    acc_p = 0.0
    acc_q = 0.0
    @inbounds for (k, (ga, za, gb, zb)) in enumerate(pairs)
        pa = partitions.parts[za]
        pb = partitions.parts[zb]
        la = pa.branch_idx_full_to_zone[ga]
        lb = pb.branch_idx_full_to_zone[gb]
        if la != 0
            p_expr, q_expr = _pq_flow_from_exprs(pa, la)
            p_val = JuMP.value(p_expr)
            q_val = JuMP.value(q_expr)
            isfinite(p_val) && (acc_p += (p_val - p_consensus[k])^2)
            isfinite(q_val) && (acc_q += (q_val - q_consensus[k])^2)
        end
        if lb != 0
            p_expr, q_expr = _pq_flow_from_exprs(pb, lb)
            p_val = JuMP.value(p_expr)
            q_val = JuMP.value(q_expr)
            isfinite(p_val) && (acc_p += (p_val - p_consensus[k])^2)
            isfinite(q_val) && (acc_q += (q_val - q_consensus[k])^2)
        end
    end
    return sqrt(acc_p), sqrt(acc_q)
end

function _dual_residual_branch_pairs(
    consensus::AbstractVector{Float64},
    consensus_prev::AbstractVector{Float64},
    pairs::Vector{NTuple{4, Int}},
    rho::Float64,
)
    isempty(pairs) && return 0.0
    s = 0.0
    @inbounds for k in eachindex(pairs)
        d = consensus[k] - consensus_prev[k]
        s += 2 * d^2
    end
    return rho * sqrt(s)
end

# ============================================================
# DC
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function compute_residuals!(ws::DCStateWorkspace)
    r = _primal_residual_bus_single(
        ws.partitions,
        ws.vars.theta.consensus,
        snapshot_theta!,
    )
    s = _dual_residual_bus_single(
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.theta.consensus_prev,
        ws.vars.theta.rho[],
    )
    ws.r_primal = r
    ws.r_dual = s
    return nothing
end

# ------------------------------------------------------------
# Power — TieLine
# ------------------------------------------------------------

function compute_residuals!(ws::DCPowerWorkspace{TieLine})
    rho = ws.vars.p_flow.rho[]
    ws.r_primal = _primal_residual_branch_single(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.interzonal,
        _dc_flow_value,
    )
    ws.r_dual = _dual_residual_branch_single(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.p_flow.consensus_prev,
        ws.interzonal,
        rho,
    )
    return nothing
end

# ------------------------------------------------------------
# Power — DummyBus
# ------------------------------------------------------------

function compute_residuals!(ws::DCPowerWorkspace{DummyBus})
    rho = ws.vars.p_dummy.rho[]
    ws.r_primal = _primal_residual_branch_pairs(
        ws.partitions,
        ws.vars.p_dummy.consensus,
        ws.branch_pairs,
        _dc_flow_value,
    )
    ws.r_dual = _dual_residual_branch_pairs(
        ws.vars.p_dummy.consensus,
        ws.vars.p_dummy.consensus_prev,
        ws.branch_pairs,
        rho,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function compute_residuals!(ws::DCHybridWorkspace{TieLine})
    ws.r_primal_bus = _primal_residual_bus_single(
        ws.partitions,
        ws.vars.theta.consensus,
        snapshot_theta!,
    )
    ws.r_dual_bus = _dual_residual_bus_single(
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.theta.consensus_prev,
        ws.vars.theta.rho[],
    )

    rho = ws.vars.p_flow.rho[]
    ws.r_primal_branch = _primal_residual_branch_single(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.interzonal,
        _dc_flow_value,
    )
    ws.r_dual_branch = _dual_residual_branch_single(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.p_flow.consensus_prev,
        ws.interzonal,
        rho,
    )

    ws.r_primal = hypot(ws.r_primal_bus, ws.r_primal_branch)
    ws.r_dual = hypot(ws.r_dual_bus, ws.r_dual_branch)
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function compute_residuals!(ws::DCHybridWorkspace{DummyBus})
    ws.r_primal_bus = _primal_residual_bus_single(
        ws.partitions,
        ws.vars.theta.consensus,
        snapshot_theta!,
    )
    ws.r_dual_bus = _dual_residual_bus_single(
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.theta.consensus_prev,
        ws.vars.theta.rho[],
    )

    rho = ws.vars.p_dummy.rho[]
    ws.r_primal_branch = _primal_residual_branch_pairs(
        ws.partitions,
        ws.vars.p_dummy.consensus,
        ws.branch_pairs,
        _dc_flow_value,
    )
    ws.r_dual_branch = _dual_residual_branch_pairs(
        ws.vars.p_dummy.consensus,
        ws.vars.p_dummy.consensus_prev,
        ws.branch_pairs,
        rho,
    )

    ws.r_primal = hypot(ws.r_primal_bus, ws.r_primal_branch)
    ws.r_dual = hypot(ws.r_dual_bus, ws.r_dual_branch)
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function compute_residuals!(ws::DCHybridWorkspace{SplitBus})
    ws.r_primal_bus = _primal_residual_bus_single(
        ws.partitions,
        ws.vars.theta.consensus,
        snapshot_theta!,
    )
    ws.r_dual_bus = _dual_residual_bus_single(
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.theta.consensus_prev,
        ws.vars.theta.rho[],
    )

    rho_pi = ws.vars.p_inj.rho[]
    ws.r_primal_bus_pair = _primal_residual_inj_groups(
        ws.partitions,
        ws.vars.p_inj.consensus,
        ws.inj_groups,
        ws.inj_copy_zone,
        ws.p_inj_vars,
    )
    ws.r_dual_bus_pair = _dual_residual_inj_groups(
        ws.vars.p_inj.consensus,
        ws.vars.p_inj.consensus_prev,
        ws.inj_groups,
        rho_pi,
    )

    ws.r_primal = hypot(ws.r_primal_bus, ws.r_primal_bus_pair)
    ws.r_dual = hypot(ws.r_dual_bus, ws.r_dual_bus_pair)
    return nothing
end

# ============================================================
# AC
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function compute_residuals!(ws::ACStateWorkspace)
    ws.r_primal_theta, ws.r_primal_v_mag = _primal_residual_bus_ac_components(
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
    )
    ws.r_dual_theta, ws.r_dual_v_mag = _dual_residual_bus_ac_components(
        ws.partitions,
        ws.vars.theta.consensus, ws.vars.v_mag.consensus,
        ws.vars.theta.consensus_prev, ws.vars.v_mag.consensus_prev,
        ws.vars.theta.rho[], ws.vars.v_mag.rho[],
    )

    ws.r_primal = hypot(ws.r_primal_theta, ws.r_primal_v_mag)
    ws.r_dual = hypot(ws.r_dual_theta, ws.r_dual_v_mag)
    return nothing
end

# ------------------------------------------------------------
# Power — TieLine
# ------------------------------------------------------------

function compute_residuals!(ws::ACPowerWorkspace{TieLine})
    rho_p = ws.vars.p_flow.rho[]
    rho_q = ws.vars.q_flow.rho[]

    ws.r_primal_p, ws.r_primal_q = _primal_residual_branch_ac_components(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.interzonal,
    )
    ws.r_dual_p, ws.r_dual_q = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.vars.p_flow.consensus_prev,
        ws.vars.q_flow.consensus_prev,
        ws.interzonal, rho_p, rho_q,
    )

    ws.r_primal = hypot(ws.r_primal_p, ws.r_primal_q)
    ws.r_dual = hypot(ws.r_dual_p, ws.r_dual_q)
    return nothing
end

# ------------------------------------------------------------
# Power — DummyBus
# ------------------------------------------------------------

function compute_residuals!(ws::ACPowerWorkspace{DummyBus})
    rho_p = ws.vars.p_dummy.rho[]
    rho_q = ws.vars.q_dummy.rho[]

    ws.r_primal_p, ws.r_primal_q = _primal_residual_branch_pairs_ac(
        ws.partitions,
        ws.vars.p_dummy.consensus,
        ws.vars.q_dummy.consensus,
        ws.branch_pairs,
    )
    ws.r_dual_p = _dual_residual_branch_pairs(
        ws.vars.p_dummy.consensus,
        ws.vars.p_dummy.consensus_prev,
        ws.branch_pairs, rho_p,
    )
    ws.r_dual_q = _dual_residual_branch_pairs(
        ws.vars.q_dummy.consensus,
        ws.vars.q_dummy.consensus_prev,
        ws.branch_pairs, rho_q,
    )

    ws.r_primal = hypot(ws.r_primal_p, ws.r_primal_q)
    ws.r_dual = hypot(ws.r_dual_p, ws.r_dual_q)
    return nothing
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function compute_residuals!(ws::ACHybridWorkspace{TieLine})
    rho_p = ws.vars.p_flow.rho[]
    rho_q = ws.vars.q_flow.rho[]

    ws.r_primal_theta, ws.r_primal_v_mag = _primal_residual_bus_ac_components(
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
    )
    ws.r_dual_theta, ws.r_dual_v_mag = _dual_residual_bus_ac_components(
        ws.partitions,
        ws.vars.theta.consensus, ws.vars.v_mag.consensus,
        ws.vars.theta.consensus_prev, ws.vars.v_mag.consensus_prev,
        ws.vars.theta.rho[], ws.vars.v_mag.rho[],
    )

    ws.r_primal_p, ws.r_primal_q = _primal_residual_branch_ac_components(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.interzonal,
    )
    ws.r_dual_p, ws.r_dual_q = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.vars.p_flow.consensus_prev,
        ws.vars.q_flow.consensus_prev,
        ws.interzonal, rho_p, rho_q,
    )

    ws.r_primal = sqrt(
        ws.r_primal_theta^2 + ws.r_primal_v_mag^2 +
        ws.r_primal_p^2 + ws.r_primal_q^2,
    )
    ws.r_dual = sqrt(
        ws.r_dual_theta^2 + ws.r_dual_v_mag^2 +
        ws.r_dual_p^2 + ws.r_dual_q^2,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function compute_residuals!(ws::ACHybridWorkspace{DummyBus})
    rho_p = ws.vars.p_dummy.rho[]
    rho_q = ws.vars.q_dummy.rho[]

    ws.r_primal_theta, ws.r_primal_v_mag = _primal_residual_bus_ac_components(
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
    )
    ws.r_dual_theta, ws.r_dual_v_mag = _dual_residual_bus_ac_components(
        ws.partitions,
        ws.vars.theta.consensus, ws.vars.v_mag.consensus,
        ws.vars.theta.consensus_prev, ws.vars.v_mag.consensus_prev,
        ws.vars.theta.rho[], ws.vars.v_mag.rho[],
    )

    ws.r_primal_p, ws.r_primal_q = _primal_residual_branch_pairs_ac(
        ws.partitions,
        ws.vars.p_dummy.consensus,
        ws.vars.q_dummy.consensus,
        ws.branch_pairs,
    )
    ws.r_dual_p = _dual_residual_branch_pairs(
        ws.vars.p_dummy.consensus,
        ws.vars.p_dummy.consensus_prev,
        ws.branch_pairs, rho_p,
    )
    ws.r_dual_q = _dual_residual_branch_pairs(
        ws.vars.q_dummy.consensus,
        ws.vars.q_dummy.consensus_prev,
        ws.branch_pairs, rho_q,
    )

    ws.r_primal = sqrt(
        ws.r_primal_theta^2 + ws.r_primal_v_mag^2 +
        ws.r_primal_p^2 + ws.r_primal_q^2,
    )
    ws.r_dual = sqrt(
        ws.r_dual_theta^2 + ws.r_dual_v_mag^2 +
        ws.r_dual_p^2 + ws.r_dual_q^2,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function compute_residuals!(ws::ACHybridWorkspace{SplitBus})
    ws.r_primal_theta, ws.r_primal_v_mag = _primal_residual_bus_ac_components(
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
    )
    ws.r_dual_theta, ws.r_dual_v_mag = _dual_residual_bus_ac_components(
        ws.partitions,
        ws.vars.theta.consensus, ws.vars.v_mag.consensus,
        ws.vars.theta.consensus_prev, ws.vars.v_mag.consensus_prev,
        ws.vars.theta.rho[], ws.vars.v_mag.rho[],
    )

    ws.r_primal_p_inj = _primal_residual_inj_groups(
        ws.partitions, ws.vars.p_inj.consensus, ws.inj_groups, ws.inj_copy_zone,
        ws.p_inj_vars,
    )
    ws.r_dual_p_inj = _dual_residual_inj_groups(
        ws.vars.p_inj.consensus, ws.vars.p_inj.consensus_prev,
        ws.inj_groups, ws.vars.p_inj.rho[],
    )
    ws.r_primal_q_inj = _primal_residual_inj_groups(
        ws.partitions, ws.vars.q_inj.consensus, ws.inj_groups, ws.inj_copy_zone,
        ws.q_inj_vars,
    )
    ws.r_dual_q_inj = _dual_residual_inj_groups(
        ws.vars.q_inj.consensus, ws.vars.q_inj.consensus_prev,
        ws.inj_groups, ws.vars.q_inj.rho[],
    )

    ws.r_primal = sqrt(
        ws.r_primal_theta^2 + ws.r_primal_v_mag^2 +
        ws.r_primal_p_inj^2 + ws.r_primal_q_inj^2,
    )
    ws.r_dual = sqrt(
        ws.r_dual_theta^2 + ws.r_dual_v_mag^2 +
        ws.r_dual_p_inj^2 + ws.r_dual_q_inj^2,
    )
    return nothing
end

# ============================================================
# SOCWR
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function compute_residuals!(ws::SOCWRStateWorkspace)
    ws.r_primal_bus = _primal_residual_bus_single(
        ws.partitions,
        ws.vars.v2.consensus,
        snapshot_v2!,
    )
    ws.r_dual_bus = _dual_residual_bus_single(
        ws.partitions,
        ws.vars.v2.consensus,
        ws.vars.v2.consensus_prev,
        ws.vars.v2.rho[],
    )
    ws.r_primal_r, ws.r_primal_t = _primal_residual_branch_paired_components(
        ws.partitions,
        ws.vars.r.consensus,
        ws.vars.t.consensus,
        ws.interzonal,
        _socwr_rt_values,
    )
    ws.r_dual_r, ws.r_dual_t = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.vars.r.consensus, ws.vars.t.consensus,
        ws.vars.r.consensus_prev, ws.vars.t.consensus_prev,
        ws.interzonal, ws.vars.r.rho[], ws.vars.t.rho[],
    )

    ws.r_primal = sqrt(ws.r_primal_bus^2 + ws.r_primal_r^2 + ws.r_primal_t^2)
    ws.r_dual = sqrt(ws.r_dual_bus^2 + ws.r_dual_r^2 + ws.r_dual_t^2)
    return nothing
end

# ------------------------------------------------------------
# Power — TieLine
# ------------------------------------------------------------

function compute_residuals!(ws::SOCWRPowerWorkspace{TieLine})
    rho_p = ws.vars.p_flow.rho[]
    rho_q = ws.vars.q_flow.rho[]

    ws.r_primal_p, ws.r_primal_q = _primal_residual_branch_paired_components(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.interzonal,
        _socwr_pq_values,
    )
    ws.r_dual_p, ws.r_dual_q = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.vars.p_flow.consensus_prev,
        ws.vars.q_flow.consensus_prev,
        ws.interzonal, rho_p, rho_q,
    )

    ws.r_primal = hypot(ws.r_primal_p, ws.r_primal_q)
    ws.r_dual = hypot(ws.r_dual_p, ws.r_dual_q)
    return nothing
end

# ------------------------------------------------------------
# Power — DummyBus
# ------------------------------------------------------------

function compute_residuals!(ws::SOCWRPowerWorkspace{DummyBus})
    rho_p = ws.vars.p_dummy.rho[]
    rho_q = ws.vars.q_dummy.rho[]

    ws.r_primal_p = _primal_residual_branch_pairs(
        ws.partitions,
        ws.vars.p_dummy.consensus,
        ws.branch_pairs,
        _socwr_p_value,
    )
    ws.r_primal_q = _primal_residual_branch_pairs(
        ws.partitions,
        ws.vars.q_dummy.consensus,
        ws.branch_pairs,
        _socwr_q_value,
    )
    ws.r_dual_p = _dual_residual_branch_pairs(
        ws.vars.p_dummy.consensus,
        ws.vars.p_dummy.consensus_prev,
        ws.branch_pairs, rho_p,
    )
    ws.r_dual_q = _dual_residual_branch_pairs(
        ws.vars.q_dummy.consensus,
        ws.vars.q_dummy.consensus_prev,
        ws.branch_pairs, rho_q,
    )

    ws.r_primal = hypot(ws.r_primal_p, ws.r_primal_q)
    ws.r_dual = hypot(ws.r_dual_p, ws.r_dual_q)
    return nothing
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function compute_residuals!(ws::SOCWRHybridWorkspace{TieLine})
    rho_p = ws.vars.p_flow.rho[]
    rho_q = ws.vars.q_flow.rho[]

    ws.r_primal_bus = _primal_residual_bus_single(
        ws.partitions,
        ws.vars.v2.consensus,
        snapshot_v2!,
    )
    ws.r_dual_bus = _dual_residual_bus_single(
        ws.partitions,
        ws.vars.v2.consensus,
        ws.vars.v2.consensus_prev,
        ws.vars.v2.rho[],
    )

    ws.r_primal_p, ws.r_primal_q = _primal_residual_branch_paired_components(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.interzonal,
        _socwr_pq_values,
    )
    ws.r_dual_p, ws.r_dual_q = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.vars.p_flow.consensus_prev,
        ws.vars.q_flow.consensus_prev,
        ws.interzonal, rho_p, rho_q,
    )

    ws.r_primal_r, ws.r_primal_t = _primal_residual_branch_paired_components(
        ws.partitions,
        ws.vars.r.consensus,
        ws.vars.t.consensus,
        ws.interzonal,
        _socwr_rt_values,
    )
    ws.r_dual_r, ws.r_dual_t = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.vars.r.consensus, ws.vars.t.consensus,
        ws.vars.r.consensus_prev, ws.vars.t.consensus_prev,
        ws.interzonal, ws.vars.r.rho[], ws.vars.t.rho[],
    )

    ws.r_primal = sqrt(
        ws.r_primal_bus^2 +
        ws.r_primal_p^2 + ws.r_primal_q^2 +
        ws.r_primal_r^2 + ws.r_primal_t^2,
    )
    ws.r_dual = sqrt(
        ws.r_dual_bus^2 +
        ws.r_dual_p^2 + ws.r_dual_q^2 +
        ws.r_dual_r^2 + ws.r_dual_t^2,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function compute_residuals!(ws::SOCWRHybridWorkspace{DummyBus})
    rho_p = ws.vars.p_dummy.rho[]
    rho_q = ws.vars.q_dummy.rho[]

    ws.r_primal_bus = _primal_residual_bus_single(
        ws.partitions,
        ws.vars.v2.consensus,
        snapshot_v2!,
    )
    ws.r_dual_bus = _dual_residual_bus_single(
        ws.partitions,
        ws.vars.v2.consensus,
        ws.vars.v2.consensus_prev,
        ws.vars.v2.rho[],
    )

    ws.r_primal_p = _primal_residual_branch_pairs(
        ws.partitions,
        ws.vars.p_dummy.consensus,
        ws.branch_pairs,
        _socwr_p_value,
    )
    ws.r_primal_q = _primal_residual_branch_pairs(
        ws.partitions,
        ws.vars.q_dummy.consensus,
        ws.branch_pairs,
        _socwr_q_value,
    )
    ws.r_dual_p = _dual_residual_branch_pairs(
        ws.vars.p_dummy.consensus,
        ws.vars.p_dummy.consensus_prev,
        ws.branch_pairs, rho_p,
    )
    ws.r_dual_q = _dual_residual_branch_pairs(
        ws.vars.q_dummy.consensus,
        ws.vars.q_dummy.consensus_prev,
        ws.branch_pairs, rho_q,
    )

    ws.r_primal = sqrt(
        ws.r_primal_bus^2 +
        ws.r_primal_p^2 + ws.r_primal_q^2,
    )
    ws.r_dual = sqrt(
        ws.r_dual_bus^2 +
        ws.r_dual_p^2 + ws.r_dual_q^2,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function compute_residuals!(ws::SOCWRHybridWorkspace{SplitBus})
    ws.r_primal_bus = _primal_residual_bus_single(
        ws.partitions,
        ws.vars.v2.consensus,
        snapshot_v2!,
    )
    ws.r_dual_bus = _dual_residual_bus_single(
        ws.partitions,
        ws.vars.v2.consensus,
        ws.vars.v2.consensus_prev,
        ws.vars.v2.rho[],
    )

    ws.r_primal_p_inj = _primal_residual_inj_groups(
        ws.partitions, ws.vars.p_inj.consensus, ws.inj_groups, ws.inj_copy_zone,
        ws.p_inj_vars,
    )
    ws.r_dual_p_inj = _dual_residual_inj_groups(
        ws.vars.p_inj.consensus, ws.vars.p_inj.consensus_prev,
        ws.inj_groups, ws.vars.p_inj.rho[],
    )
    ws.r_primal_q_inj = _primal_residual_inj_groups(
        ws.partitions, ws.vars.q_inj.consensus, ws.inj_groups, ws.inj_copy_zone,
        ws.q_inj_vars,
    )
    ws.r_dual_q_inj = _dual_residual_inj_groups(
        ws.vars.q_inj.consensus, ws.vars.q_inj.consensus_prev,
        ws.inj_groups, ws.vars.q_inj.rho[],
    )

    ws.r_primal = sqrt(
        ws.r_primal_bus^2 +
        ws.r_primal_p_inj^2 + ws.r_primal_q_inj^2,
    )
    ws.r_dual = sqrt(
        ws.r_dual_bus^2 +
        ws.r_dual_p_inj^2 + ws.r_dual_q_inj^2,
    )
    return nothing
end
