# ============================================================
# ADMM Dual Update — unified for DC · AC · SOCWR
# ============================================================

# ------------------------------------------------------------
# Bus duals: single variable, linear (DC theta, SOCWR V²)
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

# ------------------------------------------------------------
# Bus duals: AC (angle-wrapped theta + linear |V|)
# ------------------------------------------------------------
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
# Branch duals: single variable (DC P)
# ------------------------------------------------------------
function _update_branch_duals_single!(
    duals::AbstractMatrix{Float64},
    partitions::Partitions,
    consensus::AbstractVector{Float64},
    rho::Float64,
    interzonal::Vector{Int},
    extract_fn,    # (part, l) -> Float64
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

# ------------------------------------------------------------
# Branch pair duals (DummyBus Method 2)
# ------------------------------------------------------------
function _update_branch_pairs_duals_single!(
    duals::AbstractMatrix{Float64},         # rows = pair region, cols = zones
    partitions::Partitions,
    consensus::AbstractVector{Float64},     # length == n_pairs
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

# ------------------------------------------------------------
# Bus pair duals (SplitBus Method 3): sign-flipped consensus.
# Zone A: λ_A += ρ · (p_A - z)
# Zone B: λ_B += ρ · (p_B + z)   (target is -z ⇒ residual is p_B + z)
# ------------------------------------------------------------
function _update_bus_pairs_duals_signed!(
    duals::AbstractMatrix{Float64},       # rows = bus_pairs, cols = zones
    partitions::Partitions,
    consensus::AbstractVector{Float64},   # length == n_bus_pairs
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

# ------------------------------------------------------------
# Branch duals: paired variables (SOCWR PQ, SOCWR RT)
# ------------------------------------------------------------
function _update_branch_duals_paired!(
    duals_a::AbstractMatrix{Float64},
    duals_b::AbstractMatrix{Float64},
    partitions::Partitions,
    consensus_a::AbstractVector{Float64},
    consensus_b::AbstractVector{Float64},
    rho_a::Float64,
    rho_b::Float64,
    interzonal::Vector{Int},
    extract_fn,    # (part, l) -> (Float64, Float64)
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

# ------------------------------------------------------------
# Branch duals: AC P/Q (with isfinite checks for NLP)
# ------------------------------------------------------------
# Specialized variant because AC flow expressions can return NaN; a bad
# P must NOT skip a good Q update (and vice versa).
# NaN-safe branch_pairs dual update for AC (DummyBus).
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

# ============================================================
# Dispatch: dual_update!(ws)
# ============================================================

# ---- DC ----

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

function dual_update!(ws::DCPowerWorkspace)
    n_inter = length(ws.interzonal)
    n_pairs = length(ws.branch_pairs)
    if n_inter > 0
        _update_branch_duals_single!(
            view(ws.vars.p_flow.duals, 1:n_inter, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, 1:n_inter),
            ws.vars.p_flow.rho[],
            ws.interzonal,
            _dc_flow_value,
        )
    end
    if n_pairs > 0
        _update_branch_pairs_duals_single!(
            view(ws.vars.p_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, n_inter+1:n_inter+n_pairs),
            ws.vars.p_flow.rho[],
            ws.branch_pairs,
            _dc_flow_value,
        )
    end
    return nothing
end

function dual_update!(ws::DCHybridWorkspace{<:Union{TieLine, DummyBus}})
    _update_bus_duals_single!(
        ws.vars.theta.duals,
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.theta.rho[],
        snapshot_theta_shared!,
    )
    n_inter = length(ws.interzonal)
    n_pairs = length(ws.branch_pairs)
    if n_inter > 0
        _update_branch_duals_single!(
            view(ws.vars.p_flow.duals, 1:n_inter, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, 1:n_inter),
            ws.vars.p_flow.rho[],
            ws.interzonal,
            _dc_flow_value,
        )
    end
    if n_pairs > 0
        _update_branch_pairs_duals_single!(
            view(ws.vars.p_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, n_inter+1:n_inter+n_pairs),
            ws.vars.p_flow.rho[],
            ws.branch_pairs,
            _dc_flow_value,
        )
    end
    return nothing
end

function dual_update!(ws::DCHybridWorkspace{SplitBus})
    _update_bus_duals_single!(
        ws.vars.theta.duals,
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.theta.rho[],
        snapshot_theta_shared!,
    )
    n_inter = length(ws.interzonal)
    if n_inter > 0
        _update_branch_duals_single!(
            view(ws.vars.p_flow.duals, 1:n_inter, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, 1:n_inter),
            ws.vars.p_flow.rho[],
            ws.interzonal,
            _dc_flow_value,
        )
    end
    _update_bus_pairs_duals_signed!(
        ws.vars.p_inj.duals,
        ws.partitions,
        ws.vars.p_inj.consensus,
        ws.vars.p_inj.rho[],
        ws.bus_pairs,
        ws.p_inj_vars,
    )
    return nothing
end

# ---- AC ----

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

function dual_update!(ws::ACPowerWorkspace)
    n_inter = length(ws.interzonal)
    n_pairs = length(ws.branch_pairs)
    rho_p = ws.vars.p_flow.rho[]
    rho_q = ws.vars.q_flow.rho[]

    if n_inter > 0
        _update_branch_duals_ac!(
            view(ws.vars.p_flow.duals, 1:n_inter, :),
            view(ws.vars.q_flow.duals, 1:n_inter, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, 1:n_inter),
            view(ws.vars.q_flow.consensus, 1:n_inter),
            rho_p, rho_q,
            ws.interzonal,
        )
    end
    if n_pairs > 0
        _update_branch_pairs_duals_ac!(
            view(ws.vars.p_flow.duals, n_inter+1:n_inter+n_pairs, :),
            view(ws.vars.q_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, n_inter+1:n_inter+n_pairs),
            view(ws.vars.q_flow.consensus, n_inter+1:n_inter+n_pairs),
            rho_p, rho_q,
            ws.branch_pairs,
        )
    end
    return nothing
end

function dual_update!(ws::ACHybridWorkspace{<:Union{TieLine, DummyBus}})
    n_inter = length(ws.interzonal)
    n_pairs = length(ws.branch_pairs)
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
    if n_inter > 0
        _update_branch_duals_ac!(
            view(ws.vars.p_flow.duals, 1:n_inter, :),
            view(ws.vars.q_flow.duals, 1:n_inter, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, 1:n_inter),
            view(ws.vars.q_flow.consensus, 1:n_inter),
            rho_p, rho_q,
            ws.interzonal,
        )
    end
    if n_pairs > 0
        _update_branch_pairs_duals_ac!(
            view(ws.vars.p_flow.duals, n_inter+1:n_inter+n_pairs, :),
            view(ws.vars.q_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, n_inter+1:n_inter+n_pairs),
            view(ws.vars.q_flow.consensus, n_inter+1:n_inter+n_pairs),
            rho_p, rho_q,
            ws.branch_pairs,
        )
    end
    return nothing
end

function dual_update!(ws::ACHybridWorkspace{SplitBus})
    n_inter = length(ws.interzonal)
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
    if n_inter > 0
        _update_branch_duals_ac!(
            view(ws.vars.p_flow.duals, 1:n_inter, :),
            view(ws.vars.q_flow.duals, 1:n_inter, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, 1:n_inter),
            view(ws.vars.q_flow.consensus, 1:n_inter),
            rho_p, rho_q,
            ws.interzonal,
        )
    end
    _update_bus_pairs_duals_signed!(
        ws.vars.p_inj.duals, ws.partitions,
        ws.vars.p_inj.consensus, ws.vars.p_inj.rho[],
        ws.bus_pairs, ws.p_inj_vars,
    )
    _update_bus_pairs_duals_signed!(
        ws.vars.q_inj.duals, ws.partitions,
        ws.vars.q_inj.consensus, ws.vars.q_inj.rho[],
        ws.bus_pairs, ws.q_inj_vars,
    )
    return nothing
end

# ---- SOCWR ----

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

function dual_update!(ws::SOCWRPowerWorkspace)
    n_inter = length(ws.interzonal)
    n_pairs = length(ws.branch_pairs)
    rho_p = ws.vars.p_flow.rho[]
    rho_q = ws.vars.q_flow.rho[]

    if n_inter > 0
        _update_branch_duals_paired!(
            view(ws.vars.p_flow.duals, 1:n_inter, :),
            view(ws.vars.q_flow.duals, 1:n_inter, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, 1:n_inter),
            view(ws.vars.q_flow.consensus, 1:n_inter),
            rho_p, rho_q,
            ws.interzonal,
            _socwr_pq_values,
        )
    end
    if n_pairs > 0
        _update_branch_pairs_duals_single!(
            view(ws.vars.p_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, n_inter+1:n_inter+n_pairs),
            rho_p,
            ws.branch_pairs,
            _socwr_p_value,
        )
        _update_branch_pairs_duals_single!(
            view(ws.vars.q_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.partitions,
            view(ws.vars.q_flow.consensus, n_inter+1:n_inter+n_pairs),
            rho_q,
            ws.branch_pairs,
            _socwr_q_value,
        )
    end
    return nothing
end

function dual_update!(ws::SOCWRHybridWorkspace{<:Union{TieLine, DummyBus}})
    n_inter = length(ws.interzonal)
    n_pairs = length(ws.branch_pairs)
    rho_p = ws.vars.p_flow.rho[]
    rho_q = ws.vars.q_flow.rho[]

    _update_bus_duals_single!(
        ws.vars.v2.duals,
        ws.partitions,
        ws.vars.v2.consensus,
        ws.vars.v2.rho[],
        snapshot_v2_shared!,
    )
    if n_inter > 0
        _update_branch_duals_paired!(
            view(ws.vars.p_flow.duals, 1:n_inter, :),
            view(ws.vars.q_flow.duals, 1:n_inter, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, 1:n_inter),
            view(ws.vars.q_flow.consensus, 1:n_inter),
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
    end
    if n_pairs > 0
        _update_branch_pairs_duals_single!(
            view(ws.vars.p_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, n_inter+1:n_inter+n_pairs),
            rho_p,
            ws.branch_pairs,
            _socwr_p_value,
        )
        _update_branch_pairs_duals_single!(
            view(ws.vars.q_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.partitions,
            view(ws.vars.q_flow.consensus, n_inter+1:n_inter+n_pairs),
            rho_q,
            ws.branch_pairs,
            _socwr_q_value,
        )
    end
    return nothing
end

function dual_update!(ws::SOCWRHybridWorkspace{SplitBus})
    n_inter = length(ws.interzonal)
    rho_p = ws.vars.p_flow.rho[]
    rho_q = ws.vars.q_flow.rho[]

    _update_bus_duals_single!(
        ws.vars.v2.duals,
        ws.partitions,
        ws.vars.v2.consensus,
        ws.vars.v2.rho[],
        snapshot_v2_shared!,
    )
    if n_inter > 0
        _update_branch_duals_paired!(
            view(ws.vars.p_flow.duals, 1:n_inter, :),
            view(ws.vars.q_flow.duals, 1:n_inter, :),
            ws.partitions,
            view(ws.vars.p_flow.consensus, 1:n_inter),
            view(ws.vars.q_flow.consensus, 1:n_inter),
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
    end
    _update_bus_pairs_duals_signed!(
        ws.vars.p_inj.duals, ws.partitions,
        ws.vars.p_inj.consensus, ws.vars.p_inj.rho[],
        ws.bus_pairs, ws.p_inj_vars,
    )
    _update_bus_pairs_duals_signed!(
        ws.vars.q_inj.duals, ws.partitions,
        ws.vars.q_inj.consensus, ws.vars.q_inj.rho[],
        ws.bus_pairs, ws.q_inj_vars,
    )
    return nothing
end
