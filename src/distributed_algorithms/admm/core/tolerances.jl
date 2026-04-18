# ============================================================
# ADMM Tolerances (Boyd et al.)
# ============================================================

# ------------------------------------------------------------
# Bus tolerance helpers
# ------------------------------------------------------------

function _eps_primal_bus(
    partitions::Partitions,
    consensus_vecs,
    snapshot_fns,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    norm_x_sq = 0.0
    norm_z_sq = 0.0
    for part in partitions
        bufs = map(fn -> fn(part), snapshot_fns)
        n += length(snapshot_fns) * length(part.zone_shared_bus_ids)
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            li = part.shared_bus_local_idx[i]
            cj = part.shared_bus_consensus_idx[i]
            for k in eachindex(bufs)
                norm_x_sq += bufs[k][li]^2
                norm_z_sq += consensus_vecs[k][cj]^2
            end
        end
    end
    n = max(n, 1)

    return sqrt(n) * eps_abs + eps_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

function _eps_dual_bus(
    duals_list,
    partitions::Partitions,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    norm_y_sq = 0.0
    for part in partitions
        p = part.position
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = part.shared_bus_consensus_idx[i]
            for y in duals_list
                n += 1
                norm_y_sq += y[j, p]^2
            end
        end
    end
    n = max(n, 1)
    return sqrt(n) * eps_abs + eps_rel * sqrt(norm_y_sq)
end

# ------------------------------------------------------------
# Branch tolerance helpers
# ------------------------------------------------------------

function _eps_primal_branch(
    partitions::Partitions,
    consensus_vecs,
    var_symbols::NTuple{N, Symbol} where {N},
    interzonal::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 0
    norm_x_sq = 0.0
    norm_z_sq = 0.0
    for part in partitions
        m = part.model
        has_vals = JuMP.has_values(m)
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            for (vi, sym) in enumerate(var_symbols)
                dim += 1
                if has_vals
                    xv = JuMP.value(m[sym][l, 1])
                    isfinite(xv) && (norm_x_sq += xv * xv)
                end
                norm_z_sq += consensus_vecs[vi][k]^2
            end
        end
    end

    return sqrt(max(1, dim)) * eps_abs +
           eps_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq), 1e-12)
end

# ------------------------------------------------------------
# Branch-pair tolerance helpers (DummyBus Method 2)
# ------------------------------------------------------------
function _eps_primal_branch_pairs(
    partitions::Partitions,
    consensus::AbstractVector{Float64},
    pairs::Vector{NTuple{4, Int}},
    sym::Symbol,
    eps_abs::Float64,
    eps_rel::Float64,
)
    isempty(pairs) && return 0.0
    dim = 0
    norm_x_sq = 0.0
    norm_z_sq = 0.0
    @inbounds for (k, (ga, za, gb, zb)) in enumerate(pairs)
        pa = partitions.parts[za]
        la = pa.branch_idx_full_to_zone[ga]
        if la != 0
            dim += 1
            m = pa.model
            if JuMP.has_values(m)
                xv = JuMP.value(m[sym][la, 1])
                isfinite(xv) && (norm_x_sq += xv * xv)
            end
            norm_z_sq += consensus[k]^2
        end
        pb = partitions.parts[zb]
        lb = pb.branch_idx_full_to_zone[gb]
        if lb != 0
            dim += 1
            m = pb.model
            if JuMP.has_values(m)
                xv = JuMP.value(m[sym][lb, 1])
                isfinite(xv) && (norm_x_sq += xv * xv)
            end
            norm_z_sq += consensus[k]^2
        end
    end
    return sqrt(max(1, dim)) * eps_abs +
           eps_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq), 1e-12)
end

function _eps_dual_branch_pairs(
    duals::AbstractMatrix{Float64},
    pairs::Vector{NTuple{4, Int}},
    eps_abs::Float64,
    eps_rel::Float64,
)
    isempty(pairs) && return 0.0
    n = 0
    norm_y_sq = 0.0
    @inbounds for (k, (_, za, _, zb)) in enumerate(pairs)
        n += 2
        norm_y_sq += duals[k, za]^2 + duals[k, zb]^2
    end
    return sqrt(max(1, n)) * eps_abs + eps_rel * sqrt(norm_y_sq)
end

# ------------------------------------------------------------
# Bus-pair tolerance helpers (SplitBus Method 3)
# ------------------------------------------------------------
function _eps_primal_bus_pairs_signed(
    partitions::Partitions,
    consensus::AbstractVector{Float64},
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
    eps_abs::Float64,
    eps_rel::Float64,
)
    isempty(bus_pairs) && return 0.0
    dim = 0
    norm_x_sq = 0.0
    norm_z_sq = 0.0
    @inbounds for (k, (bid, za, zb)) in enumerate(bus_pairs)
        va = _p_inj_value(partitions.parts[za], p_inj_vars[za], bid)
        vb = _p_inj_value(partitions.parts[zb], p_inj_vars[zb], bid)
        dim += 2
        norm_x_sq += va^2 + vb^2
        # Zone A targets +z, zone B targets -z: both contribute z² to the target norm.
        norm_z_sq += 2 * consensus[k]^2
    end
    return sqrt(max(1, dim)) * eps_abs +
           eps_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq), 1e-12)
end

function _eps_dual_bus_pairs_signed(
    duals::AbstractMatrix{Float64},
    bus_pairs::Vector{NTuple{3, Int}},
    eps_abs::Float64,
    eps_rel::Float64,
)
    isempty(bus_pairs) && return 0.0
    n = 0
    norm_y_sq = 0.0
    @inbounds for (k, (_, za, zb)) in enumerate(bus_pairs)
        n += 2
        norm_y_sq += duals[k, za]^2 + duals[k, zb]^2
    end
    return sqrt(max(1, n)) * eps_abs + eps_rel * sqrt(norm_y_sq)
end

function _eps_dual_branch(
    duals_list,         # tuple of (Abstract)Matrix{Float64} (n_inter × n_zones)
    partitions::Partitions,
    interzonal::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    norm_y_sq = 0.0
    for part in partitions
        p = part.position
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            for y in duals_list
                n += 1
                norm_y_sq += y[k, p]^2
            end
        end
    end
    n = max(n, 1)
    return sqrt(n) * eps_abs + eps_rel * sqrt(norm_y_sq)
end

# ############################################################
#  DISPATCH: compute_tolerances!(ws, eps_abs, eps_rel)
# ############################################################

# ---- DC ----

function compute_tolerances!(ws::DCStateWorkspace, eps_abs::Float64, eps_rel::Float64)
    eps_r = _eps_primal_bus(
        ws.partitions,
        (ws.vars.theta.consensus,),
        (snapshot_theta!,),
        eps_abs,
        eps_rel,
    )
    eps_s = _eps_dual_bus(
        (ws.vars.theta.duals,),
        ws.partitions,
        eps_abs,
        eps_rel,
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end

function compute_tolerances!(ws::DCPowerWorkspace, eps_abs::Float64, eps_rel::Float64)
    n_inter = length(ws.interzonal)
    n_pairs = length(ws.branch_pairs)

    eps_r_inter =
        n_inter > 0 ?
        _eps_primal_branch(
            ws.partitions,
            (view(ws.vars.p_flow.consensus, 1:n_inter),),
            (:p_flow_branches,),
            ws.interzonal,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_inter =
        n_inter > 0 ?
        _eps_dual_branch(
            (view(ws.vars.p_flow.duals, 1:n_inter, :),),
            ws.partitions,
            ws.interzonal,
            eps_abs,
            eps_rel,
        ) : 0.0

    eps_r_pair =
        n_pairs > 0 ?
        _eps_primal_branch_pairs(
            ws.partitions,
            view(ws.vars.p_flow.consensus, n_inter+1:n_inter+n_pairs),
            ws.branch_pairs,
            :p_flow_branches,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_pair =
        n_pairs > 0 ?
        _eps_dual_branch_pairs(
            view(ws.vars.p_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.branch_pairs,
            eps_abs,
            eps_rel,
        ) : 0.0

    ws.eps_primal = hypot(eps_r_inter, eps_r_pair)
    ws.eps_dual = hypot(eps_s_inter, eps_s_pair)
    return nothing
end

function compute_tolerances!(ws::DCHybridWorkspace, eps_abs::Float64, eps_rel::Float64)
    n_inter = length(ws.interzonal)
    n_pairs = length(ws.branch_pairs)

    eps_r_bus = _eps_primal_bus(
        ws.partitions,
        (ws.vars.theta.consensus,),
        (snapshot_theta!,),
        eps_abs,
        eps_rel,
    )
    eps_s_bus = _eps_dual_bus(
        (ws.vars.theta.duals,),
        ws.partitions,
        eps_abs,
        eps_rel,
    )

    eps_r_inter =
        n_inter > 0 ?
        _eps_primal_branch(
            ws.partitions,
            (view(ws.vars.p_flow.consensus, 1:n_inter),),
            (:p_flow_branches,),
            ws.interzonal,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_inter =
        n_inter > 0 ?
        _eps_dual_branch(
            (view(ws.vars.p_flow.duals, 1:n_inter, :),),
            ws.partitions,
            ws.interzonal,
            eps_abs,
            eps_rel,
        ) : 0.0

    eps_r_pair =
        n_pairs > 0 ?
        _eps_primal_branch_pairs(
            ws.partitions,
            view(ws.vars.p_flow.consensus, n_inter+1:n_inter+n_pairs),
            ws.branch_pairs,
            :p_flow_branches,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_pair =
        n_pairs > 0 ?
        _eps_dual_branch_pairs(
            view(ws.vars.p_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.branch_pairs,
            eps_abs,
            eps_rel,
        ) : 0.0

    eps_r_bpair =
        !isempty(ws.bus_pairs) ?
        _eps_primal_bus_pairs_signed(
            ws.partitions,
            ws.vars.p_inj.consensus,
            ws.bus_pairs,
            ws.p_inj_vars,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_bpair =
        !isempty(ws.bus_pairs) ?
        _eps_dual_bus_pairs_signed(
            ws.vars.p_inj.duals,
            ws.bus_pairs,
            eps_abs,
            eps_rel,
        ) : 0.0

    ws.eps_primal =
        sqrt(eps_r_bus^2 + eps_r_inter^2 + eps_r_pair^2 + eps_r_bpair^2)
    ws.eps_dual =
        sqrt(eps_s_bus^2 + eps_s_inter^2 + eps_s_pair^2 + eps_s_bpair^2)
    return nothing
end

# ---- AC ----

function compute_tolerances!(ws::ACStateWorkspace, eps_abs::Float64, eps_rel::Float64)
    eps_r = _eps_primal_bus(
        ws.partitions,
        (ws.vars.theta.consensus, ws.vars.v_mag.consensus),
        (snapshot_theta!, snapshot_v_mag!),
        eps_abs,
        eps_rel,
    )
    eps_s = _eps_dual_bus(
        (ws.vars.theta.duals, ws.vars.v_mag.duals),
        ws.partitions,
        eps_abs,
        eps_rel,
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end

function compute_tolerances!(ws::ACPowerWorkspace, eps_abs::Float64, eps_rel::Float64)
    n_inter = length(ws.interzonal)
    n_pairs = length(ws.branch_pairs)

    eps_r_inter =
        n_inter > 0 ?
        _eps_primal_branch(
            ws.partitions,
            (view(ws.vars.p_flow.consensus, 1:n_inter),
                view(ws.vars.q_flow.consensus, 1:n_inter)),
            (:p_flow_branches, :q_flow_branches),
            ws.interzonal,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_inter =
        n_inter > 0 ?
        _eps_dual_branch(
            (view(ws.vars.p_flow.duals, 1:n_inter, :),
                view(ws.vars.q_flow.duals, 1:n_inter, :)),
            ws.partitions,
            ws.interzonal,
            eps_abs,
            eps_rel,
        ) : 0.0

    eps_r_p_pair =
        n_pairs > 0 ?
        _eps_primal_branch_pairs(
            ws.partitions,
            view(ws.vars.p_flow.consensus, n_inter+1:n_inter+n_pairs),
            ws.branch_pairs,
            :p_flow_branches,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_r_q_pair =
        n_pairs > 0 ?
        _eps_primal_branch_pairs(
            ws.partitions,
            view(ws.vars.q_flow.consensus, n_inter+1:n_inter+n_pairs),
            ws.branch_pairs,
            :q_flow_branches,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_p_pair =
        n_pairs > 0 ?
        _eps_dual_branch_pairs(
            view(ws.vars.p_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.branch_pairs,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_q_pair =
        n_pairs > 0 ?
        _eps_dual_branch_pairs(
            view(ws.vars.q_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.branch_pairs,
            eps_abs,
            eps_rel,
        ) : 0.0

    ws.eps_primal = sqrt(eps_r_inter^2 + eps_r_p_pair^2 + eps_r_q_pair^2)
    ws.eps_dual = sqrt(eps_s_inter^2 + eps_s_p_pair^2 + eps_s_q_pair^2)
    return nothing
end

function compute_tolerances!(ws::ACHybridWorkspace, eps_abs::Float64, eps_rel::Float64)
    n_inter = length(ws.interzonal)
    n_pairs = length(ws.branch_pairs)

    eps_r_bus = _eps_primal_bus(
        ws.partitions,
        (ws.vars.theta.consensus, ws.vars.v_mag.consensus),
        (snapshot_theta!, snapshot_v_mag!),
        eps_abs,
        eps_rel,
    )
    eps_s_bus = _eps_dual_bus(
        (ws.vars.theta.duals, ws.vars.v_mag.duals),
        ws.partitions,
        eps_abs,
        eps_rel,
    )

    eps_r_branch_inter =
        n_inter > 0 ?
        _eps_primal_branch(
            ws.partitions,
            (view(ws.vars.p_flow.consensus, 1:n_inter),
                view(ws.vars.q_flow.consensus, 1:n_inter)),
            (:p_flow_branches, :q_flow_branches),
            ws.interzonal,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_branch_inter =
        n_inter > 0 ?
        _eps_dual_branch(
            (view(ws.vars.p_flow.duals, 1:n_inter, :),
                view(ws.vars.q_flow.duals, 1:n_inter, :)),
            ws.partitions,
            ws.interzonal,
            eps_abs,
            eps_rel,
        ) : 0.0

    eps_r_p_pair =
        n_pairs > 0 ?
        _eps_primal_branch_pairs(
            ws.partitions,
            view(ws.vars.p_flow.consensus, n_inter+1:n_inter+n_pairs),
            ws.branch_pairs,
            :p_flow_branches,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_r_q_pair =
        n_pairs > 0 ?
        _eps_primal_branch_pairs(
            ws.partitions,
            view(ws.vars.q_flow.consensus, n_inter+1:n_inter+n_pairs),
            ws.branch_pairs,
            :q_flow_branches,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_p_pair =
        n_pairs > 0 ?
        _eps_dual_branch_pairs(
            view(ws.vars.p_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.branch_pairs,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_q_pair =
        n_pairs > 0 ?
        _eps_dual_branch_pairs(
            view(ws.vars.q_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.branch_pairs,
            eps_abs,
            eps_rel,
        ) : 0.0

    eps_r_branch = sqrt(
        eps_r_branch_inter^2 + eps_r_p_pair^2 + eps_r_q_pair^2,
    )
    eps_s_branch = sqrt(
        eps_s_branch_inter^2 + eps_s_p_pair^2 + eps_s_q_pair^2,
    )

    if !isempty(ws.bus_pairs)
        eps_r_p_inj = _eps_primal_bus_pairs_signed(
            ws.partitions, ws.vars.p_inj.consensus, ws.bus_pairs,
            ws.p_inj_vars, eps_abs, eps_rel,
        )
        eps_s_p_inj = _eps_dual_bus_pairs_signed(
            ws.vars.p_inj.duals, ws.bus_pairs, eps_abs, eps_rel,
        )
        eps_r_q_inj = _eps_primal_bus_pairs_signed(
            ws.partitions, ws.vars.q_inj.consensus, ws.bus_pairs,
            ws.q_inj_vars, eps_abs, eps_rel,
        )
        eps_s_q_inj = _eps_dual_bus_pairs_signed(
            ws.vars.q_inj.duals, ws.bus_pairs, eps_abs, eps_rel,
        )
        ws.eps_primal = sqrt(
            eps_r_bus^2 + eps_r_branch^2 + eps_r_p_inj^2 + eps_r_q_inj^2,
        )
        ws.eps_dual = sqrt(
            eps_s_bus^2 + eps_s_branch^2 + eps_s_p_inj^2 + eps_s_q_inj^2,
        )
    else
        ws.eps_primal = hypot(eps_r_bus, eps_r_branch)
        ws.eps_dual = hypot(eps_s_bus, eps_s_branch)
    end
    return nothing
end

# ---- SOCWR ----

function compute_tolerances!(ws::SOCWRStateWorkspace, eps_abs::Float64, eps_rel::Float64)
    eps_r = hypot(
        _eps_primal_bus(
            ws.partitions,
            (ws.vars.v2.consensus,),
            (snapshot_v2!,),
            eps_abs,
            eps_rel,
        ),
        _eps_primal_branch(
            ws.partitions,
            (ws.vars.r.consensus, ws.vars.t.consensus),
            (:R, :T_var),
            ws.interzonal,
            eps_abs,
            eps_rel,
        ),
    )
    eps_s = hypot(
        _eps_dual_bus(
            (ws.vars.v2.duals,),
            ws.partitions,
            eps_abs,
            eps_rel,
        ),
        _eps_dual_branch(
            (ws.vars.r.duals, ws.vars.t.duals),
            ws.partitions,
            ws.interzonal,
            eps_abs,
            eps_rel,
        ),
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end

function compute_tolerances!(ws::SOCWRPowerWorkspace, eps_abs::Float64, eps_rel::Float64)
    n_inter = length(ws.interzonal)
    n_pairs = length(ws.branch_pairs)

    eps_r_inter =
        n_inter > 0 ?
        _eps_primal_branch(
            ws.partitions,
            (view(ws.vars.p_flow.consensus, 1:n_inter),
                view(ws.vars.q_flow.consensus, 1:n_inter)),
            (:p_flow_branches, :q_flow_branches),
            ws.interzonal,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_inter =
        n_inter > 0 ?
        _eps_dual_branch(
            (view(ws.vars.p_flow.duals, 1:n_inter, :),
                view(ws.vars.q_flow.duals, 1:n_inter, :)),
            ws.partitions,
            ws.interzonal,
            eps_abs,
            eps_rel,
        ) : 0.0

    eps_r_p_pair =
        n_pairs > 0 ?
        _eps_primal_branch_pairs(
            ws.partitions,
            view(ws.vars.p_flow.consensus, n_inter+1:n_inter+n_pairs),
            ws.branch_pairs,
            :p_flow_branches,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_r_q_pair =
        n_pairs > 0 ?
        _eps_primal_branch_pairs(
            ws.partitions,
            view(ws.vars.q_flow.consensus, n_inter+1:n_inter+n_pairs),
            ws.branch_pairs,
            :q_flow_branches,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_p_pair =
        n_pairs > 0 ?
        _eps_dual_branch_pairs(
            view(ws.vars.p_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.branch_pairs,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_q_pair =
        n_pairs > 0 ?
        _eps_dual_branch_pairs(
            view(ws.vars.q_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.branch_pairs,
            eps_abs,
            eps_rel,
        ) : 0.0

    ws.eps_primal = sqrt(eps_r_inter^2 + eps_r_p_pair^2 + eps_r_q_pair^2)
    ws.eps_dual = sqrt(eps_s_inter^2 + eps_s_p_pair^2 + eps_s_q_pair^2)
    return nothing
end

function compute_tolerances!(ws::SOCWRHybridWorkspace, eps_abs::Float64, eps_rel::Float64)
    n_inter = length(ws.interzonal)
    n_pairs = length(ws.branch_pairs)

    eps_r_bus = _eps_primal_bus(
        ws.partitions,
        (ws.vars.v2.consensus,),
        (snapshot_v2!,),
        eps_abs,
        eps_rel,
    )
    eps_s_bus = _eps_dual_bus(
        (ws.vars.v2.duals,),
        ws.partitions,
        eps_abs,
        eps_rel,
    )

    # Interzonal P/Q (TieLine)
    eps_r_pq_inter =
        n_inter > 0 ?
        _eps_primal_branch(
            ws.partitions,
            (view(ws.vars.p_flow.consensus, 1:n_inter),
                view(ws.vars.q_flow.consensus, 1:n_inter)),
            (:p_flow_branches, :q_flow_branches),
            ws.interzonal,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_pq_inter =
        n_inter > 0 ?
        _eps_dual_branch(
            (view(ws.vars.p_flow.duals, 1:n_inter, :),
                view(ws.vars.q_flow.duals, 1:n_inter, :)),
            ws.partitions,
            ws.interzonal,
            eps_abs,
            eps_rel,
        ) : 0.0

    # Branch pairs P/Q (DummyBus)
    eps_r_p_pair =
        n_pairs > 0 ?
        _eps_primal_branch_pairs(
            ws.partitions,
            view(ws.vars.p_flow.consensus, n_inter+1:n_inter+n_pairs),
            ws.branch_pairs,
            :p_flow_branches,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_r_q_pair =
        n_pairs > 0 ?
        _eps_primal_branch_pairs(
            ws.partitions,
            view(ws.vars.q_flow.consensus, n_inter+1:n_inter+n_pairs),
            ws.branch_pairs,
            :q_flow_branches,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_p_pair =
        n_pairs > 0 ?
        _eps_dual_branch_pairs(
            view(ws.vars.p_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.branch_pairs,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_q_pair =
        n_pairs > 0 ?
        _eps_dual_branch_pairs(
            view(ws.vars.q_flow.duals, n_inter+1:n_inter+n_pairs, :),
            ws.branch_pairs,
            eps_abs,
            eps_rel,
        ) : 0.0

    eps_r_pq = sqrt(eps_r_pq_inter^2 + eps_r_p_pair^2 + eps_r_q_pair^2)
    eps_s_pq = sqrt(eps_s_pq_inter^2 + eps_s_p_pair^2 + eps_s_q_pair^2)

    eps_r_rt =
        n_inter > 0 ?
        _eps_primal_branch(
            ws.partitions,
            (ws.vars.r.consensus, ws.vars.t.consensus),
            (:R, :T_var),
            ws.interzonal,
            eps_abs,
            eps_rel,
        ) : 0.0
    eps_s_rt =
        n_inter > 0 ?
        _eps_dual_branch(
            (ws.vars.r.duals, ws.vars.t.duals),
            ws.partitions,
            ws.interzonal,
            eps_abs,
            eps_rel,
        ) : 0.0

    if !isempty(ws.bus_pairs)
        eps_r_p_inj = _eps_primal_bus_pairs_signed(
            ws.partitions, ws.vars.p_inj.consensus, ws.bus_pairs,
            ws.p_inj_vars, eps_abs, eps_rel,
        )
        eps_s_p_inj = _eps_dual_bus_pairs_signed(
            ws.vars.p_inj.duals, ws.bus_pairs, eps_abs, eps_rel,
        )
        eps_r_q_inj = _eps_primal_bus_pairs_signed(
            ws.partitions, ws.vars.q_inj.consensus, ws.bus_pairs,
            ws.q_inj_vars, eps_abs, eps_rel,
        )
        eps_s_q_inj = _eps_dual_bus_pairs_signed(
            ws.vars.q_inj.duals, ws.bus_pairs, eps_abs, eps_rel,
        )
        ws.eps_primal = sqrt(
            eps_r_bus^2 + eps_r_pq^2 + eps_r_rt^2 +
            eps_r_p_inj^2 + eps_r_q_inj^2,
        )
        ws.eps_dual = sqrt(
            eps_s_bus^2 + eps_s_pq^2 + eps_s_rt^2 +
            eps_s_p_inj^2 + eps_s_q_inj^2,
        )
    else
        ws.eps_primal = sqrt(eps_r_bus^2 + eps_r_pq^2 + eps_r_rt^2)
        ws.eps_dual = sqrt(eps_s_bus^2 + eps_s_pq^2 + eps_s_rt^2)
    end
    return nothing
end
