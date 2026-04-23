# ============================================================
# SOCWR OBJECTIVE
# ============================================================

# ------------------------------------------------------------
# QUADRATIC TERMS
# ------------------------------------------------------------

function build_quad_bus_expr!(part::PartitionSOCWR, rho_v2::Float64)
    v2 = part.v2_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    quad = zero(JuMP.QuadExpr)

    @inbounds for b in part.zone_shared_bus_ids
        quad += 0.5 * rho_v2 * v2[b, 1]^2
    end

    part.quad_bus_expr = quad
    return quad
end

function build_quad_branch_rt_expr!(
    part::PartitionSOCWR,
    rho_r::Float64,
    rho_t::Float64,
    zone_edges::Vector{Int},
)
    m = part.model::JuMP.Model
    quad = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue
        r, t = _rt_from_exprs(m, l)
        quad += 0.5 * rho_r * r^2 + 0.5 * rho_t * t^2
    end

    part.quad_rt_expr = quad
    return quad
end

function build_quad_branch_pq_expr!(
    part::PartitionSOCWR,
    rho_p::Float64,
    rho_q::Float64,
    zone_edges::Vector{Int},
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    m = part.model::JuMP.Model
    quad = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue
        p, q = _pq_flow_from_exprs(m, l)
        quad += 0.5 * rho_p * p^2 + 0.5 * rho_q * q^2
    end
    pos = part.position
    @inbounds for (ga, za, gb, zb) in branch_pairs
        if za == pos
            l = part.branch_idx_full_to_zone[ga]
            if l != 0
                p, q = _pq_flow_from_exprs(m, l)
                quad += 0.5 * rho_p * p^2 + 0.5 * rho_q * q^2
            end
        end
        if zb == pos
            l = part.branch_idx_full_to_zone[gb]
            if l != 0
                p, q = _pq_flow_from_exprs(m, l)
                quad += 0.5 * rho_p * p^2 + 0.5 * rho_q * q^2
            end
        end
    end

    part.quad_pq_expr = quad
    return quad
end

# ------------------------------------------------------------
# LINEAR TERMS
# ------------------------------------------------------------

function build_linear_bus_term(
    part::PartitionSOCWR,
    v2_consensus::Vector{Float64},
    v2_duals::AbstractVector{Float64},
    rho_v2::Float64,
)::JuMP.AffExpr
    v2 = part.v2_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)

    @inbounds for i in eachindex(part.zone_shared_bus_ids)
        b = part.zone_shared_bus_ids[i]
        j = part.shared_bus_consensus_idx[i]
        lin += (v2_duals[j] - rho_v2 * v2_consensus[j]) * v2[b, 1]
    end

    return lin
end

function build_linear_branch_rt_term(
    part::PartitionSOCWR,
    zone_edges::Vector{Int},
    r_consensus::Vector{Float64},
    t_consensus::Vector{Float64},
    duals_r::AbstractVector{Float64},
    duals_t::AbstractVector{Float64},
    interzonal::Vector{Int},
    rho_r::Float64,
    rho_t::Float64,
)::JuMP.AffExpr
    m = part.model::JuMP.Model
    lin = zero(JuMP.AffExpr)

    @inbounds for (k, e) in enumerate(interzonal)
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue
        r, t = _rt_from_exprs(m, l)

        r_cons = r_consensus[k]
        t_cons = t_consensus[k]
        r_dual = duals_r[k]
        t_dual = duals_t[k]

        lin += (r_dual - rho_r * r_cons) * r + (t_dual - rho_t * t_cons) * t
    end

    return lin
end

function build_linear_branch_pq_term(
    part::PartitionSOCWR,
    zone_edges::Vector{Int},
    p_consensus::Vector{Float64},
    q_consensus::Vector{Float64},
    duals_p::AbstractVector{Float64},
    duals_q::AbstractVector{Float64},
    interzonal::Vector{Int},
    rho_p::Float64,
    rho_q::Float64,
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)::JuMP.AffExpr
    m = part.model::JuMP.Model
    lin = zero(JuMP.AffExpr)
    n_inter = length(interzonal)

    # Interzonal (TieLine): indices 1..n_inter in consensus/duals.
    @inbounds for (k, e) in enumerate(interzonal)
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue
        p, q = _pq_flow_from_exprs(m, l)

        p_cons = p_consensus[k]
        q_cons = q_consensus[k]
        p_dual = duals_p[k]
        q_dual = duals_q[k]

        lin += (p_dual - rho_p * p_cons) * p + (q_dual - rho_q * q_cons) * q
    end

    # Branch pairs (DummyBus): indices n_inter+1..n_inter+n_pairs.
    pos = part.position
    @inbounds for (k, (ga, za, gb, zb)) in enumerate(branch_pairs)
        idx = n_inter + k
        p_cons = p_consensus[idx]
        q_cons = q_consensus[idx]
        p_dual = duals_p[idx]
        q_dual = duals_q[idx]
        if za == pos
            l = part.branch_idx_full_to_zone[ga]
            if l != 0
                p, q = _pq_flow_from_exprs(m, l)
                lin += (p_dual - rho_p * p_cons) * p +
                       (q_dual - rho_q * q_cons) * q
            end
        end
        if zb == pos
            l = part.branch_idx_full_to_zone[gb]
            if l != 0
                p, q = _pq_flow_from_exprs(m, l)
                lin += (p_dual - rho_p * p_cons) * p +
                       (q_dual - rho_q * q_cons) * q
            end
        end
    end

    return lin
end

# ------------------------------------------------------------
# UPDATE OBJECTIVE
# ------------------------------------------------------------

function update_objective!(part::PartitionSOCWR, ws::SOCWRStateWorkspace)
    p = part.position

    lin_bus = build_linear_bus_term(
        part,
        ws.vars.v2.consensus,
        view(ws.vars.v2.duals, :, p),
        ws.vars.v2.rho[],
    )

    lin_rt = build_linear_branch_rt_term(
        part,
        part.zone_edges,
        ws.vars.r.consensus,
        ws.vars.t.consensus,
        view(ws.vars.r.duals, :, p),
        view(ws.vars.t.duals, :, p),
        ws.interzonal,
        ws.vars.r.rho[],
        ws.vars.t.rho[],
    )

    @objective(part.model, Min, part.obj_static + lin_bus + lin_rt)
end

function update_objective!(part::PartitionSOCWR, ws::SOCWRPowerWorkspace)
    p = part.position

    lin = build_linear_branch_pq_term(
        part,
        part.zone_edges,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        view(ws.vars.p_flow.duals, :, p),
        view(ws.vars.q_flow.duals, :, p),
        ws.interzonal,
        ws.vars.p_flow.rho[],
        ws.vars.q_flow.rho[],
        ws.branch_pairs,
    )

    @objective(part.model, Min, part.obj_static + lin)
end

function update_objective!(part::PartitionSOCWR, ws::SOCWRHybridWorkspace)
    p = part.position

    lin_bus = build_linear_bus_term(
        part,
        ws.vars.v2.consensus,
        view(ws.vars.v2.duals, :, p),
        ws.vars.v2.rho[],
    )

    lin_pq = build_linear_branch_pq_term(
        part,
        part.zone_edges,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        view(ws.vars.p_flow.duals, :, p),
        view(ws.vars.q_flow.duals, :, p),
        ws.interzonal,
        ws.vars.p_flow.rho[],
        ws.vars.q_flow.rho[],
        ws.branch_pairs,
    )

    lin_rt = build_linear_branch_rt_term(
        part,
        part.zone_edges,
        ws.vars.r.consensus,
        ws.vars.t.consensus,
        view(ws.vars.r.duals, :, p),
        view(ws.vars.t.duals, :, p),
        ws.interzonal,
        ws.vars.r.rho[],
        ws.vars.t.rho[],
    )

    # Fast path: TieLine has no split buses, skip p_inj/q_inj machinery.
    if isempty(ws.bus_pairs)
        @objective(part.model, Min, part.obj_static + lin_bus + lin_pq + lin_rt)
        return nothing
    end

    lin_p_inj = build_linear_inj_term_socwr(
        part,
        ws.bus_pairs,
        ws.p_inj_vars,
        ws.vars.p_inj.consensus,
        view(ws.vars.p_inj.duals, :, p),
        ws.vars.p_inj.rho[],
    )
    lin_q_inj = build_linear_inj_term_socwr(
        part,
        ws.bus_pairs,
        ws.q_inj_vars,
        ws.vars.q_inj.consensus,
        view(ws.vars.q_inj.duals, :, p),
        ws.vars.q_inj.rho[],
    )

    @objective(
        part.model, Min,
        part.obj_static + lin_bus + lin_pq + lin_rt + lin_p_inj + lin_q_inj,
    )
    return nothing
end

# ------------------------------------------------------------
# Signed linear term for virtual injection at split buses
# (p_inj or q_inj). Zone A targets +z: (λ - ρz)·v.
# Zone B targets -z:                   (λ + ρz)·v.
# ------------------------------------------------------------
function build_linear_inj_term_socwr(
    part::PartitionSOCWR,
    bus_pairs::Vector{NTuple{3, Int}},
    inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
    consensus::Vector{Float64},
    duals::AbstractVector{Float64},
    rho::Float64,
)::JuMP.AffExpr
    lin = zero(JuMP.AffExpr)
    isempty(bus_pairs) && return lin
    pos = part.position
    store = inj_vars[pos]
    isempty(store) && return lin
    @inbounds for (k, (bid, za, zb)) in enumerate(bus_pairs)
        z = consensus[k]
        lam = duals[k]
        if za == pos
            v = get(store, bid, nothing)
            v === nothing && continue
            lin += (lam - rho * z) * v
        elseif zb == pos
            v = get(store, bid, nothing)
            v === nothing && continue
            lin += (lam + rho * z) * v
        end
    end
    return lin
end
