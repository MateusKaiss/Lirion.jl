# ============================================================
# SOCP OBJECTIVE
# ============================================================

# ------------------------------------------------------------
# QUADRATIC TERMS
# ------------------------------------------------------------

function build_quad_bus_expr!(part::PartitionSOCP, rho_v2::Float64)
    v2 = part.v2_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    quad = zero(JuMP.QuadExpr)

    @inbounds for b in part.zone_shared_bus_ids
        quad += 0.5 * rho_v2 * v2[b, 1]^2
    end

    part.quad_bus_expr_socp = quad
    return quad
end

function build_quad_branch_rt_expr!(
    part::PartitionSOCP,
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

    part.quad_rt_expr_socp = quad
    return quad
end

function build_quad_branch_pq_expr!(
    part::PartitionSOCP,
    rho_p::Float64,
    rho_q::Float64,
    zone_edges::Vector{Int},
)
    m = part.model::JuMP.Model
    quad = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue
        p, q = _pq_flow_from_exprs(m, l)
        quad += 0.5 * rho_p * p^2 + 0.5 * rho_q * q^2
    end

    part.quad_pq_expr_socp = quad
    return quad
end

# ------------------------------------------------------------
# LINEAR TERMS
# ------------------------------------------------------------

function build_linear_bus_term(
    part::PartitionSOCP,
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
    part::PartitionSOCP,
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
    part::PartitionSOCP,
    zone_edges::Vector{Int},
    p_consensus::Vector{Float64},
    q_consensus::Vector{Float64},
    duals_p::AbstractVector{Float64},
    duals_q::AbstractVector{Float64},
    interzonal::Vector{Int},
    rho_p::Float64,
    rho_q::Float64,
)::JuMP.AffExpr
    m = part.model::JuMP.Model
    lin = zero(JuMP.AffExpr)

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

    return lin
end

# ------------------------------------------------------------
# UPDATE OBJECTIVE
# ------------------------------------------------------------

function update_objective!(part::PartitionSOCP, ws::SOCPStateWorkspace)
    p = part.position

    lin_bus = build_linear_bus_term(
        part,
        ws.v2_consensus,
        view(ws.duals_v2, :, p),
        ws.rho_v2[],
    )

    lin_rt = build_linear_branch_rt_term(
        part,
        part.zone_edges,
        ws.r_consensus,
        ws.t_consensus,
        view(ws.duals_r, :, p),
        view(ws.duals_t, :, p),
        ws.interzonal,
        ws.rho_r[],
        ws.rho_t[],
    )

    @objective(part.model, Min, part.obj_static + lin_bus + lin_rt)
end

function update_objective!(part::PartitionSOCP, ws::SOCPPowerWorkspace)
    p = part.position

    lin = build_linear_branch_pq_term(
        part,
        part.zone_edges,
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        view(ws.duals_p, :, p),
        view(ws.duals_q, :, p),
        ws.interzonal,
        ws.rho_p[],
        ws.rho_q[],
    )

    @objective(part.model, Min, part.obj_static + lin)
end

function update_objective!(part::PartitionSOCP, ws::SOCPHybridWorkspace)
    p = part.position

    lin_bus = build_linear_bus_term(
        part,
        ws.v2_consensus,
        view(ws.duals_v2, :, p),
        ws.rho_v2[],
    )

    lin_pq = build_linear_branch_pq_term(
        part,
        part.zone_edges,
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        view(ws.duals_p, :, p),
        view(ws.duals_q, :, p),
        ws.interzonal,
        ws.rho_p[],
        ws.rho_q[],
    )

    lin_rt = build_linear_branch_rt_term(
        part,
        part.zone_edges,
        ws.r_consensus,
        ws.t_consensus,
        view(ws.duals_r, :, p),
        view(ws.duals_t, :, p),
        ws.interzonal,
        ws.rho_r[],
        ws.rho_t[],
    )

    @objective(part.model, Min, part.obj_static + lin_bus + lin_pq + lin_rt)
end
