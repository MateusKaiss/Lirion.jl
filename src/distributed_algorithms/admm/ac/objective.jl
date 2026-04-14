# ============================================================
# AC OBJECTIVE
# ============================================================

# ------------------------------------------------------------
# QUADRATIC TERMS
# ------------------------------------------------------------

function build_quad_bus_expr!(
    part::PartitionAC,
    rho_theta::Float64,
    rho_v_mag::Float64,
)
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    v_mag = part.v_mag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    quad = zero(JuMP.QuadExpr)

    @inbounds for b in part.zone_shared_bus_ids
        quad += 0.5 * rho_theta * theta[b, 1]^2 +
                0.5 * rho_v_mag * v_mag[b, 1]^2
    end

    part.quad_bus_expr = quad
    return quad
end

function build_quad_branch_expr!(
    part::PartitionAC,
    rho_p::Float64,
    rho_q::Float64,
    zone_edges::Vector{Int},
)
    quad = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue

        p, q = _pq_flow_from_exprs(part, l)
        quad += 0.5 * rho_p * p^2 + 0.5 * rho_q * q^2
    end

    part.quad_branch_expr = quad
    return quad
end

# ------------------------------------------------------------
# LINEAR TERMS
# ------------------------------------------------------------

function build_linear_bus_term(
    part::PartitionAC,
    theta_consensus::Vector{Float64},
    v_mag_consensus::Vector{Float64},
    theta_duals::AbstractVector{Float64},
    v_mag_duals::AbstractVector{Float64},
    rho_theta::Float64,
    rho_v_mag::Float64,
)::JuMP.AffExpr
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    v_mag = part.v_mag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)

    @inbounds for i in eachindex(part.zone_shared_bus_ids)
        b = part.zone_shared_bus_ids[i]
        j = part.shared_bus_consensus_idx[i]
        lin += (theta_duals[j] - rho_theta * theta_consensus[j]) * theta[b, 1]
        lin += (v_mag_duals[j] - rho_v_mag * v_mag_consensus[j]) * v_mag[b, 1]
    end

    return lin
end

function build_linear_branch_term(
    part::PartitionAC,
    zone_edges::Vector{Int},
    p_consensus::Vector{Float64},
    q_consensus::Vector{Float64},
    duals_p::AbstractVector{Float64},
    duals_q::AbstractVector{Float64},
    interzonal::Vector{Int},
    rho_p::Float64,
    rho_q::Float64,
)
    lin = zero(JuMP.AffExpr)

    @inbounds for (k, e) in enumerate(interzonal)
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue

        p, q = _pq_flow_from_exprs(part, l)

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

function update_objective!(part::PartitionAC, ws::ACStateWorkspace)
    lin = build_linear_bus_term(
        part,
        ws.theta_consensus,
        ws.v_mag_consensus,
        view(ws.duals_theta, :, part.position),
        view(ws.duals_v_mag, :, part.position),
        ws.rho_theta[],
        ws.rho_v[],
    )

    @objective(part.model, Min, part.obj_base + part.quad_bus_expr + lin)
end

function update_objective!(part::PartitionAC, ws::ACPowerWorkspace)
    lin = build_linear_branch_term(
        part,
        part.zone_edges,
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        view(ws.duals_p_flow, :, part.position),
        view(ws.duals_q_flow, :, part.position),
        ws.interzonal,
        ws.rho_p[],
        ws.rho_q[],
    )

    @objective(part.model, Min, part.obj_base + part.quad_branch_expr + lin)
end

function update_objective!(part::PartitionAC, ws::ACHybridWorkspace)
    lin_bus = build_linear_bus_term(
        part,
        ws.theta_consensus,
        ws.v_mag_consensus,
        view(ws.duals_theta, :, part.position),
        view(ws.duals_v_mag, :, part.position),
        ws.rho_theta[],
        ws.rho_v[],
    )

    lin_branch = build_linear_branch_term(
        part,
        part.zone_edges,
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        view(ws.duals_p_flow, :, part.position),
        view(ws.duals_q_flow, :, part.position),
        ws.interzonal,
        ws.rho_p[],
        ws.rho_q[],
    )

    @objective(
        part.model,
        Min,
        part.obj_base + part.quad_bus_expr + part.quad_branch_expr +
        lin_bus + lin_branch,
    )
end
