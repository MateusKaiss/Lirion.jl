# ============================================================
# DC OBJECTIVE
# ============================================================

# ------------------------------------------------------------
# QUADRATIC TERMS
# ------------------------------------------------------------

function build_quad_bus_expr!(
    part::PartitionDC,
    rho_theta::Float64,
)
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    quad = zero(JuMP.QuadExpr)

    @inbounds for b in part.zone_shared_bus_ids
        quad += 0.5 * rho_theta * theta[b, 1]^2
    end

    part.quad_bus_expr = quad
    return quad
end

function build_quad_branch_expr!(
    part::PartitionDC,
    rho_p::Float64,
    zone_edges::Vector{Int},
)
    m = part.model::JuMP.Model
    quad = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue

        p = m[:p_flow_branches][l, 1]
        quad += 0.5 * rho_p * p^2
    end

    part.quad_branch_expr = quad
    return quad
end

# ------------------------------------------------------------
# LINEAR TERMS
# ------------------------------------------------------------

function build_linear_bus_term(
    part::PartitionDC,
    theta_consensus::Vector{Float64},
    theta_duals::AbstractVector{Float64},   # column view over n_shared
    rho_theta::Float64,
)::JuMP.AffExpr
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)

    @inbounds for i in eachindex(part.zone_shared_bus_ids)
        b = part.zone_shared_bus_ids[i]
        j = part.shared_bus_consensus_idx[i]
        lin += (theta_duals[j] - rho_theta * theta_consensus[j]) * theta[b, 1]
    end
    return lin
end

function build_linear_branch_term(
    part::PartitionDC,
    zone_edges::Vector{Int},
    p_flow_consensus::Vector{Float64},
    duals_p_flow::AbstractVector{Float64},   # column view over n_inter
    interzonal::Vector{Int},
    rho_p::Float64,
)::JuMP.AffExpr
    m = part.model::JuMP.Model
    lin = zero(JuMP.AffExpr)

    # Build a lookup from full-branch index to position in `interzonal`.
    # We iterate enumerate(interzonal) and check if this zone cuts it.
    @inbounds for (k, e) in enumerate(interzonal)
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue
        p = m[:p_flow_branches][l, 1]
        p_cons = p_flow_consensus[k]
        p_dual = duals_p_flow[k]
        lin += (p_dual - rho_p * p_cons) * p
    end
    return lin
end

# ------------------------------------------------------------
# UPDATE OBJECTIVE
# ------------------------------------------------------------

function update_objective!(part::PartitionDC, ws::DCStateWorkspace)
    lin = build_linear_bus_term(
        part,
        ws.theta_consensus,
        view(ws.duals_theta, :, part.position),
        ws.rho_theta[],
    )

    @objective(
        part.model,
        Min,
        part.obj_base + part.quad_bus_expr + lin,
    )
end

function update_objective!(part::PartitionDC, ws::DCPowerWorkspace)
    lin = build_linear_branch_term(
        part,
        part.zone_edges,
        ws.p_flow_consensus,
        view(ws.duals_p_flow, :, part.position),
        ws.interzonal,
        ws.rho_p[],
    )

    @objective(
        part.model,
        Min,
        part.obj_base + part.quad_branch_expr + lin,
    )
end

function update_objective!(part::PartitionDC, ws::DCHybridWorkspace)
    lin_bus = build_linear_bus_term(
        part,
        ws.theta_consensus,
        view(ws.duals_theta, :, part.position),
        ws.rho_theta[],
    )

    lin_branch = build_linear_branch_term(
        part,
        part.zone_edges,
        ws.p_flow_consensus,
        view(ws.duals_p_flow, :, part.position),
        ws.interzonal,
        ws.rho_p[],
    )

    @objective(
        part.model,
        Min,
        part.obj_base +
        part.quad_bus_expr +
        part.quad_branch_expr +
        lin_bus +
        lin_branch,
    )
end
