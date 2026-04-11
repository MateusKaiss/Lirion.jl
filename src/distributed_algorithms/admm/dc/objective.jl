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
    consensus_bus_idx_map::Dict{Int, Int},
    theta_consensus::Vector{Float64},
    theta_duals::Vector{Float64},
    rho_theta::Float64,
)::JuMP.AffExpr
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)

    @inbounds for b in part.zone_shared_bus_ids
        j = consensus_bus_idx_map[b]
        lin += (theta_duals[j] - rho_theta * theta_consensus[j]) * theta[b, 1]
    end
    return lin
end

function build_linear_branch_term(
    part::PartitionDC,
    zone_edges::Vector{Int},
    p_flow_consensus::Dict{Int, Float64},
    duals_p_flow::Dict{Int, Float64},
    rho_p::Float64,
)::JuMP.AffExpr
    m = part.model::JuMP.Model
    lin = zero(JuMP.AffExpr)

    @inbounds for e in zone_edges
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue

        p = m[:p_flow_branches][l, 1]
        p_cons = get(p_flow_consensus, e, 0.0)
        p_dual = get(duals_p_flow, e, 0.0)

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
        ws.consensus_bus_idx_map,
        ws.theta_consensus,
        ws.duals_theta[part.id],
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
        ws.duals_p_flow[part.id],
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
        ws.consensus_bus_idx_map,
        ws.theta_consensus,
        ws.duals_theta[part.id],
        ws.rho_theta[],
    )

    lin_branch = build_linear_branch_term(
        part,
        part.zone_edges,
        ws.p_flow_consensus,
        ws.duals_p_flow[part.id],
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
