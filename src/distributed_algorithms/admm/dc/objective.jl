# ============================================================
# Quadratic terms
# ============================================================

# ------------------------------------------------------------
# State
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

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function build_quad_branch_expr!(
    part::PartitionDC,
    rho_p::Float64,
    zone_edges::Vector{Int},
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    m = part.model::JuMP.Model
    quad = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue

        p = m[:p_flow_branches][l, 1]
        quad += 0.5 * rho_p * p^2
    end

    pos = part.position
    @inbounds for (ga, za, gb, zb) in branch_pairs
        if za == pos
            l = part.branch_idx_full_to_zone[ga]
            if l != 0
                p = m[:p_flow_branches][l, 1]
                quad += 0.5 * rho_p * p^2
            end
        end
        if zb == pos
            l = part.branch_idx_full_to_zone[gb]
            if l != 0
                p = m[:p_flow_branches][l, 1]
                quad += 0.5 * rho_p * p^2
            end
        end
    end

    part.quad_branch_expr = quad
    return quad
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function build_quad_p_inj_expr!(
    part::PartitionDC,
    rho::Float64,
    p_inj_store::Union{Nothing, Dict{Int, JuMP.VariableRef}},
)
    quad = zero(JuMP.QuadExpr)
    if p_inj_store !== nothing
        for (_, v) in p_inj_store
            quad += 0.5 * rho * v^2
        end
    end
    part.quad_p_inj_expr = quad
    return quad
end

# ============================================================
# Linear terms
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function build_linear_bus_term(
    part::PartitionDC,
    theta_consensus::Vector{Float64},
    theta_duals::AbstractVector{Float64},
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

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function build_linear_branch_term(
    part::PartitionDC,
    zone_edges::Vector{Int},
    p_flow_consensus::Vector{Float64},
    duals_p_flow::AbstractVector{Float64},
    interzonal::Vector{Int},
    rho_p::Float64,
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)::JuMP.AffExpr
    m = part.model::JuMP.Model
    lin = zero(JuMP.AffExpr)
    n_inter = length(interzonal)
    @inbounds for (k, e) in enumerate(interzonal)
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue
        p = m[:p_flow_branches][l, 1]
        p_cons = p_flow_consensus[k]
        p_dual = duals_p_flow[k]
        lin += (p_dual - rho_p * p_cons) * p
    end

    pos = part.position
    @inbounds for (k, (ga, za, gb, zb)) in enumerate(branch_pairs)
        idx = n_inter + k
        p_cons = p_flow_consensus[idx]
        p_dual = duals_p_flow[idx]
        if za == pos
            l = part.branch_idx_full_to_zone[ga]
            if l != 0
                p = m[:p_flow_branches][l, 1]
                lin += (p_dual - rho_p * p_cons) * p
            end
        end
        if zb == pos
            l = part.branch_idx_full_to_zone[gb]
            if l != 0
                p = m[:p_flow_branches][l, 1]
                lin += (p_dual - rho_p * p_cons) * p
            end
        end
    end
    return lin
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function build_linear_p_inj_term(
    part::PartitionDC,
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
    p_inj_consensus::Vector{Float64},
    p_inj_duals::AbstractVector{Float64},
    rho::Float64,
)::JuMP.AffExpr
    lin = zero(JuMP.AffExpr)
    isempty(bus_pairs) && return lin

    pos = part.position
    store = p_inj_vars[pos]
    isempty(store) && return lin

    @inbounds for (k, (bid, za, zb)) in enumerate(bus_pairs)
        z = p_inj_consensus[k]
        lam = p_inj_duals[k]
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

# ============================================================
# Update objective
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function update_objective!(part::PartitionDC, ws::DCStateWorkspace)
    lin = build_linear_bus_term(
        part,
        ws.vars.theta.consensus,
        view(ws.vars.theta.duals, :, part.position),
        ws.vars.theta.rho[],
    )

    @objective(part.model, Min, part.obj_static + lin)
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function update_objective!(part::PartitionDC, ws::DCPowerWorkspace)
    lin = build_linear_branch_term(
        part,
        part.zone_edges,
        ws.vars.p_flow.consensus,
        view(ws.vars.p_flow.duals, :, part.position),
        ws.interzonal,
        ws.vars.p_flow.rho[],
        ws.branch_pairs,
    )

    @objective(part.model, Min, part.obj_static + lin)
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function update_objective!(
    part::PartitionDC,
    ws::DCHybridWorkspace{<:Union{TieLine, DummyBus}},
)
    lin_bus = build_linear_bus_term(
        part,
        ws.vars.theta.consensus,
        view(ws.vars.theta.duals, :, part.position),
        ws.vars.theta.rho[],
    )

    lin_branch = build_linear_branch_term(
        part,
        part.zone_edges,
        ws.vars.p_flow.consensus,
        view(ws.vars.p_flow.duals, :, part.position),
        ws.interzonal,
        ws.vars.p_flow.rho[],
        ws.branch_pairs,
    )

    @objective(part.model, Min, part.obj_static + lin_bus + lin_branch)
    return nothing
end

function update_objective!(
    part::PartitionDC,
    ws::DCHybridWorkspace{SplitBus},
)
    lin_bus = build_linear_bus_term(
        part,
        ws.vars.theta.consensus,
        view(ws.vars.theta.duals, :, part.position),
        ws.vars.theta.rho[],
    )

    lin_branch = build_linear_branch_term(
        part,
        part.zone_edges,
        ws.vars.p_flow.consensus,
        view(ws.vars.p_flow.duals, :, part.position),
        ws.interzonal,
        ws.vars.p_flow.rho[],
        ws.branch_pairs,
    )

    lin_p_inj = build_linear_p_inj_term(
        part,
        ws.bus_pairs,
        ws.p_inj_vars,
        ws.vars.p_inj.consensus,
        view(ws.vars.p_inj.duals, :, part.position),
        ws.vars.p_inj.rho[],
    )

    @objective(
        part.model, Min,
        part.obj_static + lin_bus + lin_branch + lin_p_inj,
    )
    return nothing
end
