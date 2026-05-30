# ============================================================
# Quadratic terms
# ============================================================

# ------------------------------------------------------------
# State
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

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function build_quad_branch_expr!(
    part::PartitionAC,
    rho_p::Float64,
    rho_q::Float64,
    zone_edges::Vector{Int},
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    quad = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue

        p, q = _pq_flow_from_exprs(part, l)
        quad += 0.5 * rho_p * p^2 + 0.5 * rho_q * q^2
    end

    pos = part.position
    @inbounds for (ga, za, gb, zb) in branch_pairs
        if za == pos
            l = part.branch_idx_full_to_zone[ga]
            if l != 0
                p, q = _pq_flow_from_exprs(part, l)
                quad += 0.5 * rho_p * p^2 + 0.5 * rho_q * q^2
            end
        end
        if zb == pos
            l = part.branch_idx_full_to_zone[gb]
            if l != 0
                p, q = _pq_flow_from_exprs(part, l)
                quad += 0.5 * rho_p * p^2 + 0.5 * rho_q * q^2
            end
        end
    end

    part.quad_branch_expr = quad
    return quad
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function build_quad_p_inj_expr_ac!(
    part::PartitionAC,
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

function build_quad_q_inj_expr_ac!(
    part::PartitionAC,
    rho::Float64,
    q_inj_store::Union{Nothing, Dict{Int, JuMP.VariableRef}},
)
    quad = zero(JuMP.QuadExpr)
    if q_inj_store !== nothing
        for (_, v) in q_inj_store
            quad += 0.5 * rho * v^2
        end
    end
    part.quad_q_inj_expr = quad
    return quad
end

# ============================================================
# Linear terms
# ============================================================

# ------------------------------------------------------------
# State
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

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function build_linear_branch_term(
    part::PartitionAC,
    p_consensus::AbstractVector{Float64},
    q_consensus::AbstractVector{Float64},
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

function build_linear_dummy_term(
    part::PartitionAC,
    p_consensus::AbstractVector{Float64},
    q_consensus::AbstractVector{Float64},
    duals_p::AbstractVector{Float64},
    duals_q::AbstractVector{Float64},
    branch_pairs::Vector{NTuple{4, Int}},
    rho_p::Float64,
    rho_q::Float64,
)
    lin = zero(JuMP.AffExpr)
    isempty(branch_pairs) && return lin

    pos = part.position
    @inbounds for (k, (ga, za, gb, zb)) in enumerate(branch_pairs)
        p_cons = p_consensus[k]
        q_cons = q_consensus[k]
        p_dual = duals_p[k]
        q_dual = duals_q[k]
        if za == pos
            l = part.branch_idx_full_to_zone[ga]
            if l != 0
                p, q = _pq_flow_from_exprs(part, l)
                lin += (p_dual - rho_p * p_cons) * p +
                       (q_dual - rho_q * q_cons) * q
            end
        end
        if zb == pos
            l = part.branch_idx_full_to_zone[gb]
            if l != 0
                p, q = _pq_flow_from_exprs(part, l)
                lin += (p_dual - rho_p * p_cons) * p +
                       (q_dual - rho_q * q_cons) * q
            end
        end
    end

    return lin
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function build_linear_inj_term_ac(
    part::PartitionAC,
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
    @inbounds for (k, (b_id, za, zb)) in enumerate(bus_pairs)
        z = consensus[k]
        lam = duals[k]
        if za == pos
            v = get(store, b_id, nothing)
            v === nothing && continue
            lin += (lam - rho * z) * v
        elseif zb == pos
            v = get(store, b_id, nothing)
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

function update_objective!(part::PartitionAC, ws::ACStateWorkspace)
    lin = build_linear_bus_term(
        part,
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
        view(ws.vars.theta.duals, :, part.position),
        view(ws.vars.v_mag.duals, :, part.position),
        ws.vars.theta.rho[],
        ws.vars.v_mag.rho[],
    )

    @objective(part.model, Min, part.obj_base + part.quad_bus_expr + lin)
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function update_objective!(part::PartitionAC, ws::ACPowerWorkspace)
    lin_branch = build_linear_branch_term(
        part,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        view(ws.vars.p_flow.duals, :, part.position),
        view(ws.vars.q_flow.duals, :, part.position),
        ws.interzonal,
        ws.vars.p_flow.rho[],
        ws.vars.q_flow.rho[],
    )
    lin_dummy = build_linear_dummy_term(
        part,
        ws.vars.p_dummy.consensus,
        ws.vars.q_dummy.consensus,
        view(ws.vars.p_dummy.duals, :, part.position),
        view(ws.vars.q_dummy.duals, :, part.position),
        ws.branch_pairs,
        ws.vars.p_dummy.rho[],
        ws.vars.q_dummy.rho[],
    )

    @objective(
        part.model,
        Min,
        part.obj_base + part.quad_branch_expr + lin_branch + lin_dummy,
    )
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function update_objective!(
    part::PartitionAC,
    ws::ACHybridWorkspace{<:Union{TieLine, DummyBus}},
)
    lin_bus = build_linear_bus_term(
        part,
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
        view(ws.vars.theta.duals, :, part.position),
        view(ws.vars.v_mag.duals, :, part.position),
        ws.vars.theta.rho[],
        ws.vars.v_mag.rho[],
    )

    lin_branch = build_linear_branch_term(
        part,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        view(ws.vars.p_flow.duals, :, part.position),
        view(ws.vars.q_flow.duals, :, part.position),
        ws.interzonal,
        ws.vars.p_flow.rho[],
        ws.vars.q_flow.rho[],
    )
    lin_dummy = build_linear_dummy_term(
        part,
        ws.vars.p_dummy.consensus,
        ws.vars.q_dummy.consensus,
        view(ws.vars.p_dummy.duals, :, part.position),
        view(ws.vars.q_dummy.duals, :, part.position),
        ws.branch_pairs,
        ws.vars.p_dummy.rho[],
        ws.vars.q_dummy.rho[],
    )

    @objective(
        part.model, Min,
        part.obj_base + part.quad_bus_expr + part.quad_branch_expr +
        lin_bus + lin_branch + lin_dummy,
    )
    return nothing
end

function update_objective!(
    part::PartitionAC,
    ws::ACHybridWorkspace{SplitBus},
)
    lin_bus = build_linear_bus_term(
        part,
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
        view(ws.vars.theta.duals, :, part.position),
        view(ws.vars.v_mag.duals, :, part.position),
        ws.vars.theta.rho[],
        ws.vars.v_mag.rho[],
    )

    lin_p_inj = build_linear_inj_term_groups(
        part.position,
        ws.inj_groups,
        ws.inj_copy_zone,
        ws.p_inj_vars,
        ws.vars.p_inj.consensus,
        view(ws.vars.p_inj.duals, :, part.position),
        ws.vars.p_inj.rho[],
    )
    lin_q_inj = build_linear_inj_term_groups(
        part.position,
        ws.inj_groups,
        ws.inj_copy_zone,
        ws.q_inj_vars,
        ws.vars.q_inj.consensus,
        view(ws.vars.q_inj.duals, :, part.position),
        ws.vars.q_inj.rho[],
    )

    @objective(
        part.model, Min,
        part.obj_base + part.quad_bus_expr +
        part.quad_p_inj_expr + part.quad_q_inj_expr +
        lin_bus + lin_p_inj + lin_q_inj,
    )
    return nothing
end
