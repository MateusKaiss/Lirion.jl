# ============================================================
# Objective terms for ADMM (quadratic caches + linear terms)
# ============================================================

# ------------------------------------------------------------
# Quadratic caches
# ------------------------------------------------------------

# ----------------------------
# AC: bus quadratic cache
# ----------------------------
function build_quad_bus_expr!(
    part::PartitionAC,
    rho_theta::Float64,
    rho_v_mag::Float64,
)
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    v_mag = part.v_mag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    quad = zero(JuMP.QuadExpr)

    @inbounds for b in part.zone_shared_bus_ids
        quad += 0.5 * rho_theta * theta[b, 1]^2 + 0.5 * rho_v_mag * v_mag[b, 1]^2
    end

    part.quad_bus_expr_ac = quad
    part.rho_bus_theta_cached = rho_theta
    part.rho_bus_v_cached = rho_v_mag
    return quad
end

# ----------------------------
# DC: bus quadratic cache
# ----------------------------
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
    part.rho_bus_cached = rho_theta
    return quad
end

# ----------------------------
# SOCP: bus quadratic cache
# ----------------------------
function build_quad_bus_expr!(
    sub::PartitionSOCP,
    rho_v2::Float64,
)
    m = sub.model::JuMP.Model
    v2 = m[:v_magnitude_sq]
    quad = zero(JuMP.QuadExpr)

    @inbounds for b in sub.zone_shared_bus_ids
        quad += 0.5 * rho_v2 * v2[b, 1]^2
    end

    sub.quad_bus_expr_socp = quad
    sub.rho_bus_v2_cached = rho_v2
    return quad
end

# ----------------------------
# AC: branch quadratic cache (P/Q)
# ----------------------------
function build_quad_branch_expr!(
    part::PartitionAC,
    rhoP::Float64,
    rhoQ::Float64,
    zone_edges::Vector{Int},
)
    quad = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l_local = part.branch_idx_full_to_zone[e]
        l_local == 0 && continue

        P_loc, Q_loc = _pq_flow_from_exprs(part, l_local)
        quad += 0.5 * rhoP * P_loc^2 + 0.5 * rhoQ * Q_loc^2
    end

    part.quad_branch_expr_ac = quad
    part.rho_branch_p_cached = rhoP
    part.rho_branch_q_cached = rhoQ
    return quad
end

# ----------------------------
# DC: branch quadratic cache
# ----------------------------
function build_quad_branch_expr!(
    part::PartitionDC,
    rho_flow::Float64,
    zone_edges::Vector{Int},
)
    m = part.model::JuMP.Model
    quad = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l_local = part.branch_idx_full_to_zone[e]
        l_local == 0 && continue

        P = m[:p_flow_branches][l_local, 1]
        quad += 0.5 * rho_flow * P^2
    end

    part.quad_branch_expr = quad
    part.rho_branch_cached = rho_flow
    return quad
end

# ----------------------------
# SOCP: branch quadratic cache (R/T)
# ----------------------------
function build_quad_branch_rt_expr!(
    sub::PartitionSOCP,
    rho_r::Float64,
    rho_t::Float64,
    zone_edges::Vector{Int},
)
    m = sub.model::JuMP.Model
    quad = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l_local = sub.branch_idx_full_to_zone[e]
        l_local == 0 && continue

        Rloc, Tloc = _rt_from_exprs(m, l_local)
        quad += 0.5 * rho_r * Rloc^2 + 0.5 * rho_t * Tloc^2
    end

    sub.quad_rt_expr_socp = quad
    sub.rho_branch_r_cached = rho_r
    sub.rho_branch_t_cached = rho_t
    return quad
end

# ----------------------------
# SOCP: branch quadratic cache (P/Q)
# ----------------------------
function build_quad_branch_pq_expr!(
    sub::PartitionSOCP,
    rhoP::Float64,
    rhoQ::Float64,
    zone_edges::Vector{Int},
)
    m = sub.model::JuMP.Model
    quad = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l_local = sub.branch_idx_full_to_zone[e]
        l_local == 0 && continue

        P_loc, Q_loc = _pq_flow_from_exprs(m, l_local)
        quad += 0.5 * rhoP * P_loc^2 + 0.5 * rhoQ * Q_loc^2
    end

    sub.quad_pq_expr_socp = quad
    sub.rho_branch_p_cached = rhoP
    sub.rho_branch_q_cached = rhoQ
    return quad
end

# ------------------------------------------------------------
# Linear terms
# ------------------------------------------------------------

# ----------------------------
# BUS linear term (AC)
# ----------------------------
function build_linear_bus_term(
    part::PartitionAC,
    shared_buses_sorted::Vector{Int},       # compat (não usado)
    shared_buses_to_pos::Dict{Int, Int},
    theta_consensus::Vector{Float64},
    v_mag_consensus::Vector{Float64},
    theta_dual_buf_z::Vector{Float64},
    v_mag_dual_buf_z::Vector{Float64},
    rho_theta::Float64,
    rho_v::Float64,
)::JuMP.AffExpr
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    v_mag = part.v_mag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)

    @inbounds for b in part.zone_shared_bus_ids
        j = shared_buses_to_pos[b]
        lin += (theta_dual_buf_z[j] - rho_theta * theta_consensus[j]) * theta[b, 1]
        lin += (v_mag_dual_buf_z[j] - rho_v * v_mag_consensus[j]) * v_mag[b, 1]
    end
    return lin
end

# ----------------------------
# BUS linear term (DC)
# ----------------------------
function build_linear_bus_term(
    part::PartitionDC,
    shared_buses_sorted::Vector{Int},       # compat (não usado)
    shared_buses_to_pos::Dict{Int, Int},
    theta_consensus::Vector{Float64},
    theta_dual_buf_z::Vector{Float64},
    rho_theta::Float64,
)::JuMP.AffExpr
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)

    @inbounds for b in part.zone_shared_bus_ids
        j = shared_buses_to_pos[b]
        lin += (theta_dual_buf_z[j] - rho_theta * theta_consensus[j]) * theta[b, 1]
    end
    return lin
end

# ----------------------------
# BUS linear term (SOCP)
# ----------------------------
function build_linear_bus_term(
    sub::PartitionSOCP,
    shared_buses_sorted::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int},
    v2_consensus::Vector{Float64},
    v2_dual_buf_z::Vector{Float64},
    rho_v2::Float64,
)::JuMP.AffExpr
    m = sub.model::JuMP.Model
    v2 = m[:v_magnitude_sq]
    lin = zero(JuMP.AffExpr)

    @inbounds for b in sub.zone_shared_bus_ids
        j = shared_buses_to_pos[b]
        lin += (v2_dual_buf_z[j] - rho_v2 * v2_consensus[j]) * v2[b, 1]
    end
    return lin
end

# ----------------------------
# BRANCH linear term (AC)
# ----------------------------
function build_linear_branch_term(
    part::PartitionAC,
    zone_edges::Vector{Int},
    p_consensus::Dict{Int, Float64},
    q_consensus::Dict{Int, Float64},
    duals_p_flow::Dict{Int, Float64},
    duals_q_flow::Dict{Int, Float64},
    rhoP::Float64,
    rhoQ::Float64,
)::JuMP.NonlinearExpr
    lin = zero(JuMP.AffExpr)

    @inbounds for e in zone_edges
        l_local = part.branch_idx_full_to_zone[e]
        l_local == 0 && continue

        P, Q = _pq_flow_from_exprs(part, l_local)

        zP = get(p_consensus, e, 0.0)
        zQ = get(q_consensus, e, 0.0)
        dP = get(duals_p_flow, e, 0.0)
        dQ = get(duals_q_flow, e, 0.0)

        lin += dP * (P - zP) + dQ * (Q - zQ) - rhoP * zP * P - rhoQ * zQ * Q
    end
    return lin
end

# ----------------------------
# BRANCH linear term (DC)
# ----------------------------
function build_linear_branch_term(
    part::PartitionDC,
    zone_edges::Vector{Int},
    flow_consensus::Dict{Int, Float64},
    duals_flow::Dict{Int, Float64},
    rho_flow::Float64,
)::JuMP.AffExpr
    m = part.model::JuMP.Model
    lin = zero(JuMP.AffExpr)

    @inbounds for e in zone_edges
        l_local = part.branch_idx_full_to_zone[e]
        l_local == 0 && continue

        P = m[:p_flow_branches][l_local, 1]
        z = get(flow_consensus, e, 0.0)
        d = get(duals_flow, e, 0.0)

        lin += (d - rho_flow * z) * P
    end
    return lin
end

# ----------------------------
# BRANCH linear term (SOCP): RT
# ----------------------------
function build_linear_branch_rt_term(
    sub::PartitionSOCP,
    zone_edges::Vector{Int},
    r_consensus::Dict{Int, Float64},
    t_consensus::Dict{Int, Float64},
    duals_r::Dict{Int, Float64},
    duals_t::Dict{Int, Float64},
    rho_r::Float64,
    rho_t::Float64,
)::JuMP.AffExpr
    m = sub.model::JuMP.Model
    lin = zero(JuMP.AffExpr)

    @inbounds for e in zone_edges
        l_local = sub.branch_idx_full_to_zone[e]
        l_local == 0 && continue

        R, T = _rt_from_exprs(m, l_local)
        lin += (get(duals_r, e, 0.0) - rho_r * get(r_consensus, e, 0.0)) * R
        lin += (get(duals_t, e, 0.0) - rho_t * get(t_consensus, e, 0.0)) * T
    end
    return lin
end

# ----------------------------
# BRANCH linear term (SOCP): PQ
# ----------------------------
function build_linear_branch_pq_term(
    sub::PartitionSOCP,
    zone_edges::Vector{Int},
    p_consensus::Dict{Int, Float64},
    q_consensus::Dict{Int, Float64},
    duals_p_flow::Dict{Int, Float64},
    duals_q_flow::Dict{Int, Float64},
    rhoP::Float64,
    rhoQ::Float64,
)::JuMP.AffExpr
    m = sub.model::JuMP.Model
    lin = zero(JuMP.AffExpr)

    @inbounds for e in zone_edges
        l_local = sub.branch_idx_full_to_zone[e]
        l_local == 0 && continue

        P, Q = _pq_flow_from_exprs(m, l_local)
        lin += (get(duals_p_flow, e, 0.0) - rhoP * get(p_consensus, e, 0.0)) * P
        lin += (get(duals_q_flow, e, 0.0) - rhoQ * get(q_consensus, e, 0.0)) * Q
    end
    return lin
end

# ------------------------------------------------------------
# SOCP extras: soft anchoring (stabilization)
# ------------------------------------------------------------
function soft_anchor_branch_expr!(
    sub::PartitionSOCP;
    κV2::Float64 = 1e-6,
    κPQ::Float64 = 1e-8,
    use_bus::Bool = true,
    use_pq::Bool = true,
)
    m = sub.model::JuMP.Model
    quad = zero(JuMP.QuadExpr)

    if use_bus && haskey(m, :v_magnitude_sq)
        v2 = m[:v_magnitude_sq]
        sub.v2_var = v2

        for I in keys(v2)
            x = v2[I]
            xref = (JuMP.start_value(x) === nothing) ? 1.0 : JuMP.start_value(x)
            quad += 0.5 * κV2 * (x - xref)^2
        end
    end

    if use_pq && !isempty(sub.zone_edges)
        @inbounds for e in sub.zone_edges
            l = sub.branch_idx_full_to_zone[e]
            l == 0 && continue

            P, Q = _pq_flow_from_exprs(m, l)
            quad += 0.5 * κPQ * (P^2 + Q^2)
        end
    end

    sub.soft_anchor_expr = quad
    return quad
end
