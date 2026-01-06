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

    @inbounds for b in part.shared_buses_zone
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

    @inbounds for b in part.shared_buses_zone
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
    V2 = m[:v_magnitude_sq]
    quad = zero(JuMP.QuadExpr)

    @inbounds for b in sub.shared_buses_zone
        quad += 0.5 * rho_v2 * V2[b, 1]^2
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
        l_local = get(part.branch_global_to_local, e, nothing)
        l_local === nothing && continue

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
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    quad = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l_local = get(part.branch_global_to_local, e, nothing)
        l_local === nothing && continue

        f = part.branch.from[l_local]
        t = part.branch.to[l_local]
        x = part.branch.x[l_local]
        @assert abs(x) > 1e-10 "x≈0 no ramo global $e (DC)."

        p = (theta[f, 1] - theta[t, 1]) / x
        quad += 0.5 * rho_flow * p^2
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
        l_local = get(sub.branch_global_to_local, e, nothing)
        l_local === nothing && continue

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
        l_local = get(sub.branch_global_to_local, e, nothing)
        l_local === nothing && continue

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
    shared_buses_sorted::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int},
    theta_consensus::Vector{Float64},
    v_mag_consensus::Vector{Float64},
    lambda_theta_buf::Vector{Float64},
    lambda_v_buf::Vector{Float64},
    rho_theta::Float64,
    rho_v::Float64,
)::JuMP.AffExpr
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    v_mag = part.v_mag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)

    @inbounds for b in part.shared_buses_zone
        j = shared_buses_to_pos[b]
        lin += (lambda_theta_buf[j] - rho_theta * theta_consensus[j]) * theta[b, 1]
        lin += (lambda_v_buf[j] - rho_v * v_mag_consensus[j]) * v_mag[b, 1]
    end
    return lin
end

# ----------------------------
# BUS linear term (DC)
# ----------------------------

function build_linear_bus_term(
    part::PartitionDC,
    shared_buses_sorted::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int},
    theta_consensus::Vector{Float64},
    theta_dual_buf::Vector{Float64},
    rho_theta::Float64,
)::JuMP.AffExpr
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)

    @inbounds for b in part.shared_buses_zone
        j = shared_buses_to_pos[b]
        lin += (theta_dual_buf[j] - rho_theta * theta_consensus[j]) * theta[b, 1]
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
    lambda_v2_buf::Vector{Float64},
    rho_v2::Float64,
)::JuMP.AffExpr
    m = sub.model::JuMP.Model
    V2 = m[:v_magnitude_sq]
    lin = zero(JuMP.AffExpr)

    @inbounds for b in sub.shared_buses_zone
        j = shared_buses_to_pos[b]
        lin += (lambda_v2_buf[j] - rho_v2 * v2_consensus[j]) * V2[b, 1]
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
    lambdaP_zone::Dict{Int, Float64},
    lambdaQ_zone::Dict{Int, Float64},
    rhoP::Float64,
    rhoQ::Float64,
)::JuMP.NonlinearExpr
    lin = zero(JuMP.AffExpr)

    @inbounds for e in zone_edges
        l_local = get(part.branch_global_to_local, e, nothing)
        l_local === nothing && continue

        P, Q = _pq_flow_from_exprs(part, l_local)

        zP = get(p_consensus, e, 0.0)
        zQ = get(q_consensus, e, 0.0)
        λP = get(lambdaP_zone, e, 0.0)
        λQ = get(lambdaQ_zone, e, 0.0)

        lin += λP * (P - zP) + λQ * (Q - zQ) - rhoP * zP * P - rhoQ * zQ * Q
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
    lambda_zone::Dict{Int, Float64},
    rho_flow::Float64,
)::JuMP.AffExpr
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)

    @inbounds for e in zone_edges
        l_local = get(part.branch_global_to_local, e, nothing)
        l_local === nothing && continue

        f = part.branch.from[l_local]
        t = part.branch.to[l_local]
        x = part.branch.x[l_local]
        @assert abs(x) > 1e-10 "x≈0 no ramo global $e (DC)."

        z = get(flow_consensus, e, 0.0)
        λ = get(lambda_zone, e, 0.0)
        coef = (λ - rho_flow * z)

        lin += coef * (theta[f, 1] - theta[t, 1]) / x
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
    lambdaR_zone::Dict{Int, Float64},
    lambdaT_zone::Dict{Int, Float64},
    rho_r::Float64,
    rho_t::Float64,
)::JuMP.AffExpr
    m = sub.model::JuMP.Model
    lin = zero(JuMP.AffExpr)

    @inbounds for e in zone_edges
        l_local = get(sub.branch_global_to_local, e, nothing)
        l_local === nothing && continue

        R, T = _rt_from_exprs(m, l_local)
        lin += (get(lambdaR_zone, e, 0.0) - rho_r * get(r_consensus, e, 0.0)) * R
        lin += (get(lambdaT_zone, e, 0.0) - rho_t * get(t_consensus, e, 0.0)) * T
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
    lambdaP_zone::Dict{Int, Float64},
    lambdaQ_zone::Dict{Int, Float64},
    rhoP::Float64,
    rhoQ::Float64,
)::JuMP.AffExpr
    m = sub.model::JuMP.Model
    lin = zero(JuMP.AffExpr)

    @inbounds for e in zone_edges
        l_local = get(sub.branch_global_to_local, e, nothing)
        l_local === nothing && continue

        P, Q = _pq_flow_from_exprs(m, l_local)
        lin += (get(lambdaP_zone, e, 0.0) - rhoP * get(p_consensus, e, 0.0)) * P
        lin += (get(lambdaQ_zone, e, 0.0) - rhoQ * get(q_consensus, e, 0.0)) * Q
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
        V2 = m[:v_magnitude_sq]
        sub.v2_var = V2

        for I in keys(V2)
            x = V2[I]
            xref = (JuMP.start_value(x) === nothing) ? 1.0 : JuMP.start_value(x)
            quad += 0.5 * κV2 * (x - xref)^2
        end
    end

    if use_pq && !isempty(sub.zone_edges)
        @inbounds for e in sub.zone_edges
            l = get(sub.branch_global_to_local, e, nothing)
            l === nothing && continue
            P, Q = _pq_flow_from_exprs(m, l)
            quad += 0.5 * κPQ * (P^2 + Q^2)
        end
    end

    sub.soft_anchor_expr = quad
    return quad
end
