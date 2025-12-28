# ============================================================
# Objective terms for ADMM (quadratic caches + linear terms)
# ============================================================

# ------------------------------------------------------------
# DC terms: quadratic caches (bus / branch)
# ------------------------------------------------------------

function _build_quad_bus_expr!(part::PartitionDC, rho::Float64)
    m = part.model::JuMP.Model
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    q = zero(JuMP.QuadExpr)

    for b in part.shared_buses_local
        q += 0.5 * rho * theta[b, 1]^2
    end

    part.quad_bus_expr = q
    part.rho_bus_cached = rho
    return q
end

function _build_quad_branch_expr!(
    part::PartitionDC,
    rho::Float64,
    zone_edges::Vector{Int},
)
    m = part.model::JuMP.Model
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    q = zero(JuMP.QuadExpr)

    # termo (rho/2)*‖p_flow(e)‖² para cada ramo de corte da zona
    for e in zone_edges
        l_local = get(part.edge_to_local, e, nothing)
        l_local === nothing && continue
        f = part.branch.from[l_local]
        t = part.branch.to[l_local]
        x = part.branch.x[l_local]
        @assert abs(x) > 1e-10 "x≈0 no ramo global $e (DC)."
        p = (theta[f, 1] - theta[t, 1]) / x
        q += 0.5 * rho * p^2
    end

    part.quad_branch_expr = q
    part.rho_branch_cached = rho
    return q
end

# ------------------------------------------------------------
# DC terms: linear terms (bus / branch)
# ------------------------------------------------------------

function _build_linear_bus_expr(
    part::PartitionDC,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    theta_consensus_vec::Vector{Float64},
    lambda_vec_zone::Vector{Float64},
    rho::Float64,
)
    m = part.model::JuMP.Model
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)

    for b in part.shared_buses_local
        j = b2pos[b]
        lin += (lambda_vec_zone[j] - rho * theta_consensus_vec[j]) * theta[b, 1]
    end

    return lin
end

function _build_linear_branch_expr(
    part::PartitionDC,
    zone_edges::Vector{Int},
    flow_consensus::Dict{Int, Float64},
    lambdae_zone::Dict{Int, Float64},
    rho::Float64,
)
    # soma_e ( lambda_e - rho z_e ) p_flow_e(theta)
    m = part.model::JuMP.Model
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)

    for e in zone_edges
        l_local = get(part.edge_to_local, e, nothing)
        l_local === nothing && continue
        f = part.branch.from[l_local]
        t = part.branch.to[l_local]
        x = part.branch.x[l_local]
        @assert abs(x) > 1e-10 "x≈0 no ramo global $e (DC)."
        coef = get(lambdae_zone, e, 0.0) - rho * flow_consensus[e]
        lin += coef * (theta[f, 1] - theta[t, 1]) / x
    end

    return lin
end

# ------------------------------------------------------------
# AC terms: quadratic caches (bus / branch)
# ------------------------------------------------------------

function _build_quad_bus_expr_ac!(part::PartitionAC, rho_theta::Float64, rhoV::Float64)
    m = part.model::JuMP.Model
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    V = part.vmag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    q = zero(JuMP.QuadExpr)

    for b in part.shared_buses_local
        q += 0.5 * rho_theta * theta[b, 1]^2 + 0.5 * rhoV * V[b, 1]^2
    end

    part.quad_bus_expr_ac = q
    part.rho_bus_theta_cached = rho_theta
    part.rho_bus_v_cached = rhoV
    return q
end

function _build_quad_branch_expr_ac!(
    part::PartitionAC,
    rhoP::Float64,
    rhoQ::Float64,
    zone_edges::Vector{Int},
)
    m = part.model::JuMP.Model
    q = zero(JuMP.QuadExpr)

    for e in zone_edges
        l_local = get(part.edge_to_local, e, nothing)
        l_local === nothing && continue
        P_loc, Q_loc = _ac_flow_from_exprs_existing(part, l_local)
        q += 0.5 * rhoP * P_loc^2 + 0.5 * rhoQ * Q_loc^2
    end

    part.quad_branch_expr_ac = q
    part.rho_branch_p_cached = rhoP
    part.rho_branch_q_cached = rhoQ
    return q
end

# ------------------------------------------------------------
# AC terms: linear terms (bus / branch)
# ------------------------------------------------------------

function _build_linear_bus_expr_ac!(
    part::PartitionAC,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    theta_consensus_vec::Vector{Float64},
    vmag_consensus_vec::Vector{Float64},
    lambdatheta_zone::Vector{Float64},
    lambdaV_zone::Vector{Float64},
    rho_theta::Float64,
    rhoV::Float64,
)
    m = part.model::JuMP.Model
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    V = part.vmag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)

    for b in part.shared_buses_local
        j = b2pos[b]
        lin += (lambdatheta_zone[j] - rho_theta * theta_consensus_vec[j]) * theta[b, 1]
        lin += (lambdaV_zone[j] - rhoV * vmag_consensus_vec[j]) * V[b, 1]
    end

    return lin
end

function _build_linear_branch_expr_ac!(
    part::PartitionAC,
    zone_edges::Vector{Int},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    dualsP_zone::Dict{Int, Float64},
    dualsQ_zone::Dict{Int, Float64},
    rhoP::Float64,
    rhoQ::Float64,
)
    m = part.model::JuMP.Model
    lin = zero(JuMP.AffExpr)

    for e in zone_edges
        l_local = get(part.edge_to_local, e, nothing)
        l_local === nothing && continue
        P_loc, Q_loc = _ac_flow_from_exprs_existing(part, l_local)
        lambdaP = get(dualsP_zone, e, 0.0)
        lambdaQ = get(dualsQ_zone, e, 0.0)
        zPe = get(zP, e, 0.0)
        zQe = get(zQ, e, 0.0)
        lin +=
            lambdaP * (P_loc - zPe) + lambdaQ * (Q_loc - zQe) -
            rhoP * zPe * P_loc - rhoQ * zQe * Q_loc
    end

    return lin
end

# ------------------------------------------------------------
# SOCP terms: quadratic caches (bus / branch RT / branch PQ)
# ------------------------------------------------------------

function _build_quad_bus_expr_socp!(sub::PartitionSOCP, rho::Float64)
    m = sub.model::JuMP.Model
    V2 = m[:v_magnitude_sq]
    q = zero(JuMP.QuadExpr)

    # penaliza APENAS barras compartilhadas (padrão ADMM/DC/AC)
    for b in sub.shared_buses_local
        q += 0.5 * rho * V2[b, 1]^2
    end

    sub.quad_bus_expr_socp = q
    sub.rho_bus_v2_cached = rho
    return sub.quad_bus_expr_socp
end

function _build_quad_branch_rt_expr!(
    sub::PartitionSOCP,
    rhoR::Float64,
    rhoT::Float64,
    zone_edges::Vector{Int},
)
    m = sub.model::JuMP.Model
    q = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l_local = get(sub.edge_to_local, e, nothing)
        l_local === nothing && continue
        Rloc, Tloc = _rt_from_exprs_socp(m, l_local)
        q += 0.5 * rhoR * Rloc^2 + 0.5 * rhoT * Tloc^2
    end

    sub.quad_rt_expr_socp = q
    sub.rho_branch_r_cached = rhoR
    sub.rho_branch_t_cached = rhoT
    return q
end

function _build_quad_branch_pq_expr!(
    sub::PartitionSOCP,
    rhoP::Float64,
    rhoQ::Float64,
    zone_edges::Vector{Int},
)
    m = sub.model::JuMP.Model
    q = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l_local = get(sub.edge_to_local, e, nothing)
        l_local === nothing && continue
        P_loc, Q_loc = _socp_flow_from_exprs_existing(m, l_local)
        q += 0.5 * rhoP * P_loc^2 + 0.5 * rhoQ * Q_loc^2
    end

    sub.quad_pq_expr_socp = q
    sub.rho_branch_p_cached = rhoP
    sub.rho_branch_q_cached = rhoQ
    return q
end

# ------------------------------------------------------------
# SOCP terms: linear terms (bus / branch RT / branch PQ)
# ------------------------------------------------------------

function _build_linear_bus_expr_socp(
    sub::PartitionSOCP,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    v2bar_vec::Vector{Float64},
    lambda_vec_zone::Vector{Float64},
    rhoV2::Float64,
)
    m = sub.model::JuMP.Model
    V2 = m[:v_magnitude_sq]
    lin = zero(JuMP.AffExpr)

    for b in sub.shared_buses_local
        j = b2pos[b]
        lin += (lambda_vec_zone[j] - rhoV2 * v2bar_vec[j]) * V2[b, 1]
    end

    return lin
end

function _build_linear_branch_rt_expr(
    sub::PartitionSOCP,
    zone_edges::Vector{Int},
    zR::Dict{Int, Float64},
    zT::Dict{Int, Float64},
    lambdaR_zone::Dict{Int, Float64},
    lambdaT_zone::Dict{Int, Float64},
    rhoR::Float64,
    rhoT::Float64,
)
    m = sub.model::JuMP.Model
    lin = zero(JuMP.AffExpr)

    @inbounds for e in zone_edges
        l_local = get(sub.edge_to_local, e, nothing)
        l_local === nothing && continue
        Rloc, Tloc = _rt_from_exprs_socp(m, l_local)
        lin += (get(lambdaR_zone, e, 0.0) - rhoR * get(zR, e, 0.0)) * Rloc
        lin += (get(lambdaT_zone, e, 0.0) - rhoT * get(zT, e, 0.0)) * Tloc
    end

    return lin
end

function _build_linear_branch_pq_expr(
    sub::PartitionSOCP,
    zone_edges::Vector{Int},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    lambdaP_zone::Dict{Int, Float64},
    lambdaQ_zone::Dict{Int, Float64},
    rhoP::Float64,
    rhoQ::Float64,
)
    m = sub.model::JuMP.Model
    lin = zero(JuMP.AffExpr)

    @inbounds for e in zone_edges
        l_local = get(sub.edge_to_local, e, nothing)
        l_local === nothing && continue
        P_loc, Q_loc = _socp_flow_from_exprs_existing(m, l_local)
        lin += (get(lambdaP_zone, e, 0.0) - rhoP * get(zP, e, 0.0)) * P_loc
        lin += (get(lambdaQ_zone, e, 0.0) - rhoQ * get(zQ, e, 0.0)) * Q_loc
    end

    return lin
end

# ------------------------------------------------------------
# SOCP extras: soft anchoring (stabilization)
# ------------------------------------------------------------
function _soft_anchor_branch_expr!(
    sub::PartitionSOCP;
    κV2::Float64 = 1e-6,
    κPQ::Float64 = 1e-8,
    use_bus::Bool = true,
    use_pq::Bool = true,
)
    m = sub.model::JuMP.Model
    q = zero(JuMP.QuadExpr)

    if use_bus && haskey(m, :v_magnitude_sq)
        V2 = m[:v_magnitude_sq]
        sub.v2_var = V2
        # ancora TODAS as barras da zona (não só compartilhadas)
        for I in keys(V2)
            x = V2[I]
            # start_value se houver, senão 1.0 (flat seguro p/ SOCP)
            xref = (JuMP.start_value(x) === nothing) ? 1.0 : JuMP.start_value(x)
            q += 0.5 * κV2 * (x - xref)^2
        end
    end

    if use_pq && !isempty(sub.zone_edges)
        for e in sub.zone_edges
            l = get(sub.edge_to_local, e, nothing)
            l === nothing && continue
            P, Q = _socp_flow_from_exprs_existing(m, l) # deve ser AffExpr!
            q += 0.5 * κPQ * (P^2 + Q^2)
        end
    end

    sub.soft_anchor_expr = q
    return q
end
