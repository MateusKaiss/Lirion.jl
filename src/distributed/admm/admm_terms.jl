function _build_quad_bus_expr!(part::PartitionDC, rho::Float64)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    q = zero(JuMP.QuadExpr)

    # termo (rho/2)*‖θ_shared‖² pra cada barra compartilhada
    for b in part.shared_buses_local
        q += 0.5 * rho * θ[b, 1]^2
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
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    q = zero(JuMP.QuadExpr)

    # termo (rho/2)*‖p_flow(e)‖² para cada ramo de corte da zona
    for e in zone_edges
        l_local = get(part.edge_to_local, e, nothing)
        l_local === nothing && continue
        f = part.branch.from[l_local]
        t = part.branch.to[l_local]
        x = part.branch.x[l_local]
        @assert abs(x) > 1e-10 "x≈0 no ramo global $e (DC)."
        p = (θ[f, 1] - θ[t, 1]) / x
        q += 0.5 * rho * p^2
    end
    part.quad_branch_expr = q
    part.rho_branch_cached = rho
    return q
end

#############################
# EXPRESSÕES LINEARES       #
#############################

function _build_linear_bus_expr(
    part::PartitionDC,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    θ̄_vec::Vector{Float64},
    lambda_vec_zone::Vector{Float64},
    rho::Float64,
)
    m = part.model::JuMP.Model 
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr) #cria expressão linear vazia
    for b in part.shared_buses_local
        j = b2pos[b]
        lin += (lambda_vec_zone[j] - rho * θ̄_vec[j]) * θ[b, 1]
    end
    return lin
end

function _build_linear_branch_expr(
    part::PartitionDC,
    zone_edges::Vector{Int},
    z_flow::Dict{Int, Float64},
    lambdae_zone::Dict{Int, Float64},
    rho::Float64,
)
    # soma_e ( lambda_e - rho z_e ) p_flow_e(θ)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)
    for e in zone_edges
        l_local = get(part.edge_to_local, e, nothing)
        l_local === nothing && continue
        f = part.branch.from[l_local]
        t = part.branch.to[l_local]
        x = part.branch.x[l_local]
        @assert abs(x) > 1e-10 "x≈0 no ramo global $e (DC)."
        coef = get(lambdae_zone, e, 0.0) - rho * z_flow[e]
        lin += coef * (θ[f, 1] - θ[t, 1]) / x
    end
    return lin
end

function _build_quad_bus_expr_ac!(part::PartitionAC, rhoθ::Float64, rhoV::Float64)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    V = part.vmag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    q = zero(JuMP.QuadExpr)
    for b in part.shared_buses_local
        q += 0.5 * rhoθ * θ[b, 1]^2 + 0.5 * rhoV * V[b, 1]^2
    end
    part.quad_bus_expr_ac = q
    part.rho_bus_theta_cached = rhoθ
    part.rho_bus_v_cached = rhoV
    return q
end

function _build_linear_bus_expr_ac!(
    part::PartitionAC,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    lambdaθ_zone::Vector{Float64},
    lambdaV_zone::Vector{Float64},
    rhoθ::Float64,
    rhoV::Float64,
)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    V = part.vmag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)
    for b in part.shared_buses_local
        j = b2pos[b]
        lin += (lambdaθ_zone[j] - rhoθ * θ̄_vec[j]) * θ[b, 1]
        lin += (lambdaV_zone[j] - rhoV * V̄_vec[j]) * V[b, 1]
    end
    return lin
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
        lin += lambdaP * (P_loc - zPe) + lambdaQ * (Q_loc - zQe) -
               rhoP * zPe * P_loc - rhoQ * zQe * Q_loc
    end
    return lin
end


"Ridge suave para estabilizar QCP no runner branch (P,Q)."
function _soft_anchor_branch_expr!(sub::PartitionSOCP;
    κV2::Float64 = 1e-6, κPQ::Float64 = 1e-8,
    use_bus::Bool = true, use_pq::Bool = true,
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
            P, Q = _socp_flow_from_exprs_existing(m, l)   # deve ser AffExpr!
            q += 0.5 * κPQ * (P^2 + Q^2)
        end
    end

    sub.soft_anchor_expr = q
    return q
end

@inline function _socp_flow_from_exprs_existing(
    model::Model, l_local::Int; t::Int = 1,
)
    return model[:p_flow_branches][l_local, t], model[:q_flow_branches][l_local, t]
end

@inline function _rt_from_exprs_socp(model::Model, l_local::Int; t::Int = 1)
    return model[:R][l_local, t], model[:T_var][l_local, t]
end

# ----------------- Construção de termos (cache de quadrático) -----------------
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

function _build_quad_branch_rt_expr!(sub::PartitionSOCP, rhoR::Float64, rhoT::Float64,
    zone_edges::Vector{Int})
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

function _build_quad_branch_pq_expr!(sub::PartitionSOCP, rhoP::Float64, rhoQ::Float64,
    zone_edges::Vector{Int})
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

function _build_linear_bus_expr_socp(sub::PartitionSOCP,
    shared_buses_sorted::Vector{Int}, b2pos::Dict{Int, Int},
    v2bar_vec::Vector{Float64}, lambda_vec_zone::Vector{Float64}, rhoV2::Float64)
    m = sub.model::JuMP.Model
    V2 = m[:v_magnitude_sq]
    lin = zero(JuMP.AffExpr)
    for b in sub.shared_buses_local
        j = b2pos[b]
        lin += (lambda_vec_zone[j] - rhoV2 * v2bar_vec[j]) * V2[b, 1]
    end
    return lin
end

function _build_linear_branch_rt_expr(sub::PartitionSOCP,
    zone_edges::Vector{Int}, zR::Dict{Int, Float64}, zT::Dict{Int, Float64},
    lambdaR_zone::Dict{Int, Float64}, lambdaT_zone::Dict{Int, Float64}, rhoR::Float64, rhoT::Float64)
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

function _build_linear_branch_pq_expr(sub::PartitionSOCP,
    zone_edges::Vector{Int}, zP::Dict{Int, Float64}, zQ::Dict{Int, Float64},
    lambdaP_zone::Dict{Int, Float64}, lambdaQ_zone::Dict{Int, Float64}, rhoP::Float64, rhoQ::Float64)
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