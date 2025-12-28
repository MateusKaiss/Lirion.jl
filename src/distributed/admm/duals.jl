# ============================================================
# Dual updates (ADMM)
# ============================================================

# ------------------------------------------------------------
# AC duals: bus variables (theta, |V|) - vector form
# ------------------------------------------------------------
@inline function update_duals_bus_ac_all!(
    duals_theta::Dict{Int, Vector{Float64}},
    duals_v::Dict{Int, Vector{Float64}},
    partitions::Partitions{PartitionAC},
    zone_ids::AbstractVector{Int},
    theta_consensus::Vector{Float64},
    vmag_consensus::Vector{Float64},
    rho_theta::Base.RefValue{Float64},
    rho_v::Base.RefValue{Float64},
    shared_buses_sorted::AbstractVector{Int},
    shared_buses_to_pos::Dict{Int, Int};
    α = nothing,
    theta_consensus_prev = nothing,
    vmag_consensus_prev = nothing,
)
    ρθ = rho_theta[]
    ρv = rho_v[]

    @inbounds for z in zone_ids
        part = getzone(partitions, z)
        update_duals_bus_ac_vec!(
            duals_theta[z],
            duals_v[z],
            part,
            theta_consensus,
            vmag_consensus,
            ρθ,
            ρv,
            shared_buses_sorted,
            shared_buses_to_pos;
            α = α,
            theta_consensus_prev = theta_consensus_prev,
            vmag_consensus_prev = vmag_consensus_prev,
        )
    end
    return nothing
end

@inline function update_duals_bus_ac_vec!(
    duals_theta_vec_zone::Vector{Float64},
    dualsV_vec_zone::Vector{Float64},
    part::PartitionAC,
    theta_consensus_vec::Vector{Float64},
    vmag_consensus_vec::Vector{Float64},
    rho_theta::Float64,
    rhoV::Float64,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int};
    α::Union{Nothing, Float64} = nothing,
    theta_consensus_prev::Union{Nothing, Vector{Float64}} = nothing,
    vmag_consensus_prev::Union{Nothing, Vector{Float64}} = nothing,
)
    thetas = snapshot_theta_ac_shared!(part)
    Vs = snapshot_vmag_shared!(part)
    lbs = part.shared_buses_local

    if α === nothing
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            dtheta = _wrap_angle_diff(thetas[i] - theta_consensus_vec[j])
            dV = (Vs[i] - vmag_consensus_vec[j])
            duals_theta_vec_zone[j] += rho_theta * dtheta
            dualsV_vec_zone[j] += rhoV * dV
        end
    else
        αv = something(α, 1.0)
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            thetarel = _wrap_angle_diff(
                αv * thetas[i] + (1 - αv) * theta_consensus_prev[j] -
                theta_consensus_vec[j],
            )
            Vrel = αv * Vs[i] + (1 - αv) * vmag_consensus_prev[j] - vmag_consensus_vec[j]
            duals_theta_vec_zone[j] += rho_theta * thetarel
            dualsV_vec_zone[j] += rhoV * Vrel
        end
    end

    return nothing
end

# ------------------------------------------------------------
# AC duals: branch variables (P, Q) - Dict form
# ------------------------------------------------------------

function update_duals_branch_ac!(
    dualsP::Dict{Int, Dict{Int, Float64}},
    dualsQ::Dict{Int, Dict{Int, Float64}},
    partitions::Partitions{PartitionAC},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    cut_edges::Vector{Int},
    rhoP::Float64,
    rhoQ::Float64;
    α::Union{Nothing, Float64} = nothing,
    zP_prev::Union{Nothing, Dict{Int, Float64}} = nothing,
    zQ_prev::Union{Nothing, Dict{Int, Float64}} = nothing,
)
    for part in partitions
        z = part.id
        m = part.model::JuMP.Model
        lambdaPz = get!(dualsP, z, Dict{Int, Float64}())
        lambdaQz = get!(dualsQ, z, Dict{Int, Float64}())

        @inbounds for e in part.zone_edges
            l_local = get(part.edge_to_local, e, 0)
            l_local == 0 && continue

            P_loc, Q_loc = _ac_flow_from_exprs_existing(part, l_local)
            pval = JuMP.value(P_loc)
            qval = JuMP.value(Q_loc)

            if α === nothing
                if isfinite(pval)
                    lambdaPz[e] = get(lambdaPz, e, 0.0) + rhoP * (pval - get(zP, e, 0.0))
                end
                if isfinite(qval)
                    lambdaQz[e] = get(lambdaQz, e, 0.0) + rhoQ * (qval - get(zQ, e, 0.0))
                end
            else
                αv = something(α, 1.0)
                if isfinite(pval)
                    p_rel = αv * pval + (1 - αv) * get(zP_prev, e, get(zP, e, 0.0))
                    lambdaPz[e] = get(lambdaPz, e, 0.0) + rhoP * (p_rel - get(zP, e, 0.0))
                end
                if isfinite(qval)
                    q_rel = αv * qval + (1 - αv) * get(zQ_prev, e, get(zQ, e, 0.0))
                    lambdaQz[e] = get(lambdaQz, e, 0.0) + rhoQ * (q_rel - get(zQ, e, 0.0))
                end
            end
        end
    end

    return nothing
end

# ------------------------------------------------------------
# DC duals: bus variables (theta) - vector form
# ------------------------------------------------------------

@inline function update_duals_theta_all!(
    dual_theta::Dict{Int, Vector{Float64}},
    partitions::Partitions{PartitionDC},
    zone_ids::AbstractVector{Int},
    theta_consensus::Vector{Float64},
    rho::Base.RefValue{Float64},
    shared_buses_sorted::AbstractVector{Int},
    shared_buses_to_pos::Dict{Int, Int},
)
    ρ = rho[]
    @inbounds for z in zone_ids
        part = getzone(partitions, z)
        update_duals_theta_vec!(
            dual_theta[z],
            part,
            theta_consensus,
            ρ,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
    end
    return nothing
end

@inline function update_duals_theta_vec!(
    dual_vec_zone::Vector{Float64},
    part::PartitionDC,
    theta_consensus_vec::Vector{Float64},
    rho::Float64,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    thetas = snapshot_theta_shared!(part)
    lbs = part.shared_buses_local

    @inbounds for i in eachindex(lbs)
        j = b2pos[lbs[i]]
        dual_vec_zone[j] += rho * (thetas[i] - theta_consensus_vec[j])
    end
    return nothing
end

# ------------------------------------------------------------
# DC duals: branch variables (flow on cut edges) - Dict form
# ------------------------------------------------------------

function update_duals_flow!(
    duals::Dict{Int, Dict{Int, Float64}},
    partitions::Partitions{PartitionDC},
    z_flow::Dict{Int, Float64},
    cut_edges::Vector{Int},
    rho::Float64,
)
    for part in partitions
        z = part.id
        m = part.model::JuMP.Model
        dz = get!(duals, z, Dict{Int, Float64}())

        @inbounds for e in part.zone_edges
            l_local = get(part.edge_to_local, e, 0)
            l_local == 0 && continue
            p_val = _dc_flow_val_model(m, part, l_local)
            dz[e] = get(dz, e, 0.0) + rho * (p_val - get(z_flow, e, 0.0))
        end
    end
    return nothing
end

# ------------------------------------------------------------
# SOCP duals: bus variables (V^2) - vector form
# ------------------------------------------------------------

@inline function update_duals_v2_vec!(
    dual_vec_zone::Vector{Float64},
    sub::PartitionSOCP,
    v2bar_vec::Vector{Float64},
    v2bar_prev_vec::Vector{Float64},
    rho::Float64,
    α::Union{Nothing, Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    v2s = snapshot_v2_shared!(sub)
    lbs = sub.shared_buses_local

    if isnothing(α)
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            dual_vec_zone[j] += rho * (v2s[i] - v2bar_vec[j])
        end
    else
        αv = something(α, 1.0)
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            v_rel = αv * v2s[i] + (1 - αv) * v2bar_prev_vec[j]
            dual_vec_zone[j] += rho * (v_rel - v2bar_vec[j])
        end
    end

    return nothing
end

@inline function update_duals_v2_all!(
    duals_v2::Dict{Int, Vector{Float64}},
    partitions::Partitions{PartitionSOCP},
    zone_ids::AbstractVector{Int},
    v2bar::Vector{Float64},
    v2bar_prev::Vector{Float64},
    rhoV2::Float64,
    α_used,
    shared_buses_sorted::AbstractVector{Int},
    shared_buses_to_pos::Dict{Int, Int},
)
    @inbounds for z in zone_ids
        sub = getzone(partitions, z)
        update_duals_v2_vec!(
            duals_v2[z],
            sub,
            v2bar,
            v2bar_prev,
            rhoV2,
            α_used,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# SOCP duals: branch variables (R,T) - Dict form
# ------------------------------------------------------------

@inline function update_duals_rt_all!(
    dualsR::Dict{Int, Dict{Int, Float64}},
    dualsT::Dict{Int, Dict{Int, Float64}},
    partitions::Partitions{PartitionSOCP},
    zone_ids::AbstractVector{Int},
    zR::Dict{Int, Float64},
    zT::Dict{Int, Float64},
    rhoR::Float64,
    rhoT::Float64,
)
    @inbounds for z in zone_ids
        sub = getzone(partitions, z)
        m = sub.model::JuMP.Model
        dR = get!(dualsR, z, Dict{Int, Float64}())
        dT = get!(dualsT, z, Dict{Int, Float64}())

        @inbounds for e in sub.zone_edges
            l = get(sub.edge_to_local, e, 0)
            l == 0 && continue
            R, T = _rt_from_exprs_socp(m, l)
            dR[e] = get(dR, e, 0.0) + rhoR * (JuMP.value(R) - get(zR, e, 0.0))
            dT[e] = get(dT, e, 0.0) + rhoT * (JuMP.value(T) - get(zT, e, 0.0))
        end
    end
    return nothing
end

# ------------------------------------------------------------
# SOCP duals: branch variables (P,Q) - Dict form
# ------------------------------------------------------------

@inline function update_duals_pq_all!(
    dualsP::Dict{Int, Dict{Int, Float64}},
    dualsQ::Dict{Int, Dict{Int, Float64}},
    partitions::Partitions{PartitionSOCP},
    zone_ids::AbstractVector{Int},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    rhoP::Float64,
    rhoQ::Float64,
)
    @inbounds for z in zone_ids
        sub = getzone(partitions, z)
        m = sub.model::JuMP.Model
        dP = get!(dualsP, z, Dict{Int, Float64}())
        dQ = get!(dualsQ, z, Dict{Int, Float64}())

        @inbounds for e in sub.zone_edges
            l = get(sub.edge_to_local, e, 0)
            l == 0 && continue
            P, Q = _socp_flow_from_exprs_existing(m, l)
            dP[e] = get(dP, e, 0.0) + rhoP * (JuMP.value(P) - get(zP, e, 0.0))
            dQ[e] = get(dQ, e, 0.0) + rhoQ * (JuMP.value(Q) - get(zQ, e, 0.0))
        end
    end
    return nothing
end
