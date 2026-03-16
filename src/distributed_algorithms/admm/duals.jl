# ============================================================
# Dual updates (ADMM)
# ============================================================

@inline function update_duals_bus_ac_vec!(
    duals_theta::Vector{Float64},
    duals_v_mag::Vector{Float64},
    part::PartitionAC,
    theta_consensus::Vector{Float64},
    v_mag_consensus_vec::Vector{Float64},
    rho_theta::Float64,
    rho_v_mag::Float64,
    shared_buses_sorted::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int};
    theta_consensus_prev::Union{Nothing, Vector{Float64}} = nothing,
    v_mag_consensus_prev::Union{Nothing, Vector{Float64}} = nothing,
)
    theta_shared_local = snapshot_theta_shared!(part)
    v_mag_shared_local = snapshot_v_mag_shared!(part)
    shared_bus_ids = part.shared_buses_zone

    @inbounds for i in eachindex(shared_bus_ids)
        j = shared_buses_to_pos[shared_bus_ids[i]]
        dtheta = _wrap_angle_diff(theta_shared_local[i] - theta_consensus[j])
        dV = (v_mag_shared_local[i] - v_mag_consensus_vec[j])
        duals_theta[j] += rho_theta * dtheta
        duals_v_mag[j] += rho_v_mag * dV
    end

    return nothing
end

# ------------------------------------------------------------
# AC duals: bus variables (theta, |V|)
# ------------------------------------------------------------

@inline function update_duals_bus_ac_all!(
    duals_theta::Dict{Int, Vector{Float64}},
    duals_v_mag::Dict{Int, Vector{Float64}},
    partitions::Partitions{PartitionAC},
    zone_ids::AbstractVector{Int},
    theta_consensus::Vector{Float64},
    v_mag_consensus::Vector{Float64},
    rho_theta::Base.RefValue{Float64},
    rho_v::Base.RefValue{Float64},
    shared_buses_sorted::AbstractVector{Int},
    shared_buses_to_pos::Dict{Int, Int};
    theta_consensus_prev = nothing,
    v_mag_consensus_prev = nothing,
)
    n = length(zone_ids)
    @inbounds for i in 1:n
        zone = zone_ids[i]
        part = partitions.parts[i]

        update_duals_bus_ac_vec!(
            duals_theta[zone],
            duals_v_mag[zone],
            part,
            theta_consensus,
            v_mag_consensus,
            rho_theta[],
            rho_v[],
            shared_buses_sorted,
            shared_buses_to_pos;
            theta_consensus_prev = theta_consensus_prev,
            v_mag_consensus_prev = v_mag_consensus_prev,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# DC duals: bus variables (theta)
# ------------------------------------------------------------

@inline function update_duals_theta_vec!(
    dual_theta::Vector{Float64},
    part::PartitionDC,
    theta_consensus::Vector{Float64},
    rho::Float64,
    shared_buses_sorted::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int},
)
    theta_shared_local = snapshot_theta_shared!(part)
    shared_bus_ids = part.shared_buses_zone

    @inbounds for i in eachindex(shared_bus_ids)
        j = shared_buses_to_pos[shared_bus_ids[i]]
        dual_theta[j] += rho * (theta_shared_local[i] - theta_consensus[j])
    end
    return nothing
end

@inline function update_duals_theta_all!(
    dual_theta::Dict{Int, Vector{Float64}},
    partitions::Partitions{PartitionDC},
    zone_ids::AbstractVector{Int},
    theta_consensus::Vector{Float64},
    rho::Base.RefValue{Float64},
    shared_buses_sorted::AbstractVector{Int},
    shared_buses_to_pos::Dict{Int, Int},
)
    n = length(zone_ids)
    @inbounds for i in 1:n
        zone = zone_ids[i]
        part = partitions.parts[i]

        update_duals_theta_vec!(
            dual_theta[zone],
            part,
            theta_consensus,
            rho[],
            shared_buses_sorted,
            shared_buses_to_pos,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# SOCP duals: bus variables (V^2)
# ------------------------------------------------------------

@inline function update_duals_v2_vec!(
    duals_v2::Vector{Float64},
    part::PartitionSOCP,
    v2_consensus::Vector{Float64},
    rho_v2::Float64,
    shared_buses_sorted::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int},
)
    v2s = snapshot_v2_shared!(part)
    shared_bus_ids = part.shared_buses_zone

    @inbounds for i in eachindex(shared_bus_ids)
        j = shared_buses_to_pos[shared_bus_ids[i]]
        duals_v2[j] += rho_v2 * (v2s[i] - v2_consensus[j])
    end

    return nothing
end

@inline function update_duals_v2_all!(
    duals_v2::Dict{Int, Vector{Float64}},
    partitions::Partitions{PartitionSOCP},
    zone_ids::AbstractVector{Int},
    v2_consensus::Vector{Float64},
    rho_v2::Float64,
    shared_buses_sorted::AbstractVector{Int},
    shared_buses_to_pos::Dict{Int, Int},
)
    n = length(zone_ids)
    @inbounds for i in 1:n
        z = zone_ids[i]
        part = partitions.parts[i]

        update_duals_v2_vec!(
            duals_v2[z],
            part,
            v2_consensus,
            rho_v2,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# AC duals: branch variables (P, Q)
# ------------------------------------------------------------

function update_duals_branch_ac!(
    dualsP::Dict{Int, Dict{Int, Float64}},
    dualsQ::Dict{Int, Dict{Int, Float64}},
    partitions::Partitions{PartitionAC},
    p_flow_consensus::Dict{Int, Float64},
    q_flow_consensus::Dict{Int, Float64},
    interzonal_branches::Vector{Int},
    rhoP::Float64,
    rhoQ::Float64;
    alpha::Union{Nothing, Float64} = nothing,
    p_flow_consensus_prev::Union{Nothing, Dict{Int, Float64}} = nothing,
    q_flow_consensus_prev::Union{Nothing, Dict{Int, Float64}} = nothing,
)
    for part in partitions
        z = part.id
        m = part.model::JuMP.Model
        lambdaPz = get!(dualsP, z, Dict{Int, Float64}())
        lambdaQz = get!(dualsQ, z, Dict{Int, Float64}())

        @inbounds for e in part.zone_edges
            l_local = part.branch_idx_full_to_zone[e]
            l_local == 0 && continue

            P_loc, Q_loc = _pq_flow_from_exprs(part, l_local)
            pval = JuMP.value(P_loc)
            qval = JuMP.value(Q_loc)

            if alpha === nothing
                if isfinite(pval)
                    lambdaPz[e] =
                        get(lambdaPz, e, 0.0) +
                        rhoP * (pval - get(p_flow_consensus, e, 0.0))
                end
                if isfinite(qval)
                    lambdaQz[e] =
                        get(lambdaQz, e, 0.0) +
                        rhoQ * (qval - get(q_flow_consensus, e, 0.0))
                end
            else
                alphav = something(alpha, 1.0)
                if isfinite(pval)
                    p_rel =
                        alphav * pval +
                        (1 - alphav) *
                        get(p_flow_consensus_prev, e, get(p_flow_consensus, e, 0.0))
                    lambdaPz[e] =
                        get(lambdaPz, e, 0.0) +
                        rhoP * (p_rel - get(p_flow_consensus, e, 0.0))
                end
                if isfinite(qval)
                    q_rel =
                        alphav * qval +
                        (1 - alphav) *
                        get(q_flow_consensus_prev, e, get(q_flow_consensus, e, 0.0))
                    lambdaQz[e] =
                        get(lambdaQz, e, 0.0) +
                        rhoQ * (q_rel - get(q_flow_consensus, e, 0.0))
                end
            end
        end
    end

    return nothing
end

# ------------------------------------------------------------
# DC duals: branch variables (flow on cut edges)
# ------------------------------------------------------------

function update_duals_flow!(
    duals::Dict{Int, Dict{Int, Float64}},
    partitions::Partitions{PartitionDC},
    flow_consensus::Dict{Int, Float64},
    interzonal_branches::Vector{Int},
    rho::Float64,
)
    for part in partitions
        z = part.id
        m = part.model::JuMP.Model
        dz = get!(duals, z, Dict{Int, Float64}())

        @inbounds for e in part.zone_edges
            l_local = part.branch_idx_full_to_zone[e]
            l_local == 0 && continue

            p_val = _dc_flow_val_model(m, part, l_local)
            dz[e] = get(dz, e, 0.0) + rho * (p_val - get(flow_consensus, e, 0.0))
        end
    end
    return nothing
end

# ------------------------------------------------------------
# SOCP duals: branch variables (R,T)
# ------------------------------------------------------------

@inline function update_duals_rt_all!(
    duals_r::Dict{Int, Dict{Int, Float64}},
    duals_t::Dict{Int, Dict{Int, Float64}},
    partitions::Partitions{PartitionSOCP},
    zone_ids::AbstractVector{Int},
    r_consensus::Dict{Int, Float64},
    t_consensus::Dict{Int, Float64},
    rho_r::Float64,
    rho_t::Float64,
)
    n = length(zone_ids)
    @inbounds for i in 1:n
        z = zone_ids[i]
        part = partitions.parts[i]

        m = part.model::JuMP.Model
        dR = get!(duals_r, z, Dict{Int, Float64}())
        dT = get!(duals_t, z, Dict{Int, Float64}())

        @inbounds for e in part.zone_edges
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue

            R, T = _rt_from_exprs(m, l)
            dR[e] = get(dR, e, 0.0) + rho_r * (JuMP.value(R) - get(r_consensus, e, 0.0))
            dT[e] = get(dT, e, 0.0) + rho_t * (JuMP.value(T) - get(t_consensus, e, 0.0))
        end
    end
    return nothing
end

# ------------------------------------------------------------
# SOCP duals: branch variables (P,Q)
# ------------------------------------------------------------

@inline function update_duals_pq_all!(
    dualsP::Dict{Int, Dict{Int, Float64}},
    dualsQ::Dict{Int, Dict{Int, Float64}},
    partitions::Partitions{PartitionSOCP},
    zone_ids::AbstractVector{Int},
    p_flow_consensus::Dict{Int, Float64},
    q_flow_consensus::Dict{Int, Float64},
    rhoP::Float64,
    rhoQ::Float64,
)
    n = length(zone_ids)
    @inbounds for i in 1:n
        z = zone_ids[i]
        part = partitions.parts[i]

        m = part.model::JuMP.Model
        dP = get!(dualsP, z, Dict{Int, Float64}())
        dQ = get!(dualsQ, z, Dict{Int, Float64}())

        @inbounds for e in part.zone_edges
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue

            P, Q = _pq_flow_from_exprs(m, l)
            dP[e] = get(dP, e, 0.0) + rhoP * (JuMP.value(P) - get(p_flow_consensus, e, 0.0))
            dQ[e] = get(dQ, e, 0.0) + rhoQ * (JuMP.value(Q) - get(q_flow_consensus, e, 0.0))
        end
    end
    return nothing
end
