# ============================================================
# Residuals and tolerances (ADMM)
# ============================================================

# ------------------------------------------------------------
# AC residuals: primal (bus / branch)
# ------------------------------------------------------------

function primal_residual_bus_ac_vec(
    partitions::Partitions{PartitionAC},
    theta_consensus_vec::Vector{Float64},
    v_mag_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int},
)
    acc = 0.0
    for part in partitions
        thetabuf = snapshot_theta_ac!(part)
        Vbuf = snapshot_v_mag!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_zone
            j = shared_buses_to_pos[b]
            acc +=
                _wrap_angle_diff(thetabuf[pos[b]] - theta_consensus_vec[j])^2 +
                (Vbuf[pos[b]] - v_mag_consensus_vec[j])^2
        end
    end
    return sqrt(acc)
end

function primal_residual_bus_ac_components(
    partitions::Partitions{PartitionAC},
    theta_consensus_vec::Vector{Float64},
    v_mag_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int},
)
    rtheta = 0.0
    rV = 0.0
    for part in partitions
        thetabuf = snapshot_theta_ac!(part)
        Vbuf = snapshot_v_mag!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_zone
            j = shared_buses_to_pos[b]
            rtheta += _wrap_angle_diff(thetabuf[pos[b]] - theta_consensus_vec[j])^2
            rV += (Vbuf[pos[b]] - v_mag_consensus_vec[j])^2
        end
    end
    return sqrt(rtheta), sqrt(rV)
end

function primal_residual_branch_ac(
    partitions::Partitions{PartitionAC},
    p_flow_consensus::Dict{Int, Float64},
    q_flow_consensus::Dict{Int, Float64},
    interzonal_branches::Vector{Int},
)
    acc = 0.0
    for part in partitions
        @inbounds for e in interzonal_branches
            l_local = get(part.branch_global_to_local, e, 0)
            l_local == 0 && continue

            P_loc, Q_loc = _pq_flow_from_exprs(part, l_local; t = 1)
            p = JuMP.value(P_loc)
            q = JuMP.value(Q_loc)

            if isfinite(p)
                acc += (p - get(p_flow_consensus, e, 0.0))^2
            end
            if isfinite(q)
                acc += (q - get(q_flow_consensus, e, 0.0))^2
            end
        end
    end
    return sqrt(acc)
end

# ------------------------------------------------------------
# DC residuals: primal (bus / branch)
# ------------------------------------------------------------

function primal_residual_bus_dc(
    partitions::Partitions{PartitionDC},
    theta_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int},
)
    acc = 0.0
    for part in partitions
        thetabuf = snapshot_theta!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_zone
            j = shared_buses_to_pos[b]
            acc += (thetabuf[pos[b]] - theta_consensus_vec[j])^2
        end
    end
    return sqrt(acc)
end

function primal_residual_branch_dc(
    partitions::Partitions{PartitionDC},
    flow_consensus::Dict{Int, Float64},
    interzonal_branches::Vector{Int},
)
    acc = 0.0
    for part in partitions
        m = part.model::JuMP.Model
        @inbounds for e in interzonal_branches
            l_local = get(part.branch_global_to_local, e, 0)
            l_local == 0 && continue
            acc += (_dc_flow_val_model(m, part, l_local) - get(flow_consensus, e, 0.0))^2
        end
    end
    return sqrt(acc)
end

# ------------------------------------------------------------
# SOCP residuals: primal (bus / branch RT / branch PQ)
# ------------------------------------------------------------

function primal_residual_bus_socp_vec(
    partitions::Partitions{PartitionSOCP},
    v2_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int}, # compat (não usado)
    shared_buses_to_pos::Dict{Int, Int},
)
    acc = 0.0
    for sub in partitions
        v2buf = snapshot_v2!(sub)
        pos = sub.bus_local_index
        @inbounds for b in sub.shared_buses_zone
            j = shared_buses_to_pos[b]
            acc += (v2buf[pos[b]] - v2_consensus_vec[j])^2
        end
    end
    return sqrt(acc)
end

function primal_residual_branch_rt_socp(
    partitions::Partitions{PartitionSOCP},
    r_consensus::Dict{Int, Float64},
    t_consensus::Dict{Int, Float64},
    interzonal_branches::Vector{Int},
)
    acc = 0.0
    for sub in partitions
        m = sub.model::JuMP.Model
        @inbounds for e in interzonal_branches
            l = get(sub.branch_global_to_local, e, 0)
            l == 0 && continue
            R, T = _rt_from_exprs(m, l, 1)
            acc +=
                (value(R) - get(r_consensus, e, 0.0))^2 +
                (value(T) - get(t_consensus, e, 0.0))^2
        end
    end
    return sqrt(acc)
end

function primal_residual_branch_pq_socp(
    partitions::Partitions{PartitionSOCP},
    p_flow_consensus::Dict{Int, Float64},
    q_flow_consensus::Dict{Int, Float64},
    interzonal_branches::Vector{Int},
)
    acc = 0.0
    for sub in partitions
        m = sub.model::JuMP.Model
        @inbounds for e in interzonal_branches
            l = get(sub.branch_global_to_local, e, 0)
            l == 0 && continue
            P, Q = _pq_flow_from_exprs(m, l, 1)
            acc +=
                (value(P) - get(p_flow_consensus, e, 0.0))^2 +
                (value(Q) - get(q_flow_consensus, e, 0.0))^2
        end
    end
    return sqrt(acc)
end

# ------------------------------------------------------------
# AC residuals eps primal
# ------------------------------------------------------------

function eps_primal_bus_ac_vec(
    partitions::Partitions{PartitionAC},
    theta_consensus_vec::Vector{Float64},
    v_mag_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 0
    for part in partitions
        dim += 2 * length(part.shared_buses_zone)
    end

    z_norm_sq = 0.0
    @inbounds for j in eachindex(theta_consensus_vec)
        z_norm_sq += theta_consensus_vec[j]^2 + v_mag_consensus_vec[j]^2
    end

    return sqrt(max(1, dim)) * eps_abs + eps_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_primal_branch_ac(
    partitions::Partitions{PartitionAC},
    p_flow_consensus::Dict{Int, Float64},
    q_flow_consensus::Dict{Int, Float64},
    interzonal_branches::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 0
    for part in partitions
        dim += 2 * length(_zone_cut_edges(part, interzonal_branches))
    end

    z_norm_sq = 0.0
    @inbounds for e in interzonal_branches
        z_norm_sq += get(p_flow_consensus, e, 0.0)^2 + get(q_flow_consensus, e, 0.0)^2
    end

    return sqrt(max(1, dim)) * eps_abs + eps_rel * max(sqrt(z_norm_sq), 1e-12)
end


# ------------------------------------------------------------
# DC tolerances: eps_primal / eps_dual (bus / branch)
# ------------------------------------------------------------

# function eps_primal_bus_dc(
#     partitions::Partitions{PartitionDC},
#     theta_consensus_vec::Vector{Float64},
#     shared_buses_sorted::Vector{Int},
#     shared_buses_to_pos::Dict{Int, Int},
#     eps_abs::Float64,
#     eps_rel::Float64,
# )
#     N = max(1, length(shared_buses_sorted))
#     norm_x_sq = 0.0

#     for part in partitions
#         thetabuf = snapshot_theta!(part)
#         pos = part.bus_local_index
#         @inbounds for b in part.shared_buses_zone
#             norm_x_sq += thetabuf[pos[b]]^2
#         end
#     end

#     norm_z_sq = sum(theta_consensus_vec[j]^2 for j in eachindex(theta_consensus_vec))
#     return sqrt(N) * eps_abs + eps_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
# end

function eps_primal_bus_dc(
    partitions::Partitions{PartitionDC},
    theta_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    norm_x_sq = 0.0

    for part in partitions
        thetabuf = snapshot_theta!(part)
        pos = part.bus_local_index
        n += length(part.shared_buses_zone)
        @inbounds for b in part.shared_buses_zone
            norm_x_sq += thetabuf[pos[b]]^2
        end
    end

    n = max(n, 1)
    norm_z_sq = sum(theta_consensus_vec[j]^2 for j in eachindex(theta_consensus_vec))
    return sqrt(n) * eps_abs + eps_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end


function eps_primal_branch_dc(
    partitions::Partitions{PartitionDC},
    flow_consensus::Dict{Int, Float64},
    interzonal_branches::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim_pri = 0
    for part in partitions
        dim_pri += length(_zone_cut_edges(part, interzonal_branches))
    end

    z_norm_sq = 0.0
    @inbounds for e in interzonal_branches
        z_norm_sq += get(flow_consensus, e, 0.0)^2
    end

    return sqrt(max(1, dim_pri)) * eps_abs + eps_rel * max(sqrt(z_norm_sq), 1e-12)
end

# ------------------------------------------------------------
# SOCP tolerances: eps_primal / eps_dual (bus / branch RT / branch PQ)
# ------------------------------------------------------------

function eps_primal_bus_socp_vec(
    partitions::Partitions{PartitionSOCP},
    v2_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    N = max(1, length(shared_buses_sorted))
    norm_x_sq = 0.0
    for sub in partitions
        v2buf = snapshot_v2!(sub)
        pos = sub.bus_local_index
        @inbounds for b in sub.shared_buses_zone
            norm_x_sq += v2buf[pos[b]]^2
        end
    end
    norm_z_sq = sum(v2_consensus_vec[j]^2 for j in eachindex(v2_consensus_vec))
    return sqrt(N) * eps_abs + eps_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

function eps_primal_branch_rt_socp(
    partitions::Partitions{PartitionSOCP},
    r_consensus::Dict{Int, Float64},
    t_consensus::Dict{Int, Float64},
    interzonal_branches::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 0
    for sub in partitions
        dim += 2 * length(_zone_cut_edges(sub, interzonal_branches))
    end
    z_norm_sq = 0.0
    @inbounds for e in interzonal_branches
        z_norm_sq += get(r_consensus, e, 0.0)^2 + get(t_consensus, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * eps_abs + eps_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_primal_branch_pq_socp(
    partitions::Partitions{PartitionSOCP},
    p_flow_consensus::Dict{Int, Float64},
    q_flow_consensus::Dict{Int, Float64},
    interzonal_branches::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 0
    for sub in partitions
        dim += 2 * length(_zone_cut_edges(sub, interzonal_branches))
    end
    z_norm_sq = 0.0
    @inbounds for e in interzonal_branches
        z_norm_sq += get(p_flow_consensus, e, 0.0)^2 + get(q_flow_consensus, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * eps_abs + eps_rel * max(sqrt(z_norm_sq), 1e-12)
end


# ------------------------------------------------------------
# AC residuals: dual (bus / branch)
# ------------------------------------------------------------

function dual_residual_bus_ac_vec(
    theta::Vector{Float64},
    vmag::Vector{Float64},
    theta_prev::Vector{Float64},
    vmag_prev::Vector{Float64},
    rho_theta::Float64,
    rho_v::Float64,
)
    sθ = 0.0
    sV = 0.0
    @inbounds for j in eachindex(theta)
        sθ += _wrap_angle_diff(theta[j] - theta_prev[j])^2
        sV += (vmag[j] - vmag_prev[j])^2
    end
    return sqrt((rho_theta^2) * sθ + (rho_v^2) * sV)
end

function dual_residual_bus_ac_components(
    theta_consensus_vec::Vector{Float64},
    v_mag_consensus_vec::Vector{Float64},
    theta_consensus_prev::Vector{Float64},
    v_mag_consensus_prev::Vector{Float64},
)
    stheta = 0.0
    sV = 0.0
    @inbounds for j in eachindex(theta_consensus_vec)
        stheta += _wrap_angle_diff(theta_consensus_vec[j] - theta_consensus_prev[j])^2
        sV += (v_mag_consensus_vec[j] - v_mag_consensus_prev[j])^2
    end
    return sqrt(stheta), sqrt(sV)
end

function dual_residual_branch_ac(
    p_flow::Dict{Int, Float64},
    q_flow::Dict{Int, Float64},
    p_prev::Dict{Int, Float64},
    q_prev::Dict{Int, Float64},
    interzonal_branches::Vector{Int},
    rho_p::Float64,
    rho_q::Float64,
)
    s = 0.0
    @inbounds for e in interzonal_branches
        dp = get(p_flow, e, 0.0) - get(p_prev, e, 0.0)
        dq = get(q_flow, e, 0.0) - get(q_prev, e, 0.0)
        s += (rho_p * dp)^2 + (rho_q * dq)^2
    end
    return sqrt(s)
end
# ------------------------------------------------------------
# DC residuals: dual (bus / branch)  [S ESCALADO POR ρ]
# ------------------------------------------------------------

dual_residual_bus_dc(
    theta::Vector{Float64},
    theta_prev::Vector{Float64},
    rho::Float64,
) = rho * sqrt(sum((theta[j] - theta_prev[j])^2 for j in eachindex(theta)))

dual_residual_branch_dc(
    flow::Dict{Int, Float64},
    flow_prev::Dict{Int, Float64},
    rho::Float64,
) = begin
    keys_all = union(keys(flow), keys(flow_prev))
    rho * sqrt(sum((get(flow, e, 0.0) - get(flow_prev, e, 0.0))^2 for e in keys_all))
end

# ------------------------------------------------------------
# SOCP residuals: dual (bus / branch RT / branch PQ)
# ------------------------------------------------------------

# ------------------------------------------------------------
# SOCP residuals: dual (bus / branch RT / branch PQ)  [ESCALADO POR ρ]
# ------------------------------------------------------------

function dual_residual_bus_socp_vec(
    v2_consensus_vec::Vector{Float64},
    v2_consensus_prev_vec::Vector{Float64},
    rho_v2::Float64,
)
    s = 0.0
    @inbounds for j in eachindex(v2_consensus_vec)
        s += (v2_consensus_vec[j] - v2_consensus_prev_vec[j])^2
    end
    return rho_v2 * sqrt(s)
end

function dual_residual_branch_rt_socp(
    r_consensus::Dict{Int, Float64},
    t_consensus::Dict{Int, Float64},
    r_consensus_prev::Dict{Int, Float64},
    t_consensus_prev::Dict{Int, Float64},
    interzonal_branches::Vector{Int},
    rho_r::Float64,
    rho_t::Float64,
)
    s = 0.0
    @inbounds for e in interzonal_branches
        dr = get(r_consensus, e, 0.0) - get(r_consensus_prev, e, 0.0)
        dt = get(t_consensus, e, 0.0) - get(t_consensus_prev, e, 0.0)
        s += (rho_r * dr)^2 + (rho_t * dt)^2
    end
    return sqrt(s)
end

function dual_residual_branch_pq_socp(
    p_flow_consensus::Dict{Int, Float64},
    q_flow_consensus::Dict{Int, Float64},
    p_flow_consensus_prev::Dict{Int, Float64},
    q_flow_consensus_prev::Dict{Int, Float64},
    interzonal_branches::Vector{Int},
    rhoP::Float64,
    rhoQ::Float64,
)
    s = 0.0
    @inbounds for e in interzonal_branches
        dp = get(p_flow_consensus, e, 0.0) - get(p_flow_consensus_prev, e, 0.0)
        dq = get(q_flow_consensus, e, 0.0) - get(q_flow_consensus_prev, e, 0.0)
        s += (rhoP * dp)^2 + (rhoQ * dq)^2
    end
    return sqrt(s)
end


# ------------------------------------------------------------
# AC tolerances: eps_dual baseado no dual "não escalado" y
# ------------------------------------------------------------

function eps_dual_bus_ac_unscaled(
    duals_theta::Dict{Int, Vector{Float64}},   # yθ (não escalado), por zona
    duals_v_mag::Dict{Int, Vector{Float64}},   # yV (não escalado), por zona
    partitions::Partitions{PartitionAC},
    zone_ids::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    y_norm_sq = 0.0

    @inbounds for zone_id in zone_ids
        part = get_partition(partitions, zone_id)
        yθ = duals_theta[zone_id]
        yV = duals_v_mag[zone_id]

        for b in part.shared_buses_zone
            j = shared_buses_to_pos[b]
            n += 2
            y_norm_sq += yθ[j]^2 + yV[j]^2
        end
    end

    n = max(n, 1)
    return sqrt(n) * eps_abs + eps_rel * sqrt(y_norm_sq)
end

function eps_dual_branch_ac_unscaled(
    duals_p_flow::Dict{Int, Dict{Int, Float64}},  # yp (não escalado)
    duals_q_flow::Dict{Int, Dict{Int, Float64}},  # yq (não escalado)
    zone_ids::Vector{Int},
    interzonal_branches::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    y_norm_sq = 0.0

    @inbounds for zone_id in zone_ids
        yp = duals_p_flow[zone_id]
        yq = duals_q_flow[zone_id]

        for e in interzonal_branches
            haskey(yp, e) || continue
            # assumindo simetria: se yp tem e, yq também deve ter (mas não dependo disso)
            n += 2
            y_norm_sq += get(yp, e, 0.0)^2 + get(yq, e, 0.0)^2
        end
    end

    n = max(n, 1)
    return sqrt(n) * eps_abs + eps_rel * sqrt(y_norm_sq)
end

# eps_dual baseado no dual "não escalado" y (bus)
function eps_dual_bus_dc_unscaled(
    dual_theta::Dict{Int, Vector{Float64}},   # y (não escalado)
    partitions::Partitions{PartitionDC},
    zone_ids::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    y_norm_sq = 0.0
    for zone_id in zone_ids
        part = get_partition(partitions, zone_id)
        y = dual_theta[zone_id]
        @inbounds for b in part.shared_buses_zone
            j = shared_buses_to_pos[b]
            n += 1
            y_norm_sq += y[j]^2
        end
    end
    n = max(n, 1)
    return sqrt(n) * eps_abs + eps_rel * sqrt(y_norm_sq)
end

# eps_dual baseado no dual "não escalado" y (branch)  <<< NOVO
function eps_dual_branch_dc_unscaled(
    duals_flow::Dict{Int, Dict{Int, Float64}},  # y (não escalado) por zona e por aresta global
    partitions::Partitions{PartitionDC},        # (mantido p/ simetria; pode não ser usado)
    zone_ids::Vector{Int},
    interzonal_branches::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    y_norm_sq = 0.0
    @inbounds for zone_id in zone_ids
        y = duals_flow[zone_id]
        for e in interzonal_branches
            haskey(y, e) || continue
            n += 1
            y_norm_sq += y[e]^2
        end
    end
    n = max(n, 1)
    return sqrt(n) * eps_abs + eps_rel * sqrt(y_norm_sq)
end

# ------------------------------------------------------------
# SOCP tolerances: eps_dual baseado no dual "não escalado" y
# ------------------------------------------------------------

function eps_dual_bus_socp_unscaled(
    duals_v2::Dict{Int, Vector{Float64}},
    partitions::Partitions{PartitionSOCP},
    zone_ids::Vector{Int},
    shared_buses_to_pos::Dict{Int, Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    y_norm_sq = 0.0

    @inbounds for zone_id in zone_ids
        part = get_partition(partitions, zone_id)
        y = duals_v2[zone_id]
        for b in part.shared_buses_zone
            j = shared_buses_to_pos[b]
            n += 1
            y_norm_sq += y[j]^2
        end
    end

    n = max(n, 1)
    return sqrt(n) * eps_abs + eps_rel * sqrt(y_norm_sq)
end

function eps_dual_branch_rt_socp_unscaled(
    duals_r::Dict{Int, Dict{Int, Float64}},   # y_r (não escalado)
    duals_t::Dict{Int, Dict{Int, Float64}},   # y_t (não escalado)
    zone_ids::Vector{Int},
    interzonal_branches::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    y_norm_sq = 0.0

    @inbounds for zone_id in zone_ids
        yr = duals_r[zone_id]
        yt = duals_t[zone_id]
        for e in interzonal_branches
            haskey(yr, e) || continue
            n += 2
            y_norm_sq += get(yr, e, 0.0)^2 + get(yt, e, 0.0)^2
        end
    end

    n = max(n, 1)
    return sqrt(n) * eps_abs + eps_rel * sqrt(y_norm_sq)
end

function eps_dual_branch_pq_socp_unscaled(
    dualsP::Dict{Int, Dict{Int, Float64}},
    dualsQ::Dict{Int, Dict{Int, Float64}},
    zone_ids::Vector{Int},
    interzonal_branches::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    y_norm_sq = 0.0

    @inbounds for zone_id in zone_ids
        yp = dualsP[zone_id]
        yq = dualsQ[zone_id]
        for e in interzonal_branches
            haskey(yp, e) || continue
            n += 2
            y_norm_sq += get(yp, e, 0.0)^2 + get(yq, e, 0.0)^2
        end
    end

    n = max(n, 1)
    return sqrt(n) * eps_abs + eps_rel * sqrt(y_norm_sq)
end