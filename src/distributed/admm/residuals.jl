# ============================================================
# Residuals and tolerances (ADMM)
# ============================================================

# ------------------------------------------------------------
# AC residuals: primal (bus / branch)
# ------------------------------------------------------------

function primal_residual_bus_ac_vec(
    partitions::Partitions{PartitionAC},
    theta_consensus_vec::Vector{Float64},
    vmag_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    acc = 0.0
    for part in partitions
        thetabuf = snapshot_theta_ac!(part)
        Vbuf = snapshot_vmag!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_local
            j = b2pos[b]
            acc +=
                _wrap_angle_diff(thetabuf[pos[b]] - theta_consensus_vec[j])^2 +
                (Vbuf[pos[b]] - vmag_consensus_vec[j])^2
        end
    end
    return sqrt(acc)
end

function primal_residual_bus_ac_components(
    partitions::Partitions{PartitionAC},
    theta_consensus_vec::Vector{Float64},
    vmag_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    rtheta = 0.0
    rV = 0.0
    for part in partitions
        thetabuf = snapshot_theta_ac!(part)
        Vbuf = snapshot_vmag!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_local
            j = b2pos[b]
            rtheta += _wrap_angle_diff(thetabuf[pos[b]] - theta_consensus_vec[j])^2
            rV += (Vbuf[pos[b]] - vmag_consensus_vec[j])^2
        end
    end
    return sqrt(rtheta), sqrt(rV)
end

function primal_residual_branch_ac(
    partitions::Partitions{PartitionAC},
    pflow_consensus::Dict{Int, Float64},
    qflow_consensus::Dict{Int, Float64},
    cut_edges::Vector{Int},
)
    acc = 0.0
    for part in partitions
        @inbounds for e in cut_edges
            l_local = get(part.edge_to_local, e, 0)
            l_local == 0 && continue

            P_loc, Q_loc = _pq_flow_from_exprs(part, l_local; t = 1)
            p = JuMP.value(P_loc)
            q = JuMP.value(Q_loc)

            if isfinite(p)
                acc += (p - get(pflow_consensus, e, 0.0))^2
            end
            if isfinite(q)
                acc += (q - get(qflow_consensus, e, 0.0))^2
            end
        end
    end
    return sqrt(acc)
end

# ------------------------------------------------------------
# AC residuals: dual (bus / branch)
# ------------------------------------------------------------

function dual_residual_bus_ac_vec(
    theta_consensus_vec::Vector{Float64},
    vmag_consensus_vec::Vector{Float64},
    theta_consensus_prev::Vector{Float64},
    vmag_consensus_prev::Vector{Float64},
)
    stheta = 0.0
    sV = 0.0
    @inbounds for j in eachindex(theta_consensus_vec)
        stheta += _wrap_angle_diff(theta_consensus_vec[j] - theta_consensus_prev[j])^2
        sV += (vmag_consensus_vec[j] - vmag_consensus_prev[j])^2
    end
    return sqrt(stheta + sV)
end

function dual_residual_bus_ac_components(
    theta_consensus_vec::Vector{Float64},
    vmag_consensus_vec::Vector{Float64},
    theta_consensus_prev::Vector{Float64},
    vmag_consensus_prev::Vector{Float64},
)
    stheta = 0.0
    sV = 0.0
    @inbounds for j in eachindex(theta_consensus_vec)
        stheta += _wrap_angle_diff(theta_consensus_vec[j] - theta_consensus_prev[j])^2
        sV += (vmag_consensus_vec[j] - vmag_consensus_prev[j])^2
    end
    return sqrt(stheta), sqrt(sV)
end

function dual_residual_branch_ac(
    pflow_consensus::Dict{Int, Float64},
    qflow_consensus::Dict{Int, Float64},
    pflow_consensus_prev::Dict{Int, Float64},
    qflow_consensus_prev::Dict{Int, Float64},
    cut_edges::Vector{Int},
)
    s = 0.0
    @inbounds for e in cut_edges
        s += (get(pflow_consensus, e, 0.0) - get(pflow_consensus_prev, e, 0.0))^2
        s += (get(qflow_consensus, e, 0.0) - get(qflow_consensus_prev, e, 0.0))^2
    end
    return sqrt(s)
end

# ------------------------------------------------------------
# AC tolerances: eps_primal / eps_dual (bus / branch)
# ------------------------------------------------------------

function eps_primal_bus_ac_vec(
    partitions::Partitions{PartitionAC},
    theta_consensus_vec::Vector{Float64},
    vmag_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 0
    for part in partitions
        dim += 2 * length(part.shared_buses_local)
    end

    z_norm_sq = 0.0
    @inbounds for j in eachindex(theta_consensus_vec)
        z_norm_sq += theta_consensus_vec[j]^2 + vmag_consensus_vec[j]^2
    end

    return sqrt(max(1, dim)) * eps_abs + eps_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_primal_branch_ac(
    partitions::Partitions{PartitionAC},
    pflow_consensus::Dict{Int, Float64},
    qflow_consensus::Dict{Int, Float64},
    cut_edges::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 0
    for part in partitions
        # mais consistente com o residual: só edges presentes (cut-local)
        dim += 2 * length(_zone_cut_edges(part, cut_edges))
    end

    z_norm_sq = 0.0
    @inbounds for e in cut_edges
        z_norm_sq += get(pflow_consensus, e, 0.0)^2 + get(qflow_consensus, e, 0.0)^2
    end

    return sqrt(max(1, dim)) * eps_abs + eps_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_bus_ac_vec(
    theta_consensus_vec::Vector{Float64},
    vmag_consensus_vec::Vector{Float64},
    theta_consensus_prev::Vector{Float64},
    vmag_consensus_prev::Vector{Float64},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 2 * max(1, length(theta_consensus_vec))
    s = dual_residual_bus_ac_vec(
        theta_consensus_vec,
        vmag_consensus_vec,
        theta_consensus_prev,
        vmag_consensus_prev,
    )
    return sqrt(dim) * eps_abs + eps_rel * s
end

function eps_dual_branch_ac(
    pflow_consensus::Dict{Int, Float64},
    qflow_consensus::Dict{Int, Float64},
    pflow_consensus_prev::Dict{Int, Float64},
    qflow_consensus_prev::Dict{Int, Float64},
    cut_edges::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_ac(
        pflow_consensus,
        qflow_consensus,
        pflow_consensus_prev,
        qflow_consensus_prev,
        cut_edges,
    )
    return sqrt(dim) * eps_abs + eps_rel * s
end

# ------------------------------------------------------------
# DC residuals: primal (bus / branch)
# ------------------------------------------------------------

function primal_residual_bus_dc(
    partitions::Partitions{PartitionDC},
    theta_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int}, # compat (não usado)
    b2pos::Dict{Int, Int},
)
    acc = 0.0
    for part in partitions
        thetabuf = snapshot_theta!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_local
            j = b2pos[b]
            acc += (thetabuf[pos[b]] - theta_consensus_vec[j])^2
        end
    end
    return sqrt(acc)
end

function primal_residual_branch_dc(
    partitions::Partitions{PartitionDC},
    flow_consensus::Dict{Int, Float64},
    cut_edges::Vector{Int},
)
    acc = 0.0
    for part in partitions
        m = part.model::JuMP.Model
        @inbounds for e in cut_edges
            l_local = get(part.edge_to_local, e, 0)
            l_local == 0 && continue
            acc += (_dc_flow_val_model(m, part, l_local) - get(flow_consensus, e, 0.0))^2
        end
    end
    return sqrt(acc)
end

# ------------------------------------------------------------
# DC residuals: dual (bus / branch)
# ------------------------------------------------------------

dual_residual_bus_dc(
    theta_consensus_vec::Vector{Float64},
    theta_consensus_prev_vec::Vector{Float64},
) = sqrt(
    sum(
        (theta_consensus_vec[j] - theta_consensus_prev_vec[j])^2 for
        j in eachindex(theta_consensus_vec)
    ),
)

dual_residual_branch_dc(
    flow_consensus::Dict{Int, Float64},
    flow_consensus_prev::Dict{Int, Float64},
) = sqrt(
    sum(
        (get(flow_consensus, e, 0.0) - get(flow_consensus_prev, e, 0.0))^2 for
        e in keys(flow_consensus)
    ),
)

# ------------------------------------------------------------
# DC tolerances: eps_primal / eps_dual (bus / branch)
# ------------------------------------------------------------

function eps_primal_bus_dc(
    partitions::Partitions{PartitionDC},
    theta_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    N = max(1, length(shared_buses_sorted))
    norm_x_sq = 0.0

    for part in partitions
        thetabuf = snapshot_theta!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_local
            norm_x_sq += thetabuf[pos[b]]^2
        end
    end

    norm_z_sq = sum(theta_consensus_vec[j]^2 for j in eachindex(theta_consensus_vec))
    return sqrt(N) * eps_abs + eps_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

function eps_primal_branch_dc(
    partitions::Partitions{PartitionDC},
    flow_consensus::Dict{Int, Float64},
    cut_edges::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim_pri = 0
    for part in partitions
        dim_pri += length(_zone_cut_edges(part, cut_edges))
    end

    z_norm_sq = 0.0
    @inbounds for e in cut_edges
        z_norm_sq += get(flow_consensus, e, 0.0)^2
    end

    return sqrt(max(1, dim_pri)) * eps_abs + eps_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_bus_dc(
    theta_consensus_vec::Vector{Float64},
    theta_consensus_prev_vec::Vector{Float64},
    eps_abs::Float64,
    eps_rel::Float64,
)
    N = max(1, length(theta_consensus_vec))
    s = dual_residual_bus_dc(theta_consensus_vec, theta_consensus_prev_vec)
    return sqrt(N) * eps_abs + eps_rel * s
end

function eps_dual_branch_dc(
    flow_consensus::Dict{Int, Float64},
    flow_consensus_prev::Dict{Int, Float64},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim_du = max(1, length(flow_consensus))
    s = dual_residual_branch_dc(flow_consensus, flow_consensus_prev)
    return sqrt(dim_du) * eps_abs + eps_rel * s
end

# ------------------------------------------------------------
# SOCP residuals: primal (bus / branch RT / branch PQ)
# ------------------------------------------------------------

function primal_residual_bus_socp_vec(
    partitions::Partitions{PartitionSOCP},
    v2_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int}, # compat (não usado)
    b2pos::Dict{Int, Int},
)
    acc = 0.0
    for sub in partitions
        v2buf = snapshot_v2!(sub)
        pos = sub.bus_local_index
        @inbounds for b in sub.shared_buses_local
            j = b2pos[b]
            acc += (v2buf[pos[b]] - v2_consensus_vec[j])^2
        end
    end
    return sqrt(acc)
end

function primal_residual_branch_rt_socp(
    partitions::Partitions{PartitionSOCP},
    rflow_consensus::Dict{Int, Float64},
    tflow_consensus::Dict{Int, Float64},
    cut_edges::Vector{Int},
)
    acc = 0.0
    for sub in partitions
        m = sub.model::JuMP.Model
        @inbounds for e in cut_edges
            l = get(sub.edge_to_local, e, 0)
            l == 0 && continue
            R, T = _rt_from_exprs(m, l, 1)
            acc +=
                (value(R) - get(rflow_consensus, e, 0.0))^2 +
                (value(T) - get(tflow_consensus, e, 0.0))^2
        end
    end
    return sqrt(acc)
end

function primal_residual_branch_pq_socp(
    partitions::Partitions{PartitionSOCP},
    pflow_consensus::Dict{Int, Float64},
    qflow_consensus::Dict{Int, Float64},
    cut_edges::Vector{Int},
)
    acc = 0.0
    for sub in partitions
        m = sub.model::JuMP.Model
        @inbounds for e in cut_edges
            l = get(sub.edge_to_local, e, 0)
            l == 0 && continue
            P, Q = _pq_flow_from_exprs(m, l, 1)
            acc +=
                (value(P) - get(pflow_consensus, e, 0.0))^2 +
                (value(Q) - get(qflow_consensus, e, 0.0))^2
        end
    end
    return sqrt(acc)
end

# ------------------------------------------------------------
# SOCP residuals: dual (bus / branch RT / branch PQ)
# ------------------------------------------------------------

dual_residual_bus_socp_vec(v2_consensus_vec, v2_consensus_prev_vec) =
    sqrt(
        sum(
            (v2_consensus_vec[j] - v2_consensus_prev_vec[j])^2 for
            j in eachindex(v2_consensus_vec)
        ),
    )

dual_residual_branch_rt_socp(
    rflow_consensus,
    tflow_consensus,
    rflow_consensus_prev,
    tflow_consensus_prev,
    cut_edges,
) =
    sqrt(
        sum(
            (get(rflow_consensus, e, 0.0) - get(rflow_consensus_prev, e, 0.0))^2 +
            (get(tflow_consensus, e, 0.0) - get(tflow_consensus_prev, e, 0.0))^2
            for e in cut_edges
        ),
    )

dual_residual_branch_pq_socp(
    pflow_consensus,
    qflow_consensus,
    pflow_consensus_prev,
    qflow_consensus_prev,
    cut_edges,
) =
    sqrt(
        sum(
            (get(pflow_consensus, e, 0.0) - get(pflow_consensus_prev, e, 0.0))^2 +
            (get(qflow_consensus, e, 0.0) - get(qflow_consensus_prev, e, 0.0))^2
            for e in cut_edges
        ),
    )

# ------------------------------------------------------------
# SOCP tolerances: eps_primal / eps_dual (bus / branch RT / branch PQ)
# ------------------------------------------------------------

function eps_primal_bus_socp_vec(
    partitions::Partitions{PartitionSOCP},
    v2_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    N = max(1, length(shared_buses_sorted))
    norm_x_sq = 0.0
    for sub in partitions
        v2buf = snapshot_v2!(sub)
        pos = sub.bus_local_index
        @inbounds for b in sub.shared_buses_local
            norm_x_sq += v2buf[pos[b]]^2
        end
    end
    norm_z_sq = sum(v2_consensus_vec[j]^2 for j in eachindex(v2_consensus_vec))
    return sqrt(N) * eps_abs + eps_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

function eps_primal_branch_rt_socp(
    partitions::Partitions{PartitionSOCP},
    rflow_consensus::Dict{Int, Float64},
    tflow_consensus::Dict{Int, Float64},
    cut_edges::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 0
    for sub in partitions
        dim += 2 * length(_zone_cut_edges(sub, cut_edges))
    end
    z_norm_sq = 0.0
    @inbounds for e in cut_edges
        z_norm_sq += get(rflow_consensus, e, 0.0)^2 + get(tflow_consensus, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * eps_abs + eps_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_primal_branch_pq_socp(
    partitions::Partitions{PartitionSOCP},
    pflow_consensus::Dict{Int, Float64},
    qflow_consensus::Dict{Int, Float64},
    cut_edges::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 0
    for sub in partitions
        dim += 2 * length(_zone_cut_edges(sub, cut_edges))
    end
    z_norm_sq = 0.0
    @inbounds for e in cut_edges
        z_norm_sq += get(pflow_consensus, e, 0.0)^2 + get(qflow_consensus, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * eps_abs + eps_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_bus_socp_vec(v2_consensus_vec, v2_consensus_prev_vec, eps_abs, eps_rel)
    N = max(1, length(v2_consensus_vec))
    s = dual_residual_bus_socp_vec(v2_consensus_vec, v2_consensus_prev_vec)
    return sqrt(N) * eps_abs + eps_rel * s
end

function eps_dual_branch_rt_socp(
    rflow_consensus,
    tflow_consensus,
    rflow_consensus_prev,
    tflow_consensus_prev,
    cut_edges,
    eps_abs,
    eps_rel,
)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_rt_socp(
        rflow_consensus,
        tflow_consensus,
        rflow_consensus_prev,
        tflow_consensus_prev,
        cut_edges,
    )
    return sqrt(dim) * eps_abs + eps_rel * s
end

function eps_dual_branch_pq_socp(
    pflow_consensus,
    qflow_consensus,
    pflow_consensus_prev,
    qflow_consensus_prev,
    cut_edges,
    eps_abs,
    eps_rel,
)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_pq_socp(
        pflow_consensus,
        qflow_consensus,
        pflow_consensus_prev,
        qflow_consensus_prev,
        cut_edges,
    )
    return sqrt(dim) * eps_abs + eps_rel * s
end

# ------------------------------------------------------------
# Misc utilities
# ------------------------------------------------------------

@inline function gather_zone_lambda!(
    lambda_vec_zone::Vector{Float64},
    duals_dense_zone::Vector{Float64},
)
    lambda_vec_zone .= duals_dense_zone
    return lambda_vec_zone
end
