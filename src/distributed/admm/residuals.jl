
function primal_residual_bus_ac_vec(
    partitions::Dict{Int, PartitionAC},
    theta_consensus_vec::Vector{Float64},
    vmag_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    acc = 0.0
    for (_, part) in partitions
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
    partitions::Dict{Int, PartitionAC},
    theta_consensus_vec::Vector{Float64},
    vmag_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    rtheta = 0.0
    rV = 0.0
    for (_, part) in partitions
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
    partitions::Dict{Int, PartitionAC},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    cut_edges::Vector{Int},
)
    acc = 0.0
    for (_, part) in partitions
        m = part.model::JuMP.Model
        for e in part.zone_edges
            l_local = get(part.edge_to_local, e, nothing)
            l_local === nothing && continue
            P_loc, Q_loc = _ac_flow_from_exprs_existing(part, l_local)
            p = JuMP.value(P_loc)
            q = JuMP.value(Q_loc)
            if isfinite(p)
                acc += (p - get(zP, e, 0.0))^2
            end
            if isfinite(q)
                acc += (q - get(zQ, e, 0.0))^2
            end
        end
    end
    return sqrt(acc)
end


function primal_residual_bus_dc(
    partitions::Dict{Int, PartitionDC},
    theta_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    acc = 0.0
    for (_, part) in partitions
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
    partitions::Dict{Int, PartitionDC},
    z_flow::Dict{Int, Float64},
    cut_edges::Vector{Int},
)
    acc = 0.0
    for (_, part) in partitions
        m = part.model::JuMP.Model
        @inbounds for e in cut_edges
            l_local = get(part.edge_to_local, e, nothing)
            if l_local !== nothing
                acc += (_dc_flow_val_model(m, part, l_local) - z_flow[e])^2
            end
        end
    end
    return sqrt(acc)
end


# ----------------- Resíduos / tolerâncias (bus em forma densa) -----------------
function primal_residual_bus_socp_vec(partitions, v2bar_vec, shared_buses_sorted, b2pos)
    acc = 0.0
    for (_, sub) in partitions
        v2buf = snapshot_v2!(sub)
        pos = sub.bus_local_index
        @inbounds for b in sub.shared_buses_local
            j = b2pos[b]
            acc += (v2buf[pos[b]] - v2bar_vec[j])^2
        end
    end
    return sqrt(acc)
end


function primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
    acc = 0.0
    for (_, sub) in partitions
        m = sub.model::JuMP.Model
        @inbounds for e in cut_edges
            l = get(sub.edge_to_local, e, nothing)
            l === nothing && continue
            P, Q = _socp_flow_from_exprs_existing(m, l)
            acc += (value(P) - zP[e])^2 + (value(Q) - zQ[e])^2
        end
    end
    return sqrt(acc)
end


# ----------------- Resíduos / tolerâncias (edges) -----------------
function primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
    acc = 0.0
    for (_, sub) in partitions
        m = sub.model::JuMP.Model
        @inbounds for e in cut_edges
            l = get(sub.edge_to_local, e, nothing)
            l === nothing && continue
            R, T = _rt_from_exprs_socp(m, l)
            acc += (value(R) - zR[e])^2 + (value(T) - zT[e])^2
        end
    end
    return sqrt(acc)
end

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
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    zP_prev::Dict{Int, Float64},
    zQ_prev::Dict{Int, Float64},
    cut_edges::Vector{Int},
)
    s = 0.0
    for e in cut_edges
        s += (get(zP, e, 0.0) - get(zP_prev, e, 0.0))^2
        s += (get(zQ, e, 0.0) - get(zQ_prev, e, 0.0))^2
    end
    return sqrt(s)
end

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
    z_flow::Dict{Int, Float64},
    z_flow_prev::Dict{Int, Float64},
) = sqrt(sum((z_flow[e] - get(z_flow_prev, e, 0.0))^2 for e in keys(z_flow)))

dual_residual_bus_socp_vec(v2bar_vec, v2bar_prev_vec) =
    sqrt(sum((v2bar_vec[j] - v2bar_prev_vec[j])^2 for j in eachindex(v2bar_vec)))

    dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges) =
    sqrt(
        sum(
            (get(zR, e, 0.0) - get(zR_prev, e, 0.0))^2 +
            (get(zT, e, 0.0) - get(zT_prev, e, 0.0))^2
            for e in cut_edges
        ),
    )

dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges) =
    sqrt(
        sum(
            (get(zP, e, 0.0) - get(zP_prev, e, 0.0))^2 +
            (get(zQ, e, 0.0) - get(zQ_prev, e, 0.0))^2
            for e in cut_edges
        ),
    )



function eps_primal_bus_ac_vec(
    partitions::Dict{Int, PartitionAC},
    theta_consensus_vec::Vector{Float64},
    vmag_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 0
    for (_, part) in partitions
        dim += 2 * length(part.shared_buses_local)
    end
    z_norm_sq = 0.0
    @inbounds for j in eachindex(shared_buses_sorted)
        z_norm_sq += theta_consensus_vec[j]^2 + vmag_consensus_vec[j]^2
    end
    return sqrt(max(1, dim)) * eps_abs +
           eps_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_primal_branch_ac(
    partitions::Dict{Int, PartitionAC},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    cut_edges::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 0
    for (_, part) in partitions
        dim += 2 * length(part.zone_edges)
    end
    z_norm_sq = 0.0
    for e in cut_edges
        z_norm_sq += get(zP, e, 0.0)^2 + get(zQ, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * eps_abs +
           eps_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_primal_bus_dc(
    partitions::Dict{Int, PartitionDC},
    theta_consensus_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    N = max(1, length(shared_buses_sorted))
    norm_x_sq = 0.0
    for (_, part) in partitions
        thetabuf = snapshot_theta!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_local
            norm_x_sq += (thetabuf[pos[b]])^2
        end
    end
    norm_z_sq = sum((theta_consensus_vec[j])^2 for j in eachindex(theta_consensus_vec))

    return sqrt(N) * eps_abs + eps_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

function eps_primal_branch_dc(
    partitions::Dict{Int, PartitionDC},
    z_flow::Dict{Int, Float64},
    cut_edges::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim_pri = sum(length(_zone_cut_edges(part, cut_edges)) for (_, part) in partitions)
    return sqrt(max(1, dim_pri)) * eps_abs +
           eps_rel * max(sqrt(sum(z_flow[e]^2 for e in cut_edges)), 1e-12)
end

function eps_dual_branch_dc(
    z_flow::Dict{Int, Float64},
    z_flow_prev::Dict{Int, Float64},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim_du = max(1, length(z_flow))
    return sqrt(dim_du) * eps_abs +
           eps_rel *
           sqrt(sum((z_flow[e] - get(z_flow_prev, e, 0.0))^2 for e in keys(z_flow)))
end

##########################
# DUAIS (FLOW)           #
##########################

function eps_primal_bus_socp_vec(
    partitions,
    v2bar_vec,
    shared_buses_sorted,
    b2pos,
    eps_abs,
    eps_rel,
)
    N = max(1, length(shared_buses_sorted))
    norm_x_sq = 0.0
    for (_, sub) in partitions
        v2buf = snapshot_v2!(sub)
        pos = sub.bus_local_index
        @inbounds for b in sub.shared_buses_local
            norm_x_sq += (v2buf[pos[b]])^2
        end
    end
    norm_z_sq = sum((v2bar_vec[j])^2 for j in eachindex(v2bar_vec))
    return sqrt(N) * eps_abs + eps_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

function eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, eps_abs, eps_rel)
    dim = 0
    for (_, sub) in partitions
        dim += 2 * length(_zone_cut_edges(sub, cut_edges))
    end
    z_norm_sq = 0.0
    @inbounds for e in cut_edges
        z_norm_sq += get(zR, e, 0.0)^2 + get(zT, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * eps_abs + eps_rel * max(sqrt(z_norm_sq), 1e-12)
end


function eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, eps_abs, eps_rel)
    dim = 0
    for (_, sub) in partitions
        dim += 2 * length(_zone_cut_edges(sub, cut_edges))
    end
    z_norm_sq = 0.0
    @inbounds for e in cut_edges
        z_norm_sq += get(zP, e, 0.0)^2 + get(zQ, e, 0.0)^2
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

#############################
# RESÍDUOS BRANCH AC        #
#############################


function eps_dual_branch_ac(
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    zP_prev::Dict{Int, Float64},
    zQ_prev::Dict{Int, Float64},
    cut_edges::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
    return sqrt(dim) * eps_abs + eps_rel * s
end


#############################
# DUAIS EM theta                #
#############################
@inline function gather_zone_lambda!(
    lambda_vec_zone::Vector{Float64},
    duals_dense_zone::Vector{Float64},
)
    lambda_vec_zone .= duals_dense_zone
    return lambda_vec_zone
end

#############################
# RESÍDUOS (BUS)            #
#############################
function eps_dual_bus_dc(
    theta_consensus_vec::Vector{Float64},
    theta_consensus_prev_vec::Vector{Float64},
    eps_abs::Float64,
    eps_rel::Float64,
)
    N = max(1, length(theta_consensus_vec))
    return sqrt(N) * eps_abs +
           eps_rel * sqrt(
        sum(
            (theta_consensus_vec[j] - theta_consensus_prev_vec[j])^2 for
            j in eachindex(theta_consensus_vec)
        ),
    )
end

#############################
# RESÍDUOS (BRANCH)         #
#############################



function eps_dual_bus_socp_vec(v2bar_vec, v2bar_prev_vec, eps_abs, eps_rel)
    N = max(1, length(v2bar_vec))
    return sqrt(N) * eps_abs +
           eps_rel * sqrt(sum((v2bar_vec[j] - v2bar_prev_vec[j])^2
                    for j in eachindex(v2bar_vec)))
end


function eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, eps_abs, eps_rel)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
    return sqrt(dim) * eps_abs + eps_rel * s
end


function eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, eps_abs, eps_rel)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
    return sqrt(dim) * eps_abs + eps_rel * s
end