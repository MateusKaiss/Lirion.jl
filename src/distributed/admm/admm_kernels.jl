@inline _wrap_angle_diff(Δ::Float64) = atan(sin(Δ), cos(Δ))


function compute_ac_bus_consensus_vec!(
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    partitions::Dict{Int, PartitionAC},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    fill!(θ̄_vec, 0.0)
    fill!(V̄_vec, 0.0)
    sums_sin = zeros(Float64, length(shared_buses_sorted))
    sums_cos = zeros(Float64, length(shared_buses_sorted))
    counts = zeros(Int, length(shared_buses_sorted))

    for (_, part) in partitions
        θs = snapshot_theta_ac_shared!(part)
        Vs = snapshot_vmag_shared!(part)
        lbs = part.shared_buses_local
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            s, c = sincos(θs[i])
            sums_sin[j] += s
            sums_cos[j] += c
            V̄_vec[j] += Vs[i]
            counts[j] += 1
        end
    end

    @inbounds for j in eachindex(θ̄_vec)
        if counts[j] != 0
            θ̄_vec[j] = atan(sums_sin[j], sums_cos[j])
            V̄_vec[j] /= counts[j]
        end
    end
end

function compute_ac_branch_consensus!(
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    partitions::Dict{Int, PartitionAC},
    cut_edges::Vector{Int},
)
    sumsP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    sumsQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    cntP = Dict{Int, Int}(e => 0 for e in cut_edges)
    cntQ = Dict{Int, Int}(e => 0 for e in cut_edges)

    for (_, part) in partitions
        m = part.model::JuMP.Model
        for e in part.zone_edges
            l_local = get(part.edge_to_local, e, nothing)
            l_local === nothing && continue
            P_loc, Q_loc = _ac_flow_from_exprs_existing(part, l_local)
            pval = JuMP.value(P_loc)
            qval = JuMP.value(Q_loc)
            if isfinite(pval)
                sumsP[e] += pval
                cntP[e] += 1
            end
            if isfinite(qval)
                sumsQ[e] += qval
                cntQ[e] += 1
            end
        end
    end

    for e in cut_edges
        if cntP[e] > 0
            zP[e] = sumsP[e] / cntP[e]
        end
        if cntQ[e] > 0
            zQ[e] = sumsQ[e] / cntQ[e]
        end
    end
    return nothing
end

@inline function update_duals_bus_ac_vec!(
    dualsθ_vec_zone::Vector{Float64},
    dualsV_vec_zone::Vector{Float64},
    part::PartitionAC,
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    rhoθ::Float64,
    rhoV::Float64,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int};
    α::Union{Nothing, Float64} = nothing,
    θ̄_prev::Union{Nothing, Vector{Float64}} = nothing,
    V̄_prev::Union{Nothing, Vector{Float64}} = nothing,
)
    θs = snapshot_theta_ac_shared!(part)
    Vs = snapshot_vmag_shared!(part)
    lbs = part.shared_buses_local

    if α === nothing
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            dθ = _wrap_angle_diff(θs[i] - θ̄_vec[j])
            dV = (Vs[i] - V̄_vec[j])
            dualsθ_vec_zone[j] += rhoθ * dθ
            dualsV_vec_zone[j] += rhoV * dV
        end
    else
        αv = something(α, 1.0)
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            θrel = _wrap_angle_diff(αv * θs[i] + (1 - αv) * θ̄_prev[j] - θ̄_vec[j])
            Vrel = αv * Vs[i] + (1 - αv) * V̄_prev[j] - V̄_vec[j]
            dualsθ_vec_zone[j] += rhoθ * θrel
            dualsV_vec_zone[j] += rhoV * Vrel
        end
    end
    return nothing
end

function update_duals_branch_ac!(
    dualsP::Dict{Int, Dict{Int, Float64}},
    dualsQ::Dict{Int, Dict{Int, Float64}},
    partitions::Dict{Int, PartitionAC},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    cut_edges::Vector{Int},
    rhoP::Float64,
    rhoQ::Float64;
    α::Union{Nothing, Float64} = nothing,
    zP_prev::Union{Nothing, Dict{Int, Float64}} = nothing,
    zQ_prev::Union{Nothing, Dict{Int, Float64}} = nothing,
)
    for (z, part) in partitions
        m = part.model::JuMP.Model
        lambdaPz = get!(dualsP, z, Dict{Int, Float64}())
        lambdaQz = get!(dualsQ, z, Dict{Int, Float64}())

        for e in part.zone_edges
            l_local = get(part.edge_to_local, e, nothing)
            l_local === nothing && continue
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

#############################
# RESÍDUOS BUS AC           #
#############################

function primal_residual_bus_ac_vec(
    partitions::Dict{Int, PartitionAC},
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    acc = 0.0
    for (_, part) in partitions
        θbuf = snapshot_theta_ac!(part)
        Vbuf = snapshot_vmag!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_local
            j = b2pos[b]
            acc +=
                _wrap_angle_diff(θbuf[pos[b]] - θ̄_vec[j])^2 +
                (Vbuf[pos[b]] - V̄_vec[j])^2
        end
    end
    return sqrt(acc)
end

function primal_residual_bus_ac_components(
    partitions::Dict{Int, PartitionAC},
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    rθ = 0.0
    rV = 0.0
    for (_, part) in partitions
        θbuf = snapshot_theta_ac!(part)
        Vbuf = snapshot_vmag!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_local
            j = b2pos[b]
            rθ += _wrap_angle_diff(θbuf[pos[b]] - θ̄_vec[j])^2
            rV += (Vbuf[pos[b]] - V̄_vec[j])^2
        end
    end
    return sqrt(rθ), sqrt(rV)
end

function eps_primal_bus_ac_vec(
    partitions::Dict{Int, PartitionAC},
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    dim = 0
    for (_, part) in partitions
        dim += 2 * length(part.shared_buses_local)
    end
    z_norm_sq = 0.0
    @inbounds for j in eachindex(shared_buses_sorted)
        z_norm_sq += θ̄_vec[j]^2 + V̄_vec[j]^2
    end
    return sqrt(max(1, dim)) * ϵ_abs +
           ϵ_rel * max(sqrt(z_norm_sq), 1e-12)
end

function dual_residual_bus_ac_vec(
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    θ̄_prev::Vector{Float64},
    V̄_prev::Vector{Float64},
)
    sθ = 0.0
    sV = 0.0
    @inbounds for j in eachindex(θ̄_vec)
        sθ += _wrap_angle_diff(θ̄_vec[j] - θ̄_prev[j])^2
        sV += (V̄_vec[j] - V̄_prev[j])^2
    end
    return sqrt(sθ + sV)
end

function dual_residual_bus_ac_components(
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    θ̄_prev::Vector{Float64},
    V̄_prev::Vector{Float64},
)
    sθ = 0.0
    sV = 0.0
    @inbounds for j in eachindex(θ̄_vec)
        sθ += _wrap_angle_diff(θ̄_vec[j] - θ̄_prev[j])^2
        sV += (V̄_vec[j] - V̄_prev[j])^2
    end
    return sqrt(sθ), sqrt(sV)
end

function eps_dual_bus_ac_vec(
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    θ̄_prev::Vector{Float64},
    V̄_prev::Vector{Float64},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    dim = 2 * max(1, length(θ̄_vec))
    s = dual_residual_bus_ac_vec(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev)
    return sqrt(dim) * ϵ_abs + ϵ_rel * s
end

#############################
# RESÍDUOS BRANCH AC        #
#############################

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

function eps_primal_branch_ac(
    partitions::Dict{Int, PartitionAC},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    cut_edges::Vector{Int},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    dim = 0
    for (_, part) in partitions
        dim += 2 * length(part.zone_edges)
    end
    z_norm_sq = 0.0
    for e in cut_edges
        z_norm_sq += get(zP, e, 0.0)^2 + get(zQ, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * ϵ_abs +
           ϵ_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_branch_ac(
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    zP_prev::Dict{Int, Float64},
    zQ_prev::Dict{Int, Float64},
    cut_edges::Vector{Int},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
    return sqrt(dim) * ϵ_abs + ϵ_rel * s
end


function compute_theta_consensus_vec!(
    θ̄_vec::Vector{Float64},
    partitions::Dict{Int, PartitionDC},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    fill!(θ̄_vec, 0.0)
    counts = zeros(Int, length(shared_buses_sorted))

    # média das θ nas barras compartilhadas
    for (_, part) in partitions
        θs = snapshot_theta_shared!(part)
        lbs = part.shared_buses_local
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            θ̄_vec[j] += θs[i]
            counts[j] += 1
        end
    end

    @inbounds for j in eachindex(θ̄_vec)
        counts[j] != 0 && (θ̄_vec[j] /= counts[j])
    end
end

function compute_z_flow_consensus!(
    z_flow::Dict{Int, Float64},
    partitions::Dict{Int, PartitionDC},
    cut_edges::Vector{Int},
)
    sums = Dict{Int, Float64}()
    cnts = Dict{Int, Int}()

    # média dos fluxos p_e em cada ramo de corte
    for (_, part) in partitions
        m = part.model::JuMP.Model
        @inbounds for e in cut_edges
            l_local = get(part.edge_to_local, e, nothing)
            if l_local !== nothing
                p = _dc_flow_val_model(m, part, l_local)
                sums[e] = get(sums, e, 0.0) + p
                cnts[e] = get(cnts, e, 0) + 1
            end
        end
    end

    @inbounds for e in cut_edges
        n = get(cnts, e, 0)
        n > 0 && (z_flow[e] = sums[e] / n)
    end
end

#############################
# DUAIS EM θ                #
#############################

@inline function update_duals_theta_vec!(
    dual_vec_zone::Vector{Float64},
    part::PartitionDC,
    θ̄_vec::Vector{Float64},
    rho::Float64,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    θs = snapshot_theta_shared!(part)
    lbs = part.shared_buses_local

    @inbounds for i in eachindex(lbs)
        j = b2pos[lbs[i]]
        dual_vec_zone[j] += rho * (θs[i] - θ̄_vec[j])
    end
    return nothing
end

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

function primal_residual_bus_dc(
    partitions::Dict{Int, PartitionDC},
    θ̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    acc = 0.0
    for (_, part) in partitions
        θbuf = snapshot_theta!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_local
            j = b2pos[b]
            acc += (θbuf[pos[b]] - θ̄_vec[j])^2
        end
    end
    return sqrt(acc)
end

dual_residual_bus_dc(
    θ̄_vec::Vector{Float64},
    θ̄_prev_vec::Vector{Float64},
) = sqrt(sum((θ̄_vec[j] - θ̄_prev_vec[j])^2 for j in eachindex(θ̄_vec)))

function eps_primal_bus_dc(
    partitions::Dict{Int, PartitionDC},
    θ̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    N = max(1, length(shared_buses_sorted))
    norm_x_sq = 0.0
    for (_, part) in partitions
        θbuf = snapshot_theta!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_local
            norm_x_sq += (θbuf[pos[b]])^2
        end
    end
    norm_z_sq = sum((θ̄_vec[j])^2 for j in eachindex(θ̄_vec))

    return sqrt(N) * ϵ_abs + ϵ_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

function eps_dual_bus_dc(
    θ̄_vec::Vector{Float64},
    θ̄_prev_vec::Vector{Float64},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    N = max(1, length(θ̄_vec))
    return sqrt(N) * ϵ_abs +
           ϵ_rel * sqrt(sum((θ̄_vec[j] - θ̄_prev_vec[j])^2 for j in eachindex(θ̄_vec)))
end

#############################
# RESÍDUOS (BRANCH)         #
#############################


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

dual_residual_branch_dc(
    z_flow::Dict{Int, Float64},
    z_flow_prev::Dict{Int, Float64},
) = sqrt(sum((z_flow[e] - get(z_flow_prev, e, 0.0))^2 for e in keys(z_flow)))

function eps_primal_branch_dc(
    partitions::Dict{Int, PartitionDC},
    z_flow::Dict{Int, Float64},
    cut_edges::Vector{Int},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    dim_pri = sum(length(_zone_cut_edges(part, cut_edges)) for (_, part) in partitions)
    return sqrt(max(1, dim_pri)) * ϵ_abs +
           ϵ_rel * max(sqrt(sum(z_flow[e]^2 for e in cut_edges)), 1e-12)
end

function eps_dual_branch_dc(
    z_flow::Dict{Int, Float64},
    z_flow_prev::Dict{Int, Float64},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    dim_du = max(1, length(z_flow))
    return sqrt(dim_du) * ϵ_abs +
           ϵ_rel * sqrt(sum((z_flow[e] - get(z_flow_prev, e, 0.0))^2 for e in keys(z_flow)))
end

##########################
# DUAIS (FLOW)           #
##########################

function update_duals_flow!(
    duals::Dict{Int, Dict{Int, Float64}},
    partitions::Dict{Int, PartitionDC},
    z_flow::Dict{Int, Float64},
    cut_edges::Vector{Int},
    rho::Float64,
)
    for (z, part) in partitions
        m = part.model::JuMP.Model
        @inbounds for e in cut_edges
            l_local = get(part.edge_to_local, e, nothing)
            if l_local !== nothing
                p_val = _dc_flow_val_model(m, part, l_local)
                duals[z][e] = get(duals[z], e, 0.0) + rho * (p_val - z_flow[e])
            end
        end
    end
    return nothing
end



function compute_v2_consensus_vec!(v2bar_vec::Vector{Float64},
    partitions,
    shared_buses_sorted::Vector{Int}, b2pos::Dict{Int, Int})
    fill!(v2bar_vec, 0.0)
    counts = zeros(Int, length(shared_buses_sorted))
    for (_, sub) in partitions
        v2s = snapshot_v2_shared!(sub)
        lbs = sub.shared_buses_local
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            v2bar_vec[j] += v2s[i]
            counts[j] += 1
        end
    end
    @inbounds for j in eachindex(v2bar_vec)
        counts[j] != 0 && (v2bar_vec[j] /= counts[j])
    end
end

# versões baseadas em Dict (mantidas para edges)
function compute_rt_consensus!(zR::Dict{Int, Float64}, zT::Dict{Int, Float64},
    partitions, cut_edges::Vector{Int})
    sumsR = Dict{Int, Float64}()
    sumsT = Dict{Int, Float64}()
    cnts = Dict{Int, Int}()
    for (_, sub) in partitions
        m = sub.model::JuMP.Model
        @inbounds for e in cut_edges
            l = get(sub.edge_to_local, e, nothing)
            l === nothing && continue
            R, T = _rt_from_exprs_socp(m, l)
            sumsR[e] = get(sumsR, e, 0.0) + value(R)
            sumsT[e] = get(sumsT, e, 0.0) + value(T)
            cnts[e] = get(cnts, e, 0) + 1
        end
    end
    @inbounds for e in cut_edges
        n = get(cnts, e, 0)
        n > 0 && (zR[e] = sumsR[e] / n; zT[e] = sumsT[e] / n)
    end
end

function compute_pq_consensus!(zP::Dict{Int, Float64}, zQ::Dict{Int, Float64},
    partitions, cut_edges::Vector{Int})
    sumsP = Dict{Int, Float64}()
    sumsQ = Dict{Int, Float64}()
    cnts = Dict{Int, Int}()
    for (_, sub) in partitions
        m = sub.model::JuMP.Model
        @inbounds for e in cut_edges
            l = get(sub.edge_to_local, e, nothing)
            l === nothing && continue
            P, Q = _socp_flow_from_exprs_existing(m, l)
            sumsP[e] = get(sumsP, e, 0.0) + value(P)
            sumsQ[e] = get(sumsQ, e, 0.0) + value(Q)
            cnts[e] = get(cnts, e, 0) + 1
        end
    end
    @inbounds for e in cut_edges
        n = get(cnts, e, 0)
        n > 0 && (zP[e] = sumsP[e] / n; zQ[e] = sumsQ[e] / n)
    end
end

# ----------------- Duais (forma vetorial p/ BUS) -----------------
@inline function update_duals_v2_vec!(dual_vec_zone::Vector{Float64},
    sub::PartitionSOCP,
    v2bar_vec::Vector{Float64}, v2bar_prev_vec::Vector{Float64},
    rho::Float64, α::Union{Nothing, Float64},
    shared_buses_sorted::Vector{Int}, b2pos::Dict{Int, Int})
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

dual_residual_bus_socp_vec(v2bar_vec, v2bar_prev_vec) =
    sqrt(sum((v2bar_vec[j] - v2bar_prev_vec[j])^2 for j in eachindex(v2bar_vec)))

function eps_primal_bus_socp_vec(
    partitions,
    v2bar_vec,
    shared_buses_sorted,
    b2pos,
    ϵ_abs,
    ϵ_rel,
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
    return sqrt(N) * ϵ_abs + ϵ_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

function eps_dual_bus_socp_vec(v2bar_vec, v2bar_prev_vec, ϵ_abs, ϵ_rel)
    N = max(1, length(v2bar_vec))
    return sqrt(N) * ϵ_abs +
           ϵ_rel * sqrt(sum((v2bar_vec[j] - v2bar_prev_vec[j])^2
                    for j in eachindex(v2bar_vec)))
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

dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges) =
    sqrt(
        sum(
            (get(zR, e, 0.0) - get(zR_prev, e, 0.0))^2 +
            (get(zT, e, 0.0) - get(zT_prev, e, 0.0))^2
            for e in cut_edges
        ),
    )

function eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, ϵ_abs, ϵ_rel)
    dim = 0
    for (_, sub) in partitions
        dim += 2 * length(_zone_cut_edges(sub, cut_edges))
    end
    z_norm_sq = 0.0
    @inbounds for e in cut_edges
        z_norm_sq += get(zR, e, 0.0)^2 + get(zT, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * ϵ_abs + ϵ_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, ϵ_abs, ϵ_rel)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
    return sqrt(dim) * ϵ_abs + ϵ_rel * s
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

dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges) =
    sqrt(
        sum(
            (get(zP, e, 0.0) - get(zP_prev, e, 0.0))^2 +
            (get(zQ, e, 0.0) - get(zQ_prev, e, 0.0))^2
            for e in cut_edges
        ),
    )

function eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
    dim = 0
    for (_, sub) in partitions
        dim += 2 * length(_zone_cut_edges(sub, cut_edges))
    end
    z_norm_sq = 0.0
    @inbounds for e in cut_edges
        z_norm_sq += get(zP, e, 0.0)^2 + get(zQ, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * ϵ_abs + ϵ_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
    return sqrt(dim) * ϵ_abs + ϵ_rel * s
end
# ================== RUNNERS (espelhando DC/AC) ==================