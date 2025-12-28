# ============================================================
# Consensus builders (ADMM)
# ============================================================

# ------------------------------------------------------------
# DC consensus: bus variables (theta)
# ------------------------------------------------------------

function compute_theta_consensus_vec!(
    theta_consensus_vec::Vector{Float64},
    partitions::Partitions{PartitionDC},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    fill!(theta_consensus_vec, 0.0)
    counts = zeros(Int, length(shared_buses_sorted))

    # mean theta over shared buses
    for part in partitions
        thetas = snapshot_theta_shared!(part)
        lbs = part.shared_buses_local
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            theta_consensus_vec[j] += thetas[i]
            counts[j] += 1
        end
    end

    @inbounds for j in eachindex(theta_consensus_vec)
        counts[j] != 0 && (theta_consensus_vec[j] /= counts[j])
    end
    return nothing
end

# ------------------------------------------------------------
# DC consensus: branch variables (flow on cut edges)
# ------------------------------------------------------------

function compute_flow_consensus_consensus!(
    flow_consensus::Dict{Int, Float64},
    partitions::Partitions{PartitionDC},
    cut_edges::Vector{Int},
)
    sums = Dict{Int, Float64}()
    cnts = Dict{Int, Int}()

    # mean p_e over cut edges
    for part in partitions
        m = part.model::JuMP.Model
        @inbounds for e in part.zone_edges
            l_local = get(part.edge_to_local, e, 0)
            l_local == 0 && continue
            p = _dc_flow_val_model(m, part, l_local)
            sums[e] = get(sums, e, 0.0) + p
            cnts[e] = get(cnts, e, 0) + 1
        end
    end

    # escreve apenas para os cut_edges globais (mantém contrato)
    @inbounds for e in cut_edges
        n = get(cnts, e, 0)
        n > 0 && (flow_consensus[e] = sums[e] / n)
    end
    return nothing
end

# ------------------------------------------------------------
# AC consensus: bus variables (theta, |V|)
# ------------------------------------------------------------

function compute_ac_bus_consensus_vec!(
    theta_consensus_vec::Vector{Float64},
    vmag_consensus_vec::Vector{Float64},
    partitions::Partitions{PartitionAC},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    fill!(theta_consensus_vec, 0.0)
    fill!(vmag_consensus_vec, 0.0)

    sums_sin = zeros(Float64, length(shared_buses_sorted))
    sums_cos = zeros(Float64, length(shared_buses_sorted))
    counts = zeros(Int, length(shared_buses_sorted))

    for part in partitions
        thetas = snapshot_theta_ac_shared!(part)
        Vs = snapshot_vmag_shared!(part)
        lbs = part.shared_buses_local

        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            s, c = sincos(thetas[i])
            sums_sin[j] += s
            sums_cos[j] += c
            vmag_consensus_vec[j] += Vs[i]
            counts[j] += 1
        end
    end

    @inbounds for j in eachindex(theta_consensus_vec)
        if counts[j] != 0
            theta_consensus_vec[j] = atan(sums_sin[j], sums_cos[j])
            vmag_consensus_vec[j] /= counts[j]
        end
    end
    return nothing
end

# ------------------------------------------------------------
# AC consensus: branch variables (P, Q on cut edges)
# ------------------------------------------------------------

function compute_ac_branch_consensus!(
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    partitions::Partitions{PartitionAC},
    cut_edges::Vector{Int},
)
    sumsP = Dict{Int, Float64}()
    sumsQ = Dict{Int, Float64}()
    cntsP = Dict{Int, Int}()
    cntsQ = Dict{Int, Int}()

    for part in partitions
        # itera só nos edges que existem na zona
        @inbounds for e in part.zone_edges
            l_local = get(part.edge_to_local, e, 0)
            l_local == 0 && continue

            P_loc, Q_loc = _ac_flow_from_exprs_existing(part, l_local)
            pval = JuMP.value(P_loc)
            qval = JuMP.value(Q_loc)

            if isfinite(pval)
                sumsP[e] = get(sumsP, e, 0.0) + pval
                cntsP[e] = get(cntsP, e, 0) + 1
            end
            if isfinite(qval)
                sumsQ[e] = get(sumsQ, e, 0.0) + qval
                cntsQ[e] = get(cntsQ, e, 0) + 1
            end
        end
    end

    @inbounds for e in cut_edges
        nP = get(cntsP, e, 0)
        nQ = get(cntsQ, e, 0)
        nP > 0 && (zP[e] = sumsP[e] / nP)
        nQ > 0 && (zQ[e] = sumsQ[e] / nQ)
    end

    return nothing
end

# ------------------------------------------------------------
# SOCP consensus: bus variables (V^2)
# ------------------------------------------------------------

function compute_v2_consensus_vec!(
    v2bar_vec::Vector{Float64},
    partitions::Partitions{PartitionSOCP},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    fill!(v2bar_vec, 0.0)
    counts = zeros(Int, length(shared_buses_sorted))

    for sub in partitions
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
    return nothing
end

# ------------------------------------------------------------
# SOCP consensus: branch variables (Dict-based, kept for edges)
# ------------------------------------------------------------

function compute_rt_consensus!(
    zR::Dict{Int, Float64},
    zT::Dict{Int, Float64},
    partitions::Partitions{PartitionSOCP},
    cut_edges::Vector{Int},
)
    sumsR = Dict{Int, Float64}()
    sumsT = Dict{Int, Float64}()
    cnts = Dict{Int, Int}()

    for sub in partitions
        m = sub.model::JuMP.Model
        @inbounds for e in sub.zone_edges
            l = get(sub.edge_to_local, e, 0)
            l == 0 && continue
            R, T = _rt_from_exprs_socp(m, l)
            sumsR[e] = get(sumsR, e, 0.0) + JuMP.value(R)
            sumsT[e] = get(sumsT, e, 0.0) + JuMP.value(T)
            cnts[e] = get(cnts, e, 0) + 1
        end
    end

    @inbounds for e in cut_edges
        n = get(cnts, e, 0)
        if n > 0
            zR[e] = sumsR[e] / n
            zT[e] = sumsT[e] / n
        end
    end
    return nothing
end

function compute_pq_consensus!(
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    partitions::Partitions{PartitionSOCP},
    cut_edges::Vector{Int},
)
    sumsP = Dict{Int, Float64}()
    sumsQ = Dict{Int, Float64}()
    cnts = Dict{Int, Int}()

    for sub in partitions
        m = sub.model::JuMP.Model
        @inbounds for e in sub.zone_edges
            l = get(sub.edge_to_local, e, 0)
            l == 0 && continue
            P, Q = _socp_flow_from_exprs_existing(m, l)
            sumsP[e] = get(sumsP, e, 0.0) + JuMP.value(P)
            sumsQ[e] = get(sumsQ, e, 0.0) + JuMP.value(Q)
            cnts[e] = get(cnts, e, 0) + 1
        end
    end

    @inbounds for e in cut_edges
        n = get(cnts, e, 0)
        if n > 0
            zP[e] = sumsP[e] / n
            zQ[e] = sumsQ[e] / n
        end
    end
    return nothing
end
