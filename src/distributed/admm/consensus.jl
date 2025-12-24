
function compute_theta_consensus_vec!(
    theta_consensus_vec::Vector{Float64},
    partitions::Dict{Int, PartitionDC},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    fill!(theta_consensus_vec, 0.0)
    counts = zeros(Int, length(shared_buses_sorted))

    # média das theta nas barras compartilhadas
    for (_, part) in partitions
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


function compute_ac_bus_consensus_vec!(
    theta_consensus_vec::Vector{Float64},
    vmag_consensus_vec::Vector{Float64},
    partitions::Dict{Int, PartitionAC},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    fill!(theta_consensus_vec, 0.0)
    fill!(vmag_consensus_vec, 0.0)
    sums_sin = zeros(Float64, length(shared_buses_sorted))
    sums_cos = zeros(Float64, length(shared_buses_sorted))
    counts = zeros(Int, length(shared_buses_sorted))

    for (_, part) in partitions
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