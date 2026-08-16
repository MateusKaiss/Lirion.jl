# ============================================================
# SplitBus N-way injection consensus (Σ p_z = 0 over all copies)
# ============================================================


function _consensus_inj_groups!(
    target::AbstractVector{Float64},
    partitions::Partitions,
    inj_groups::Vector{Tuple{Int, Vector{Int}}},
    inj_copy_zone::Vector{Int},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)
    isempty(inj_groups) && return nothing
    fill!(target, 0.0)
    @inbounds for (b_id, cidxs) in inj_groups
        s = 0.0
        for c in cidxs
            z = inj_copy_zone[c]
            s += _p_inj_value(partitions.parts[z], p_inj_vars[z], b_id)
        end
        m = s / length(cidxs)
        for c in cidxs
            z = inj_copy_zone[c]
            v = _p_inj_value(partitions.parts[z], p_inj_vars[z], b_id)
            target[c] = v - m
        end
    end
    return nothing
end

function _update_inj_groups_duals!(
    duals::AbstractMatrix{Float64},
    partitions::Partitions,
    target::AbstractVector{Float64},
    rho::Float64,
    inj_groups::Vector{Tuple{Int, Vector{Int}}},
    inj_copy_zone::Vector{Int},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)
    isempty(inj_groups) && return nothing
    @inbounds for (b_id, cidxs) in inj_groups
        for c in cidxs
            z = inj_copy_zone[c]
            v = _p_inj_value(partitions.parts[z], p_inj_vars[z], b_id)
            duals[c, z] += rho * (v - target[c])
        end
    end
    return nothing
end

function _primal_residual_inj_groups(
    partitions::Partitions,
    target::AbstractVector{Float64},
    inj_groups::Vector{Tuple{Int, Vector{Int}}},
    inj_copy_zone::Vector{Int},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)
    isempty(inj_groups) && return 0.0
    acc = 0.0
    @inbounds for (b_id, cidxs) in inj_groups
        for c in cidxs
            z = inj_copy_zone[c]
            v = _p_inj_value(partitions.parts[z], p_inj_vars[z], b_id)
            acc += (v - target[c])^2
        end
    end
    return sqrt(acc)
end

function _dual_residual_inj_groups(
    target::AbstractVector{Float64},
    target_prev::AbstractVector{Float64},
    inj_groups::Vector{Tuple{Int, Vector{Int}}},
    rho::Float64,
)
    isempty(inj_groups) && return 0.0
    s = 0.0
    @inbounds for (_, cidxs) in inj_groups
        for c in cidxs
            d = target[c] - target_prev[c]
            s += d^2
        end
    end
    return rho * sqrt(s)
end

function _eps_primal_inj_groups(
    partitions::Partitions,
    target::AbstractVector{Float64},
    inj_groups::Vector{Tuple{Int, Vector{Int}}},
    inj_copy_zone::Vector{Int},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
    eps_abs::Float64,
    eps_rel::Float64,
)
    isempty(inj_groups) && return 0.0
    dim = 0
    norm_x_sq = 0.0
    norm_z_sq = 0.0
    @inbounds for (b_id, cidxs) in inj_groups
        for c in cidxs
            z = inj_copy_zone[c]
            v = _p_inj_value(partitions.parts[z], p_inj_vars[z], b_id)
            dim += 1
            norm_x_sq += v^2
            norm_z_sq += target[c]^2
        end
    end
    return sqrt(max(1, dim)) * eps_abs +
           eps_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq), 1e-12)
end

function _eps_dual_inj_groups(
    duals::AbstractMatrix{Float64},
    inj_groups::Vector{Tuple{Int, Vector{Int}}},
    inj_copy_zone::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    isempty(inj_groups) && return 0.0
    n = 0
    norm_y_sq = 0.0
    @inbounds for (_, cidxs) in inj_groups
        for c in cidxs
            z = inj_copy_zone[c]
            n += 1
            norm_y_sq += duals[c, z]^2
        end
    end
    return sqrt(max(1, n)) * eps_abs + eps_rel * sqrt(norm_y_sq)
end

function build_linear_inj_term_groups(
    pos::Int,
    inj_groups::Vector{Tuple{Int, Vector{Int}}},
    inj_copy_zone::Vector{Int},
    inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
    target::AbstractVector{Float64},
    duals_col::AbstractVector{Float64},
    rho::Float64,
)::JuMP.AffExpr
    lin = zero(JuMP.AffExpr)
    isempty(inj_groups) && return lin
    store = inj_vars[pos]
    isempty(store) && return lin
    @inbounds for (b_id, cidxs) in inj_groups
        for c in cidxs
            inj_copy_zone[c] == pos || continue
            v = get(store, b_id, nothing)
            v === nothing && continue
            lin += (duals_col[c] - rho * target[c]) * v
        end
    end
    return lin
end
