# Guard for models that have not been migrated to N-way SplitBus yet: error if
# any split bus is shared by more than two zones.
@inline function _assert_bilateral_split(sys, who::AbstractString)
    @inbounds for (bid, zones) in sys.bus_groups
        length(zones) > 2 && error(
            "$(who): N-way SplitBus (bus $(bid) shared by $(length(zones)) " *
            "zones) is not supported for this model yet. Use a partition where " *
            "each boundary bus borders a single neighbor zone.",
        )
    end
    return nothing
end

# ============================================================
# SplitBus N-way injection consensus (Σ p_z = 0 over all copies)
# ============================================================
#
# Generalizes the pairwise `_*_bus_pairs_signed` helpers to a boundary bus
# shared by any number of zones. Each split bus is a "group" with one copy per
# zone; the physical constraint is that the net injection summed over all copies
# is zero. The exchange-ADMM consensus target for a copy is its local value
# minus the group mean:  target_c = v_c - mean(v over the group).  Since the
# targets of a group sum to zero, this enforces Σ v_z → 0. For a 2-zone group
# this reduces exactly to the pairwise `±0.5(v_a - v_b)` formulation.
#
# Storage is flattened per copy: `inj_groups[g] = (bus_id, [copy_idx...])` and
# `inj_copy_zone[c]` gives the zone position of copy `c`. The per-copy target is
# held in the p_inj VarState `consensus` vector (so save_prev! tracks its prev),
# and the per-copy dual lives at `duals[c, zone_of_copy(c)]`.

function _consensus_inj_groups!(
    target::AbstractVector{Float64},
    partitions::Partitions,
    inj_groups::Vector{Tuple{Int, Vector{Int}}},
    inj_copy_zone::Vector{Int},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)
    isempty(inj_groups) && return nothing
    fill!(target, 0.0)
    @inbounds for (bid, cidxs) in inj_groups
        s = 0.0
        for c in cidxs
            z = inj_copy_zone[c]
            s += _p_inj_value(partitions.parts[z], p_inj_vars[z], bid)
        end
        m = s / length(cidxs)
        for c in cidxs
            z = inj_copy_zone[c]
            v = _p_inj_value(partitions.parts[z], p_inj_vars[z], bid)
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
    @inbounds for (bid, cidxs) in inj_groups
        for c in cidxs
            z = inj_copy_zone[c]
            v = _p_inj_value(partitions.parts[z], p_inj_vars[z], bid)
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
    @inbounds for (bid, cidxs) in inj_groups
        for c in cidxs
            z = inj_copy_zone[c]
            v = _p_inj_value(partitions.parts[z], p_inj_vars[z], bid)
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
    @inbounds for (bid, cidxs) in inj_groups
        for c in cidxs
            z = inj_copy_zone[c]
            v = _p_inj_value(partitions.parts[z], p_inj_vars[z], bid)
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

# Linear proximal term for the injection copies that live in `part`'s zone.
# `target` is the per-copy consensus vector; `duals_col` is the dual column for
# this partition (duals[:, part.position]).
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
    @inbounds for (bid, cidxs) in inj_groups
        for c in cidxs
            inj_copy_zone[c] == pos || continue
            v = get(store, bid, nothing)
            v === nothing && continue
            lin += (duals_col[c] - rho * target[c]) * v
        end
    end
    return lin
end
