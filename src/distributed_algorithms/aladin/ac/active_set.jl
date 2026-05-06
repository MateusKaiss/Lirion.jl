# ============================================================
# ALADIN-AC — Active-set detection
# ============================================================
#
# AC equalities and active inequalities are non-linear, so the
# Jacobian is rebuilt every iteration inside `refresh_sensitivities_ac!`
# (see `sensitivities.jl`) using `MOI.Nonlinear.Evaluator`. The cached
# matrix lives in `zone.active_jac`; here we only expose the bound
# detector and the trivial accessor.
#
# `_delete_dependent_rows` is shared with DC (defined in
# `aladin/dc/active_set.jl` and reused here).

# ============================================================
# Active bounds at current iterate
# ============================================================

function detect_active_bounds!(zone::AladinZoneDataAC; tol::Float64 = 1e-6)
    active = zone.active_bound_cols
    empty!(active)
    y = zone.y
    @inbounds for i in 1:zone.n_vars
        lo = zone.lb[i]
        hi = zone.ub[i]
        v = y[i]
        if !isnan(lo) && v - lo < tol
            push!(active, i)
        elseif !isnan(hi) && hi - v < tol
            push!(active, i)
        end
    end
    return nothing
end

# ============================================================
# Active-constraint Jacobian (cached by refresh_sensitivities_ac!)
# ============================================================

@inline build_active_jacobian(
    zone::AladinZoneDataAC;
    kwargs...,
)::SparseArrays.SparseMatrixCSC{Float64, Int} = zone.active_jac
