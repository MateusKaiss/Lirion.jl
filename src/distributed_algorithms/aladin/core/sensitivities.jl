# ============================================================
# ALADIN — Per-zone kernels shared across DC / AC / SOCWR
# ============================================================

# ------------------------------------------------------------
# Phase-1 → buffer
# ------------------------------------------------------------

function extract_solve_values!(zone::AbstractAladinZoneData, model::JuMP.Model)
    y = zone.y
    @inbounds for i in 1:zone.n_vars
        y[i] = JuMP.value(zone.all_vars[i])
    end
    return nothing
end

# ------------------------------------------------------------
# Cost gradient
# ------------------------------------------------------------

function compute_cost_gradient!(zone::AbstractAladinZoneData)
    g = zone.g
    fill!(g, 0.0)
    y = zone.y
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        g[col] = zone.pg_c2_sbase2[k] * y[col] + zone.pg_c1_sbase[k]
    end
    return nothing
end

# ------------------------------------------------------------
# Active bounds at current iterate
# ------------------------------------------------------------

const _BOUND_DETECT_MODE = Ref(:primal)

# Box-bound activity test shared by every model; recorded in aladin_runtime_manifest().
const _BOUND_DETECT_TOL = 1e-6
const _BOUND_DETECT_DUAL_TOL = 1e-3

function detect_active_bounds!(
    zone::AbstractAladinZoneData;
    tol::Float64 = _BOUND_DETECT_TOL,
    dual_tol::Float64 = _BOUND_DETECT_DUAL_TOL,
)
    active = zone.active_bound_cols
    empty!(active)
    y = zone.y
    use_duals = _BOUND_DETECT_MODE[] === :primal_dual &&
        (try; JuMP.has_duals(zone.local_model); catch; false; end)
    @inbounds for i in 1:zone.n_vars
        lo = zone.lb[i]
        hi = zone.ub[i]
        v = y[i]
        if !isnan(lo) && v - lo < tol
            push!(active, i)
        elseif !isnan(hi) && hi - v < tol
            push!(active, i)
        elseif use_duals
            vr = zone.all_vars[i]
            du = 0.0
            try
                JuMP.has_lower_bound(vr) &&
                    (du = max(du, abs(JuMP.dual(JuMP.LowerBoundRef(vr)))))
                JuMP.has_upper_bound(vr) &&
                    (du = max(du, abs(JuMP.dual(JuMP.UpperBoundRef(vr)))))
            catch
            end
            du > dual_tol && push!(active, i)
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Proximal curvature ρ·σ deposited on diagonal column classes
# ------------------------------------------------------------

@inline function _accumulate_prox_diagonal!(
    deposit!::F,
    rho::Float64,
    classes,
) where {F}
    rho > 0.0 || return nothing
    @inbounds for (cols, sigma) in classes
        s = rho * sigma
        s > 0.0 || continue
        for c in cols
            deposit!(c, s)
        end
    end
    return nothing
end
