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

function detect_active_bounds!(
    zone::AbstractAladinZoneData;
    tol::Float64 = 1e-6,
)
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
