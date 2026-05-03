# ============================================================
# ALADIN-DC — Per-iteration sensitivities
# ============================================================
#
# Hot-path kernels refreshed every ALADIN iteration:
#   - `extract_solve_values!` snapshots the Phase-1 NLP solution into
#     `zone.y` so subsequent kernels read from a contiguous buffer;
#   - `compute_cost_gradient!` writes the cost gradient ∇f(y) into
#     `zone.g` (DC: only pg columns are non-zero);
#   - `refresh_lagrangian_hessian!` is a no-op for DC because the
#     Hessian diagonal is constant and was set in `zone_data.jl` at
#     build time (`hess_diag = 1e-6 + 2c₂·sbase²` on pg columns).
#
# Active-set logic lives in `active_set.jl`; struct definition and
# static initialisation live in `zone_data.jl`.

# ============================================================
# Phase-1 → buffer
# ============================================================

function extract_solve_values!(zone::AladinZoneData, model::JuMP.Model)
    y = zone.y
    @inbounds for i in 1:zone.n_vars
        y[i] = JuMP.value(zone.all_vars[i])
    end
    return nothing
end

# ============================================================
# Cost gradient ∇f(y)
# ============================================================

function compute_cost_gradient!(zone::AladinZoneData)
    g = zone.g
    fill!(g, 0.0)
    y = zone.y
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        g[col] = zone.pg_c2_sbase2[k] * y[col] + zone.pg_c1_sbase[k]
    end
    return nothing
end

# ============================================================
# Hessian refresh
# ============================================================
#
# DC has a constant diagonal Hessian (cost is quadratic in pg only,
# and the LM floor is fixed). The per-iteration refresh is therefore
# a no-op; the workspace's `hess_diag` was populated once in
# `build_aladin_zone_data` (zone_data.jl).

@inline refresh_lagrangian_hessian!(::AladinZoneData, ::AbstractALADINWorkspace) =
    nothing
