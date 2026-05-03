# ============================================================
# ALADIN-SOCWR — Per-iteration sensitivities
# ============================================================
#
# Hot-path kernels refreshed every ALADIN iteration:
#   - `extract_solve_values!` snapshots Phase-1 SOCP solution into
#     `zone.y`;
#   - `compute_cost_gradient!` writes ∇f(y) into `zone.g`;
#   - `_refresh_hess_diag!` rebuilds `zone.hess_diag` = base + ρ·σ
#     contributions + δ floor (HFD Theorem 4 PSD enforcement);
#   - `refresh_lagrangian_hessian!` is the public hook used by the
#     coordinator pipeline.
#
# Struct definition and static initialisation live in `zone_data.jl`;
# active-set logic (bounds, conic normals, Jacobian) lives in
# `active_set.jl`.

# ============================================================
# Phase-1 → buffer
# ============================================================

function extract_solve_values!(zone::AladinZoneDataSOCWR, model::JuMP.Model)
    y = zone.y
    @inbounds for i in 1:zone.n_vars
        y[i] = JuMP.value(zone.all_vars[i])
    end
    return nothing
end

# ============================================================
# Cost gradient ∇f(y)
# ============================================================

function compute_cost_gradient!(zone::AladinZoneDataSOCWR)
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
# Proximal Hessian refresh
# ============================================================
#
# HFD 2016 Theorem 4.1 allows H_i to be any positive-definite
# approximation of the local Lagrangian Hessian. The cost term alone
# is PSD but only on pg columns, which leaves the remaining classes
# (θ, v², R, T, qg, curtailment) with zero curvature. We recover
# strict positive definiteness on the whole space by adding the
# curvature of the proximal term itself — exactly as the Phase-1
# subproblem optimises. This is the same construction used by
# ALADIN-DC and is the standard SQP PSD fix.

function _refresh_hess_diag!(
    zone::AladinZoneDataSOCWR,
    rho::Float64;
    sigma_theta::Float64,
    sigma_v2::Float64,
    sigma_r::Float64,
    sigma_t::Float64,
    sigma_pg::Float64,
    sigma_qg::Float64,
    sigma_other::Float64,
    delta::Float64 = 1e-6,
)
    copyto!(zone.hess_diag, zone.base_hess_diag)

    if rho > 0.0
        add_rho! = (col_set, σ) -> begin
            s = rho * σ
            s > 0.0 || return
            @inbounds for c in col_set
                zone.hess_diag[c] += s
            end
        end
        add_rho!(zone.theta_cols, sigma_theta)
        add_rho!(zone.v2_cols, sigma_v2)
        add_rho!(zone.r_cols, sigma_r)
        add_rho!(zone.t_cols, sigma_t)
        add_rho!(zone.pg_cols, sigma_pg)
        add_rho!(zone.qg_cols, sigma_qg)
        add_rho!(zone.curt_cols, sigma_other)
    end

    # Ensure a strictly positive floor on every column (a column with
    # zero σ would still leave a zero on the diagonal, making the QP
    # indefinite on that direction).
    @inbounds for c in 1:zone.n_vars
        if zone.hess_diag[c] < delta
            zone.hess_diag[c] = delta
        end
    end
    @assert all(isfinite, zone.hess_diag) "SOCWR Hessian diagonal has non-finite entries"
    return nothing
end

# ============================================================
# Sensitivity refresh hook (called from aladin/core/local_solve.jl)
# ============================================================

@inline function refresh_lagrangian_hessian!(
    zone::AladinZoneDataSOCWR,
    ws::AbstractALADINWorkspace,
)
    _refresh_hess_diag!(
        zone,
        ws.rho[];
        sigma_theta = ws.sigma_theta,
        sigma_v2 = ws.sigma_v2,
        sigma_r = ws.sigma_r,
        sigma_t = ws.sigma_t,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )
    return nothing
end
