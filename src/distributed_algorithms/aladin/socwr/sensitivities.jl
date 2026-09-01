# ============================================================
# ALADIN-SOCWR — Per-iteration sensitivities
# ============================================================

# ------------------------------------------------------------
# Proximal Hessian refresh
# ------------------------------------------------------------

# Scale of the ρσ proximal deposit on the coordinator-QP diagonal (1.0 = full,
# 0.0 = none); recorded in `aladin_runtime_manifest()`.
const _SOCWR_QP_PROX_SCALE = Ref{Float64}(1.0)

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

    hd = zone.hess_diag
    # hess_diag feeds only the coordinator QP; scale the QP-side ρσ deposit
    rho_qp = rho * _SOCWR_QP_PROX_SCALE[]
    _accumulate_prox_diagonal!(
        rho_qp,
        (
            (zone.theta_cols, sigma_theta),
            (zone.v2_cols, sigma_v2),
            (zone.r_cols, sigma_r),
            (zone.t_cols, sigma_t),
            (zone.pg_cols, sigma_pg),
            (zone.qg_cols, sigma_qg),
            (zone.curt_cols, sigma_other),
        ),
    ) do c, s
        return hd[c] += s
    end

    @inbounds for c in 1:zone.n_vars
        if zone.hess_diag[c] < delta
            zone.hess_diag[c] = delta
        end
    end
    @assert all(isfinite, zone.hess_diag) "SOCWR Hessian diagonal has non-finite entries"
    return nothing
end

# ============================================================
# SOCWRNL — exact Lagrangian Hessian (W-space, AC-style)
# ============================================================

# Building blocks shared with aladin/ac/sensitivities.jl; only the proximal
# curvature classes differ (W-space θ/v²/R/T/pg/qg/curt).
function _add_prox_curvature_wspace!(
    H::SparseArrays.SparseMatrixCSC{Float64, Int},
    zone::AladinZoneDataSOCWR,
    rho::Float64;
    sigma_theta::Float64,
    sigma_v2::Float64,
    sigma_r::Float64,
    sigma_t::Float64,
    sigma_pg::Float64,
    sigma_qg::Float64,
    sigma_other::Float64,
)
    _accumulate_prox_diagonal!(
        rho,
        (
            (zone.theta_cols, sigma_theta),
            (zone.v2_cols, sigma_v2),
            (zone.r_cols, sigma_r),
            (zone.t_cols, sigma_t),
            (zone.pg_cols, sigma_pg),
            (zone.qg_cols, sigma_qg),
            (zone.curt_cols, sigma_other),
        ),
    ) do c, s
        return H[c, c] += s
    end
    return nothing
end

# Prox-in-B toggle for the SOCWRNL exact-Hessian path (cf. `_AC_QP_NO_PROX`):
# true skips the ρΣ deposit so B = ∇²L + δI.
const _SOCWRNL_QP_NO_PROX = Ref(false)

function _refresh_exact_hessian_socwrnl!(
    zone::AladinZoneDataSOCWR,
    ws::AbstractALADINWorkspace,
)
    y = zone.y
    n_vars = zone.n_vars

    evaluator, rows, strong = _build_ac_evaluator(zone, _ACTIVE_DUAL_TOL[])

    # Active NLP Jacobian C̃ (strong rows); HFD step-3 gradient correction when enabled.
    zone.active_jac, g_corr =
        _extract_active_jacobian_ac(evaluator, rows, strong, y, n_vars, _ACTIVE_DEP_TOL[])
    if _MODIFIED_GRADIENT[]
        g = zone.g
        @inbounds for c in 1:n_vars
            g[c] += g_corr[c]
        end
    end

    H = _extract_lagrangian_hessian_ac(evaluator, rows, y, n_vars)
    if !_SOCWRNL_QP_NO_PROX[]
        _add_prox_curvature_wspace!(
            H, zone, ws.rho[];
            sigma_theta = ws.sigma_theta,
            sigma_v2 = ws.sigma_v2,
            sigma_r = ws.sigma_r,
            sigma_t = ws.sigma_t,
            sigma_pg = ws.sigma_pg,
            sigma_qg = ws.sigma_qg,
            sigma_other = ws.sigma_other,
        )
    end
    n_vars > 0 && _add_tikhonov_ac!(H, n_vars)
    zone.hess = H
    zone.hess = _ensure_psd_ac!(zone.hess)   # :flip/:clip return a new matrix
    return nothing
end

# ------------------------------------------------------------
# Sensitivity refresh hook
# ------------------------------------------------------------

@inline function refresh_lagrangian_hessian!(
    zone::AladinZoneDataSOCWR,
    ws::AbstractALADINWorkspace,
)
    if zone.exact_hessian
        _refresh_exact_hessian_socwrnl!(zone, ws)
        return nothing
    end
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
