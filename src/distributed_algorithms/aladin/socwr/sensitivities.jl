# ============================================================
# ALADIN-SOCWR — Per-iteration sensitivities
# ============================================================

# ------------------------------------------------------------
# Proximal Hessian refresh
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Sensitivity refresh hook
# ------------------------------------------------------------

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
