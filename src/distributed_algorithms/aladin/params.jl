# =======================================================================
# ALADIN — Parameter container (canonical HFD 2016)
# =======================================================================

struct ALADINParams
    max_iter::Int

    # Penalty parameters
    rho::Float64              # local proximal penalty
    mu::Float64               # consensus slack penalty in coordination QP

    # HFD Σ_i diagonal weights (per variable class)
    sigma_theta_shared::Float64
    sigma_theta_other::Float64
    sigma_pg::Float64
    sigma_other::Float64

    # Line-search parameters
    alpha1::Float64
    alpha2::Float64
    alpha3::Float64

    # Convergence (HFD 2016 Theorem 4.1)
    eps_abs::Float64
    eps_rel::Float64
    active_set_tol::Float64

    adaptive::Bool
    rho_growth::Float64
    mu_growth::Float64
    rho_upper::Float64
    mu_upper::Float64
end

@inline function build_aladin_params(cfg::ALADINConfig)::ALADINParams
    return ALADINParams(
        cfg.max_iter,
        cfg.rho,
        cfg.mu,
        cfg.sigma_theta_shared,
        cfg.sigma_theta_other,
        cfg.sigma_pg,
        cfg.sigma_other,
        cfg.alpha1,
        cfg.alpha2,
        cfg.alpha3,
        cfg.eps_abs,
        cfg.eps_rel,
        cfg.active_set_tol,
        cfg.adaptive,
        cfg.rho_growth,
        cfg.mu_growth,
        cfg.rho_upper,
        cfg.mu_upper,
    )
end
