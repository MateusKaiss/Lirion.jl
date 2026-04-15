# =======================================================================
# ALADIN — Parameter container
# =======================================================================
struct ALADINParams
    max_iter::Int

    eps_abs::Float64
    eps_rel::Float64

    adaptive::Bool
    rho_growth::Float64
    rho_upper::Float64
    mu_growth::Float64
    mu_upper::Float64
end

@inline function build_params(cfg::ALADINConfig)::ALADINParams
    return ALADINParams(
        cfg.max_iter,
        cfg.eps_abs,
        cfg.eps_rel,
        cfg.adaptive,
        cfg.rho_growth,
        cfg.rho_upper,
        cfg.mu_growth,
        cfg.mu_upper,
    )
end
