# =======================================================================
# ALADIN — Parameter container
# =======================================================================
struct ALADINParams
    max_iter::Int
    rho::Float64              # local proximal penalty
    mu::Float64               # coordination-QP slack penalty
    eps_abs::Float64
    eps_rel::Float64
    active_set_tol::Float64
end

@inline function build_aladin_params(cfg::ALADINConfig)::ALADINParams
    return ALADINParams(
        cfg.max_iter,
        cfg.rho,
        cfg.mu,
        cfg.eps_abs,
        cfg.eps_rel,
        cfg.active_set_tol,
    )
end
