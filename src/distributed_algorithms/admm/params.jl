Base.@kwdef struct ADMMParams
    max_iter::Int
    eps_abs::Float64
    eps_rel::Float64
    adaptive::Bool
    mu_incr::Float64
    mu_decr::Float64
    tau_incr::Float64
    tau_decr::Float64
end

function build_admm_params(config::DistributedOPFModel)
    return ADMMParams(;
        max_iter = config.max_iter,
        eps_abs = config.eps_abs,
        eps_rel = config.eps_rel,
        adaptive = config.adaptive_rho,
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
end
