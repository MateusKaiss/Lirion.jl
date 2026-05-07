# ============================================================
# Workspace hierarchy + dual state
# ============================================================

abstract type AbstractALADINWorkspace <: AbstractDistributedWorkspace end

mutable struct LambdaState
    lambda::Vector{Float64}
    lambda_prev::Vector{Float64}
    nu::Vector{Float64}
end

LambdaState(n::Int) = LambdaState(zeros(n), zeros(n), zeros(n))

# ============================================================
# Params container
# ============================================================

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

# ============================================================
# Builder
# ============================================================

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

# ============================================================
# Coordinator QP optimizer factory
# ============================================================
#
# Defaults to OSQP (primal QP, indirect method) when the user does
# not pass an explicit factory. Tolerances tuned for ALADIN's
# coordinator step — tighter than typical OSQP defaults to match the
# precision needed by Phase-2 Newton updates.

function _resolve_qp_optimizer(user_opt)
    user_opt === nothing || return user_opt
    return JuMP.optimizer_with_attributes(
        OSQP.Optimizer,
        "eps_abs" => 1e-8,
        "eps_rel" => 1e-8,
        "max_iter" => 200_000,
        "verbose" => false,
        "polish" => true,
        "polish_refine_iter" => 5,
        "linsys_solver" => "qdldl",
    )
end
