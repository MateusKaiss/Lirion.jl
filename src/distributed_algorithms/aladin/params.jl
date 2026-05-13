# ============================================================
# Workspace hierarchy + dual state
# ============================================================

abstract type AbstractALADINWorkspace <: AbstractDistributedWorkspace end
abstract type AbstractAladinZoneData end

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
    stopping_criterion::Symbol
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
        cfg.stopping_criterion,
    )
end

# ============================================================
# Coordinator QP optimizer factory
# ============================================================

const _DEFAULT_QP_FACTORY = JuMP.optimizer_with_attributes(
    OSQP.Optimizer,
    "eps_abs" => 1e-8,
    "eps_rel" => 1e-8,
    "max_iter" => 200_000,
    "verbose" => false,
    "polish" => true,
    "polish_refine_iter" => 5,
    "linsys_solver" => "qdldl",
)

function _resolve_qp_optimizer(user_opt)
    user_opt === nothing && return _DEFAULT_QP_FACTORY
    return _silent_factory(user_opt)
end
