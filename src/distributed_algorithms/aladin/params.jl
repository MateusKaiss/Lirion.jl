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
    check_primal_step::Bool
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
        cfg.check_primal_step,
    )
end

# ============================================================
# Coordinator QP optimizer factory
# ============================================================

# Coordinator QP default = CLARABEL (interior-point). The coordinator QP is convex
# but becomes ILL-CONDITIONED at scale (many coupling variables spanning wide
# magnitude ranges); OSQP (first-order) cannot reach tolerance there and reports
# the QP as unsolved ("coordinator QP failed"), while Clarabel (second-order,
# interior-point) converges in a few iterations at every scale — its iteration
# count is nearly insensitive to the condition number. On well-conditioned small
# QPs both solvers reach the SAME optimum (verified: identical ALADIN steps), so
# Clarabel is a strict robustness upgrade with no accuracy loss. Was OSQP.
const _DEFAULT_QP_FACTORY = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "verbose" => false,
    "tol_feas" => 1e-7,
    "tol_gap_abs" => 1e-7,
    "tol_gap_rel" => 1e-7,
    "max_iter" => 2000,
)

function _resolve_qp_optimizer(user_opt)
    user_opt === nothing && return _DEFAULT_QP_FACTORY
    return _silent_factory(user_opt)
end
