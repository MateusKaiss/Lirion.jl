# ======================================================================
# ALADIN Workspace Components — composable sub-structs
# ======================================================================
#
# Analogous to `admm/core/components.jl` for ADMM: extracts the
# "convergence metrics" and "iteration state" concerns out of every
# ALADIN workspace god-struct. Every concrete workspace
# (DCStateWorkspaceALADIN, ACHybridWorkspaceALADIN, ...) now composes
# `residuals::ALADINResidualMetrics` and `state::ALADINIterationState`.
#
# ALADIN uses a DIFFERENT convergence criterion from ADMM (HFD 2016
# Theorem 4.1): the stopping rule is
#
#     r_primal  < eps_primal   AND   r_progress < eps_progress
#
# where
#     r_primal   = ‖∑_i A_i y_i‖₂       (consensus violation)
#     r_progress = ‖y_i - x_iᵏ‖₂        (Phase-1 progress / KKT proxy)
#
# `lambda_step = ‖λᵏ⁺¹ - λᵏ‖₂` is retained as a diagnostic only; it is
# logged to `ALADINHistory` but does not participate in the stopping
# test.
#
# Legacy flat access (`ws.r_primal`, `ws.r_progress`, `ws.iter`, …)
# continues to work through the `getproperty` / `setproperty!`
# forwarding layer below — no existing call site needs to change.
# ======================================================================

Base.@kwdef mutable struct ALADINResidualMetrics
    r_primal::Float64 = 0.0
    r_progress::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_progress::Float64 = 0.0
end

Base.@kwdef mutable struct ALADINIterationState
    hist::ALADINHistory = ALADINHistory()
    iter::Int = 0
end

# ----------------------------------------------------------------------
# LambdaState — per-consensus-variable dual storage (Phase 5b)
# ----------------------------------------------------------------------
#
# Each consensus variable (θ, |V|, v², p_flow, q_flow, r, t) gets one
# `LambdaState` holding its current and previous Lagrange multiplier
# vectors. Workspaces carry them in a typed `vars::NamedTuple` keyed on
# the canonical symbol (`:theta`, `:v_mag`, `:v2`, `:p_flow`, …).
#
# For DC-only couplings the historical `lambda_bus` / `lambda_edge`
# names are ALIASED to `vars.theta.lambda` / `vars.p_flow.lambda`
# through the forwarding layer below — no call site needs to change.

mutable struct LambdaState
    lambda::Vector{Float64}
    lambda_prev::Vector{Float64}
end

LambdaState(n::Int) = LambdaState(zeros(n), zeros(n))

# ----------------------------------------------------------------------
# getproperty / setproperty! forwarding on AbstractALADINWorkspace
# ----------------------------------------------------------------------
#
# Every concrete subtype must declare
#
#     residuals::ALADINResidualMetrics = ALADINResidualMetrics()
#     state::ALADINIterationState      = ALADINIterationState()
#
# in addition to its algorithm-specific fields. Literal-Symbol branches
# are constant-folded by the compiler when accessed with a static
# field name (the common case for `ws.field`), so the forwarding has
# zero runtime cost.

@inline function Base.getproperty(ws::AbstractALADINWorkspace, s::Symbol)
    s === :r_primal      && return getfield(ws, :residuals).r_primal
    s === :r_progress    && return getfield(ws, :residuals).r_progress
    s === :lambda_step   && return getfield(ws, :residuals).lambda_step
    s === :eps_primal    && return getfield(ws, :residuals).eps_primal
    s === :eps_progress  && return getfield(ws, :residuals).eps_progress
    s === :iter          && return getfield(ws, :state).iter
    s === :hist          && return getfield(ws, :state).hist

    # Per-variable lambda forwarding (Phase 5b). Each concrete workspace
    # that has been migrated carries `vars::NamedTuple{..., LambdaState}`
    # keyed on canonical variable symbols (:theta, :v_mag, :v2, :p_flow,
    # :q_flow, :r, :t). The compile-time `hasfield` guards elide the
    # whole block for workspaces that have not been migrated.
    T = typeof(ws)
    if hasfield(T, :vars)
        vars = getfield(ws, :vars)
        VT = typeof(vars)

        # --- theta (DC alias: lambda_bus) ---
        if hasfield(VT, :theta)
            s === :lambda_theta      && return getfield(vars, :theta).lambda
            s === :lambda_theta_prev && return getfield(vars, :theta).lambda_prev
            s === :lambda_bus        && return getfield(vars, :theta).lambda
            s === :lambda_bus_prev   && return getfield(vars, :theta).lambda_prev
        end
        # --- v_mag (ALADIN stores as lambda_vm) ---
        if hasfield(VT, :v_mag)
            s === :lambda_vm         && return getfield(vars, :v_mag).lambda
            s === :lambda_vm_prev    && return getfield(vars, :v_mag).lambda_prev
        end
        # --- v2 (SOCWR) ---
        if hasfield(VT, :v2)
            s === :lambda_v2         && return getfield(vars, :v2).lambda
            s === :lambda_v2_prev    && return getfield(vars, :v2).lambda_prev
        end
        # --- p_flow (DC alias: lambda_edge) ---
        if hasfield(VT, :p_flow)
            s === :lambda_p          && return getfield(vars, :p_flow).lambda
            s === :lambda_p_prev     && return getfield(vars, :p_flow).lambda_prev
            s === :lambda_edge       && return getfield(vars, :p_flow).lambda
            s === :lambda_edge_prev  && return getfield(vars, :p_flow).lambda_prev
        end
        # --- q_flow ---
        if hasfield(VT, :q_flow)
            s === :lambda_q          && return getfield(vars, :q_flow).lambda
            s === :lambda_q_prev     && return getfield(vars, :q_flow).lambda_prev
        end
        # --- r (SOCWR) ---
        if hasfield(VT, :r)
            s === :lambda_r          && return getfield(vars, :r).lambda
            s === :lambda_r_prev     && return getfield(vars, :r).lambda_prev
        end
        # --- t (SOCWR) ---
        if hasfield(VT, :t)
            s === :lambda_t          && return getfield(vars, :t).lambda
            s === :lambda_t_prev     && return getfield(vars, :t).lambda_prev
        end
    end

    return getfield(ws, s)
end

@inline function Base.setproperty!(ws::AbstractALADINWorkspace, s::Symbol, v)
    if s === :r_primal
        getfield(ws, :residuals).r_primal = v
        return v
    elseif s === :r_progress
        getfield(ws, :residuals).r_progress = v
        return v
    elseif s === :lambda_step
        getfield(ws, :residuals).lambda_step = v
        return v
    elseif s === :eps_primal
        getfield(ws, :residuals).eps_primal = v
        return v
    elseif s === :eps_progress
        getfield(ws, :residuals).eps_progress = v
        return v
    elseif s === :iter
        getfield(ws, :state).iter = v
        return v
    elseif s === :hist
        getfield(ws, :state).hist = v
        return v
    end
    setfield!(ws, s, v)
    return v
end
