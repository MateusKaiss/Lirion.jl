# ======================================================================
# ADMM Workspace Components — composable sub-structs
# ======================================================================
#
# These small mutable structs factor out the orthogonal "residual
# metrics" and "iteration state" concerns that used to live as six
# flat fields on every ADMM workspace god-struct. Each concrete
# workspace (DCStateWorkspace, ACHybridWorkspace, ...) now composes
# `residuals::ResidualMetrics` and `state::IterationState`.
#
# Legacy flat access (`ws.r_primal`, `ws.iter`, etc.) continues to
# work thanks to the `getproperty` / `setproperty!` forwarding layer
# declared at the end of this file — so the 100+ existing call sites
# need no changes. New code is free to use `ws.residuals.r_primal`
# directly, which is more idiomatic and faster under type inference.
#
# Per-component residuals (r_primal_theta, r_primal_bus, etc.) stay
# on the concrete workspaces because they are coupling-specific and
# only consumed by the matching `rho_update!` dispatch. ResidualMetrics
# is reserved for the *aggregate* numbers every workspace must expose.
# ======================================================================

Base.@kwdef mutable struct ResidualMetrics
    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0
end

Base.@kwdef mutable struct IterationState
    hist::ADMMHistory = ADMMHistory()
    iter::Int = 0
end

# ----------------------------------------------------------------------
# Legacy-field forwarding for AbstractADMMWorkspace
# ----------------------------------------------------------------------
#
# Every concrete subtype must declare
#
#     residuals::ResidualMetrics = ResidualMetrics()
#     state::IterationState      = IterationState()
#
# in addition to its coupling-specific fields. The forwarding layer
# below lets existing code continue to use the flat names
# (`ws.r_primal`, `ws.iter`, …) while the new composed storage is
# used internally.
#
# The branches use literal-Symbol equality so the Julia compiler
# constant-folds the whole dispatch when the symbol is known at
# compile time — which is the case for every `ws.field` expression.
# Result: zero runtime cost compared to a direct field access.

@inline function Base.getproperty(ws::AbstractADMMWorkspace, s::Symbol)
    if s === :r_primal
        return getfield(ws, :residuals).r_primal
    elseif s === :r_dual
        return getfield(ws, :residuals).r_dual
    elseif s === :eps_primal
        return getfield(ws, :residuals).eps_primal
    elseif s === :eps_dual
        return getfield(ws, :residuals).eps_dual
    elseif s === :iter
        return getfield(ws, :state).iter
    elseif s === :hist
        return getfield(ws, :state).hist
    end
    return getfield(ws, s)
end

@inline function Base.setproperty!(ws::AbstractADMMWorkspace, s::Symbol, v)
    if s === :r_primal
        getfield(ws, :residuals).r_primal = v
        return v
    elseif s === :r_dual
        getfield(ws, :residuals).r_dual = v
        return v
    elseif s === :eps_primal
        getfield(ws, :residuals).eps_primal = v
        return v
    elseif s === :eps_dual
        getfield(ws, :residuals).eps_dual = v
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
