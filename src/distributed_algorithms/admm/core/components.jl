# ======================================================================
# ADMM Workspace Components — composable sub-structs
# ======================================================================
#
# Three small mutable structs factor the ADMM workspace "god-structs"
# into orthogonal concerns:
#
#   ResidualMetrics  — aggregate primal/dual residuals + tolerances
#                      (the 4 scalars checked by `stop_if_converged`).
#
#   IterationState   — outer-loop counter and history recorder.
#
#   VarState         — per-consensus-variable storage: consensus
#                      vector, its previous snapshot, the dual matrix,
#                      and the current/previous penalty. One VarState
#                      instance per symbol in `consensus_vars(F, C)`.
#                      Workspaces hold them in a typed `NamedTuple`.
#
# Legacy flat access (`ws.r_primal`, `ws.theta_consensus`,
# `ws.duals_p_flow`, `ws.rho_v`, …) continues to work thanks to the
# `getproperty` / `setproperty!` forwarding layer at the bottom of
# this file. The branches use literal-Symbol equality so the Julia
# compiler constant-folds the whole dispatch when the symbol is known
# at compile time (the case for every `ws.field` expression). For
# workspaces that do not have a particular variable (e.g. a power
# coupling has no theta), a compile-time `hasfield` guard elides the
# branch completely — zero runtime cost and no invalid field access.
#
# Per-component residual caches (`r_primal_theta`, `r_primal_bus`,
# etc.) stay on the concrete workspaces because they are
# coupling-specific and only consumed by the matching `rho_update!`
# dispatch.
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

mutable struct VarState
    consensus::Vector{Float64}
    consensus_prev::Vector{Float64}
    duals::Matrix{Float64}
    rho::Base.RefValue{Float64}
    prev_rho::Base.RefValue{Float64}
end

"""
    VarState(n::Int, n_zones::Int, rho_init::Float64; init_value=0.0)

Construct a `VarState` for a consensus variable of length `n` spanning
`n_zones` partitions, with initial penalty `rho_init` and consensus
values filled to `init_value` (e.g. 1.0 for voltage magnitudes, 0.0
for angles and flows).
"""
function VarState(
    n::Int,
    n_zones::Int,
    rho_init::Float64;
    init_value::Float64 = 0.0,
)
    consensus = fill(init_value, n)
    return VarState(
        consensus,
        fill(init_value, n),
        zeros(n, n_zones),
        Ref(rho_init),
        Ref(rho_init),
    )
end

# ----------------------------------------------------------------------
# getproperty / setproperty! forwarding on AbstractADMMWorkspace
# ----------------------------------------------------------------------

@inline function Base.getproperty(ws::AbstractADMMWorkspace, s::Symbol)
    # Aggregate residuals / iteration state (always present)
    s === :r_primal   && return getfield(ws, :residuals).r_primal
    s === :r_dual     && return getfield(ws, :residuals).r_dual
    s === :eps_primal && return getfield(ws, :residuals).eps_primal
    s === :eps_dual   && return getfield(ws, :residuals).eps_dual
    s === :iter       && return getfield(ws, :state).iter
    s === :hist       && return getfield(ws, :state).hist

    # Per-variable forwarding is only meaningful for workspaces that
    # carry a `vars::NamedTuple` field. The `hasfield` check is
    # compile-time for each concrete workspace type, so the entire
    # block is elided for workspaces that still use flat fields
    # (backward-compatible during incremental Phase 4b migration).
    T = typeof(ws)
    if hasfield(T, :vars)
        vars = getfield(ws, :vars)
        VT = typeof(vars)

        # --- theta ---
        if hasfield(VT, :theta)
            s === :theta_consensus      && return getfield(vars, :theta).consensus
            s === :theta_consensus_prev && return getfield(vars, :theta).consensus_prev
            s === :duals_theta          && return getfield(vars, :theta).duals
            s === :rho_theta            && return getfield(vars, :theta).rho
            s === :prev_rho_theta       && return getfield(vars, :theta).prev_rho
        end
        # --- v_mag ---
        if hasfield(VT, :v_mag)
            s === :v_mag_consensus      && return getfield(vars, :v_mag).consensus
            s === :v_mag_consensus_prev && return getfield(vars, :v_mag).consensus_prev
            s === :duals_v_mag          && return getfield(vars, :v_mag).duals
            s === :rho_v                && return getfield(vars, :v_mag).rho
            s === :prev_rho_v           && return getfield(vars, :v_mag).prev_rho
        end
        # --- v2 (SOCWR) ---
        if hasfield(VT, :v2)
            s === :v2_consensus         && return getfield(vars, :v2).consensus
            s === :v2_consensus_prev    && return getfield(vars, :v2).consensus_prev
            s === :duals_v2             && return getfield(vars, :v2).duals
            s === :rho_v2               && return getfield(vars, :v2).rho
            s === :prev_rho_v2          && return getfield(vars, :v2).prev_rho
        end
        # --- p_flow ---
        if hasfield(VT, :p_flow)
            s === :p_flow_consensus      && return getfield(vars, :p_flow).consensus
            s === :p_flow_consensus_prev && return getfield(vars, :p_flow).consensus_prev
            s === :duals_p_flow          && return getfield(vars, :p_flow).duals
            s === :rho_p                 && return getfield(vars, :p_flow).rho
            s === :prev_rho_p            && return getfield(vars, :p_flow).prev_rho
        end
        # --- q_flow ---
        if hasfield(VT, :q_flow)
            s === :q_flow_consensus      && return getfield(vars, :q_flow).consensus
            s === :q_flow_consensus_prev && return getfield(vars, :q_flow).consensus_prev
            s === :duals_q_flow          && return getfield(vars, :q_flow).duals
            s === :rho_q                 && return getfield(vars, :q_flow).rho
            s === :prev_rho_q            && return getfield(vars, :q_flow).prev_rho
        end
        # --- r (SOCWR) ---
        if hasfield(VT, :r)
            s === :r_consensus      && return getfield(vars, :r).consensus
            s === :r_consensus_prev && return getfield(vars, :r).consensus_prev
            s === :duals_r          && return getfield(vars, :r).duals
            s === :rho_r            && return getfield(vars, :r).rho
            s === :prev_rho_r       && return getfield(vars, :r).prev_rho
        end
        # --- t (SOCWR) ---
        if hasfield(VT, :t)
            s === :t_consensus      && return getfield(vars, :t).consensus
            s === :t_consensus_prev && return getfield(vars, :t).consensus_prev
            s === :duals_t          && return getfield(vars, :t).duals
            s === :rho_t            && return getfield(vars, :t).rho
            s === :prev_rho_t       && return getfield(vars, :t).prev_rho
        end
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
