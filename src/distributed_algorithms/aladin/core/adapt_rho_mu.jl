# ============================================================
# ALADIN — Adaptive ρ/μ schedule
# ============================================================
#
# Called once per outer iteration (after residuals + line search).
# Grows ρ and μ monotonically (`adaptive = true`) and, when ρ moved,
# rebuilds each zone's cached `quad_static` via `refresh_rho!`.
#
# Monotonic growth is a practical heuristic borrowed from augmented
# Lagrangian / penalty methods. Houska-Frasch-Diehl 2016 fixes ρ
# "sufficiently large" and proves local convergence; the growth here
# is compatible with that proof as long as `rho_upper < ∞` (ρ
# stabilises in finite iterations). The schedule itself is NOT
# covered by the proof — defaults to `adaptive = false` in
# `ALADINConfig`; enable explicitly with caution.
#
# Note μ never triggers `refresh_rho!`: it enters only the coordinator
# QP (Levenberg-Marquardt regularisation) and never appears inside
# `quad_static`, so bumping `ws.mu[]` is free.

# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------
#
# Multiplies `ref[]` by `factor`, clamps at `upper`, and reports
# whether the value actually moved. Callers use the return value to
# skip cache rebuilds when ρ has saturated at the upper bound.

@inline function _grow_ref!(
    ref::Base.RefValue{Float64},
    factor::Float64,
    upper::Float64,
)
    prev = ref[]
    ref[] = min(factor * prev, upper)
    return ref[] != prev
end

# ============================================================
# Dispatch
# ============================================================

# ------------------------------------------------------------
# DC State (single ρ)
# ------------------------------------------------------------

function adapt_rho_mu!(ws::DCStateWorkspaceALADIN, params::ALADINParams)
    params.adaptive || return nothing
    rho_changed = _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    rho_changed && refresh_rho!(ws)
    return nothing
end

# ------------------------------------------------------------
# DC Power (single ρ; primal or pflow proximal mode)
# ------------------------------------------------------------

function adapt_rho_mu!(ws::DCPowerWorkspaceALADIN, params::ALADINParams)
    params.adaptive || return nothing
    rho_changed = _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    rho_changed && refresh_rho!(ws)
    return nothing
end

# ------------------------------------------------------------
# DC Hybrid (ρ_bus + ρ_edge + ρ_inj)
# ------------------------------------------------------------
#
# Three separate ρ refs — one per consensus class. `quad_static`
# only caches the ρ_bus part (the other two enter the coordinator QP
# as linear couplings, never inside the per-zone quadratic), so only
# a `ρ_bus` change forces `refresh_rho!`.

function adapt_rho_mu!(ws::DCHybridWorkspaceALADIN, params::ALADINParams)
    params.adaptive || return nothing
    rho_bus_changed = _grow_ref!(ws.rho_bus, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.rho_edge, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.rho_inj, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    rho_bus_changed && refresh_rho!(ws)
    return nothing
end

# ------------------------------------------------------------
# AC (State / Power / Hybrid)
# ------------------------------------------------------------
#
# All three AC workspaces share a single ρ — variants differ only in
# which consensus pairs exist, not in how ρ is scheduled.

function adapt_rho_mu!(
    ws::Union{
        ACStateWorkspaceALADIN,
        ACPowerWorkspaceALADIN,
        ACHybridWorkspaceALADIN,
    },
    params::ALADINParams,
)
    params.adaptive || return nothing
    rho_changed = _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    rho_changed && refresh_rho!(ws)
    return nothing
end

# ------------------------------------------------------------
# SOCWR (State / Power / Hybrid)
# ------------------------------------------------------------
#
# Same pattern as AC: one ρ across all three variants.

function adapt_rho_mu!(
    ws::Union{
        SOCWRStateWorkspaceALADIN,
        SOCWRPowerWorkspaceALADIN,
        SOCWRHybridWorkspaceALADIN,
    },
    params::ALADINParams,
)
    params.adaptive || return nothing
    rho_changed = _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    rho_changed && refresh_rho!(ws)
    return nothing
end
