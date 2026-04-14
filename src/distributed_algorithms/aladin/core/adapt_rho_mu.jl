# ============================================================
# ALADIN — Adaptive ρ/μ schedule
# ------------------------------------------------------------
# Engineering extension on top of canonical HFD 2016: both penalties
# grow monotonically geometrically until their respective caps, which
# empirically accelerates convergence on hard instances without
# affecting the algorithm's correctness.
#
#     ρ ← min(rho_growth · ρ, rho_upper)
#     μ ← min(mu_growth  · μ, mu_upper)
#
# Disabled when params.adaptive == false (then both are fixed).
# ============================================================

@inline function _grow_ref!(
    ref::Base.RefValue{Float64},
    factor::Float64,
    upper::Float64,
)
    ref[] = min(factor * ref[], upper)
    return nothing
end

# ---- State / Power (single ρ) ----

function adapt_rho_mu!(
    ws::Union{DCStateWorkspaceALADIN, DCPowerWorkspaceALADIN},
    params::ALADINParams,
)
    params.adaptive || return nothing
    _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    return nothing
end

# ---- Hybrid (ρ_bus + ρ_edge) ----

function adapt_rho_mu!(ws::DCHybridWorkspaceALADIN, params::ALADINParams)
    params.adaptive || return nothing
    _grow_ref!(ws.rho_bus, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.rho_edge, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    return nothing
end
