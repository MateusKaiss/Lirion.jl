# ============================================================
# ALADIN — Adaptive ρ/μ schedule
# ============================================================

# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------

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
# Default: single ρ + μ
# ------------------------------------------------------------

function adapt_rho_mu!(ws::AbstractALADINWorkspace, params::ALADINParams)
    params.adaptive || return nothing
    rho_changed = _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    rho_changed && refresh_rho!(ws)
    return nothing
end

# ------------------------------------------------------------
# DC Hybrid
# ------------------------------------------------------------

function adapt_rho_mu!(ws::DCHybridWorkspaceALADIN, params::ALADINParams)
    params.adaptive || return nothing
    rho_bus_changed = _grow_ref!(ws.rho_bus, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.rho_edge, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.rho_inj, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    rho_bus_changed && refresh_rho!(ws)
    return nothing
end
