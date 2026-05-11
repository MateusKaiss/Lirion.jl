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
