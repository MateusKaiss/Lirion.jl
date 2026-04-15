# ============================================================
# ALADIN — Adaptive ρ/μ schedule
# ============================================================

@inline function _grow_ref!(
    ref::Base.RefValue{Float64},
    factor::Float64,
    upper::Float64,
)
    ref[] = min(factor * ref[], upper)
    return nothing
end

# ---- DC State / Power (single ρ) ----

function adapt_rho_mu!(
    ws::Union{DCStateWorkspaceALADIN, DCPowerWorkspaceALADIN},
    params::ALADINParams,
)
    params.adaptive || return nothing
    _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    return nothing
end

# ---- DC Hybrid (ρ_bus + ρ_edge) ----

function adapt_rho_mu!(ws::DCHybridWorkspaceALADIN, params::ALADINParams)
    params.adaptive || return nothing
    _grow_ref!(ws.rho_bus, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.rho_edge, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    return nothing
end

# ---- AC (State / Power / Hybrid) ----

function adapt_rho_mu!(
    ws::Union{
        ACStateWorkspaceALADIN,
        ACPowerWorkspaceALADIN,
        ACHybridWorkspaceALADIN,
    },
    params::ALADINParams,
)
    params.adaptive || return nothing
    _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    return nothing
end

# ---- SOCWR (State / Power / Hybrid) ----

function adapt_rho_mu!(
    ws::Union{
        SOCWRStateWorkspaceALADIN,
        SOCWRPowerWorkspaceALADIN,
        SOCWRHybridWorkspaceALADIN,
    },
    params::ALADINParams,
)
    params.adaptive || return nothing
    _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    return nothing
end
