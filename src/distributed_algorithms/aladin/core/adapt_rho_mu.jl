# ============================================================
# ALADIN — Adaptive ρ/μ schedule
# ============================================================
#
# EXTENSION OUTSIDE THE CANONICAL PAPER. Houska-Frasch-Diehl 2016
# specifies a single fixed ρ "sufficiently large" (Lemma 3 / Theorem 4
# hypotheses); no schedule is prescribed. The monotonic growth here is a
# practical heuristic borrowed from augmented Lagrangian / penalty
# methods. It remains compatible with the paper's local convergence
# guarantee provided ρ_upper < ∞ (so ρ stabilises in finite iterations).
#
# Default in `ALADINConfig` is `adaptive = false`. Enable explicitly
# (`adaptive = true`) with caution; the choice of `rho_growth`,
# `mu_growth` and the upper bounds is empirical and may need tuning per
# case. The adaptive schedule is NOT covered by the convergence proofs.

@inline function _grow_ref!(
    ref::Base.RefValue{Float64},
    factor::Float64,
    upper::Float64,
)
    prev = ref[]
    ref[] = min(factor * prev, upper)
    return ref[] != prev
end

# ---- DC State (single ρ) — invalidate quad_static if ρ changed ----

function adapt_rho_mu!(ws::DCStateWorkspaceALADIN, params::ALADINParams)
    params.adaptive || return nothing
    rho_changed = _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    if rho_changed
        @inbounds for (i, zone) in enumerate(ws.zones)
            _build_quad_static_primal_dc!(
                zone, ws.rho[];
                sigma_theta_shared = ws.sigma_theta_shared,
                sigma_pg = ws.sigma_pg,
                sigma_other = ws.sigma_other,
            )
            _set_obj_static_dc!(ws.partitions.parts[i], zone)
        end
    end
    return nothing
end

# ---- DC Power (single ρ, primal or pflow) ----

function adapt_rho_mu!(ws::DCPowerWorkspaceALADIN, params::ALADINParams)
    params.adaptive || return nothing
    rho_changed = _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    if rho_changed
        @inbounds for (i, zone) in enumerate(ws.zones)
            if ws.proximal_mode === :pflow
                _build_quad_static_pflow_dc!(
                    zone, ws.rho[];
                    sigma_pflow = ws.sigma_pflow,
                    sigma_pg = ws.sigma_pg,
                    sigma_other = ws.sigma_other,
                )
            else
                _build_quad_static_primal_dc!(
                    zone, ws.rho[];
                    sigma_theta_shared = ws.sigma_theta_shared,
                    sigma_pg = ws.sigma_pg,
                    sigma_other = ws.sigma_other,
                )
            end
            _set_obj_static_dc!(ws.partitions.parts[i], zone)
        end
    end
    return nothing
end

# ---- DC Hybrid (ρ_bus + ρ_edge + ρ_inj; quad_static keyed on ρ_bus) ----

function adapt_rho_mu!(ws::DCHybridWorkspaceALADIN, params::ALADINParams)
    params.adaptive || return nothing
    rho_bus_changed = _grow_ref!(ws.rho_bus, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.rho_edge, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.rho_inj, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    if rho_bus_changed
        @inbounds for (i, zone) in enumerate(ws.zones)
            _build_quad_static_primal_dc!(
                zone, ws.rho_bus[];
                sigma_theta_shared = ws.sigma_theta_shared,
                sigma_pg = ws.sigma_pg,
                sigma_other = ws.sigma_other,
            )
            _set_obj_static_dc!(ws.partitions.parts[i], zone)
        end
    end
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
    rho_changed = _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    if rho_changed
        _refresh_quad_static_ac_all!(ws.zones, ws.rho[], ws)
    end
    return nothing
end

# ---- SOCWR (State / Hybrid: primal cache; Power: depends on proximal_mode) ----

function adapt_rho_mu!(
    ws::Union{SOCWRStateWorkspaceALADIN, SOCWRHybridWorkspaceALADIN},
    params::ALADINParams,
)
    params.adaptive || return nothing
    rho_changed = _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    if rho_changed
        _refresh_quad_static_socwr_primal_all!(ws.partitions.parts, ws.zones, ws.rho[], ws)
    end
    return nothing
end

function adapt_rho_mu!(ws::SOCWRPowerWorkspaceALADIN, params::ALADINParams)
    params.adaptive || return nothing
    rho_changed = _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    if rho_changed
        if ws.proximal_mode === :pflow
            _refresh_quad_static_socwr_pflow_all!(
                ws.partitions.parts,
                ws.zones,
                ws.edge_pairs,
                ws.rho[],
                ws,
            )
        else
            _refresh_quad_static_socwr_primal_all!(
                ws.partitions.parts,
                ws.zones,
                ws.rho[],
                ws,
            )
        end
    end
    return nothing
end
