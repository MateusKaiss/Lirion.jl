# ============================================================
# ALADIN — Phase 2 primal update
# ============================================================

function primal_update!(ws::AbstractALADINWorkspace)
    a1 = ws.alpha1
    a2 = ws.alpha2
    for zone in ws.zones
        x_prev = zone.x_prev
        y = zone.y
        delta = zone.delta
        @inbounds for k in 1:zone.n_vars
            x_prev[k] = x_prev[k] + a1 * (y[k] - x_prev[k]) + a2 * delta[k]
        end
    end
    return nothing
end
