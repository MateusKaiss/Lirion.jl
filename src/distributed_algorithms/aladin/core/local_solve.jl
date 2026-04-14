# ============================================================
# ALADIN — Phase 1 local solve
# ============================================================

function aladin_local_solve!(ws::AbstractALADINWorkspace)
    parts = ws.partitions.parts

    Threads.@threads for i in eachindex(parts)
        part = parts[i]
        zone = ws.zones[i]
        model = part.model

        update_objective!(part, ws)
        JuMP.optimize!(model)
        _require_optimal!(model; iter = ws.iter, zone_id = part.id)

        part.has_solution = true
        part.generation_cost = JuMP.value(model[:generation_cost])

        extract_solve_values!(zone, model)

        if ws.iter == 1
            copyto!(zone.x_prev, zone.y)
        end

        compute_cost_gradient!(zone)
        refresh_lagrangian_hessian!(zone, ws)
        detect_active_bounds!(zone)
    end
    return nothing
end

# Default dispatch: DC zones have a constant Hessian baked into
# `hess_diag` at build time, so nothing to refresh.
@inline refresh_lagrangian_hessian!(::AladinZoneData, ::AbstractALADINWorkspace) =
    nothing
