# ============================================================
# Shared
# ============================================================

function _solve_zones_parallel!(ws, post_solve_hook)
    parts = ws.partitions.parts

    Threads.@threads for i in eachindex(parts)
        part = parts[i]
        model = part.local_model

        update_objective!(part, ws)
        JuMP.optimize!(model)
        _require_optimal!(model; iter = ws.iter, zone_id = part.id)

        part.generation_cost = JuMP.value(model[:generation_cost])

        post_solve_hook(ws, i, part)
    end

    return nothing
end

# ============================================================
# ADMM
# ============================================================

@inline _admm_post_solve!(::AbstractADMMWorkspace, ::Integer, ::Any) = nothing

local_solve!(ws::AbstractADMMWorkspace) =
    _solve_zones_parallel!(ws, _admm_post_solve!)

# ============================================================
# ALADIN
# ============================================================

function _aladin_post_solve!(ws::AbstractALADINWorkspace, i::Integer, part)
    zone = ws.zones[i]
    _note_almost_zone_solve!(part.local_model)   # count only; accept logic unchanged
    extract_solve_values!(zone, part.local_model)

    if ws.iter == 1
        copyto!(zone.x_prev, zone.y)
    end

    compute_cost_gradient!(zone)
    detect_active_bounds!(zone)
    refresh_lagrangian_hessian!(zone, ws)
    return nothing
end

local_solve!(ws::AbstractALADINWorkspace) =
    _solve_zones_parallel!(ws, _aladin_post_solve!)
