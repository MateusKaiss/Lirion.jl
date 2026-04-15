function _aladin_post_solve!(ws::AbstractALADINWorkspace, i::Integer, part)
    zone = ws.zones[i]
    extract_solve_values!(zone, part.model)

    if ws.iter == 1
        copyto!(zone.x_prev, zone.y)
    end

    compute_cost_gradient!(zone)
    refresh_lagrangian_hessian!(zone, ws)
    detect_active_bounds!(zone)
    return nothing
end

local_solve!(ws::AbstractALADINWorkspace) =
    _solve_zones_parallel!(ws, _aladin_post_solve!)
