function _solve_zones_parallel!(ws, post_solve_hook)
    parts = ws.partitions.parts

    Threads.@threads for i in eachindex(parts)
        part = parts[i]
        model = part.model

        update_objective!(part, ws)
        JuMP.optimize!(model)
        _require_optimal!(model; iter = ws.iter, zone_id = part.id)

        part.has_solution = true
        part.generation_cost = JuMP.value(model[:generation_cost])

        post_solve_hook(ws, i, part)
    end

    return nothing
end
