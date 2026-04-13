function local_solve!(ws)
    parts = ws.partitions.parts

    Threads.@threads for i in eachindex(parts)
        part = parts[i]
        model = part.model

        update_objective!(part, ws)

        optimize!(model)

        _require_optimal!(model; iter = ws.iter, zone_id = part.id)

        part.has_solution = true

        part.generation_cost = value(model[:generation_cost])
    end

    return nothing
end
