function update_shared_variables_dc!(
    shared_buses::Vector{Int},
    partitions::Dict{Int, Dict{String, Any}},
    duals::Dict{Int, Dict{Int, Float64}},
    θ_avg::Dict{Int, Float64},
    ρ::Float64,
)
    for b in shared_buses
        θ_sum = 0.0
        count = 0
        for (_, pdata) in partitions
            if b in pdata["bus"]["bus_i"]
                model = pdata["model"]  # ← CORRIGIDO
                θ_val = value(model[:theta][b, 1])
                θ_sum += θ_val
                count += 1
            end
        end
        θ_avg[b] = θ_sum / count
    end

    # Update the λ (Lagrange multipliers)
    for (zone_id, pdata) in partitions
        model = pdata["model"]
        θ = value.(model[:theta])
        bus_ids = pdata["bus"]["bus_i"]

        for b in shared_buses
            if b in bus_ids
                λ = get(duals[zone_id], b, 0.0)
                λ += ρ * (θ[b, 1] - θ_avg[b])  # ← CORRIGIDO: usa b diretamente
                duals[zone_id][b] = λ
            end
        end
    end
end
