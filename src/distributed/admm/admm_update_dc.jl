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
                model = pdata[:model]
                θ = value.(model[:θ])
                idx = findfirst(isequal(b), pdata["bus"]["bus_i"])
                θ_sum += θ[idx]
                count += 1
            end
        end
        θ_avg[b] = θ_sum / count
    end

    # Update the λ (Lagrange multipliers)
    for (zone_id, pdata) in partitions
        model = pdata[:model]
        θ = value.(model[:θ])
        bus_ids = pdata["bus"]["bus_i"]

        for b in shared_buses
            if b in bus_ids
                idx = findfirst(isequal(b), bus_ids)
                λ = get(duals[zone_id], b, 0.0)
                λ += ρ * (θ[idx] - θ_avg[b])
                duals[zone_id][b] = λ
            end
        end
    end
end