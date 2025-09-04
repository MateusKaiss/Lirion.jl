@inline function _wrap_angle_diff(Δ::Float64)
    return atan(sin(Δ), cos(Δ))
end

@inline function _circular_mean(θs::AbstractVector{<:Real})
    sum_sin = 0.0
    sum_cos = 0.0
    @inbounds @simd for θ in θs
        sum_sin += sin(θ)
        sum_cos += cos(θ)
    end
    return atan(sum_sin, sum_cos)
end

function compute_real_generation_cost(model, params)
    T = get(model.ext, :time_periods, 1)
    active_generators = findall(params.generator_status .== 1)

    total_cost = 0.0
    for t in 1:T, g in active_generators
        if haskey(params.generator_costs, t) && haskey(params.generator_costs[t], g)
            pgen = value(model[:p_generation][g, t]) * params.sbase
            coefs = params.generator_costs[t][g]
            degree = length(coefs)
            for i in 1:degree
                total_cost += coefs[i] * pgen^(degree - i)
            end
        end
    end
    return total_cost
end

@inline function _ac_flow_from_exprs_existing(
    model::Model,
    subdata::Dict{String, Any},
    l_local::Int;
    t::Int = 1,
)
    return model[:p_flow_branches][l_local, t], model[:q_flow_branches][l_local, t]
end

@inline function _socp_flow_from_exprs_existing(
    model::Model, l_local::Int; t::Int = 1,
)
    return model[:p_flow_branches][l_local, t], model[:q_flow_branches][l_local, t]
end

@inline function _rt_from_exprs_socp(model::Model, l_local::Int; t::Int = 1)
    return model[:R][l_local, t], model[:T_var][l_local, t]
end

function _dc_flow_expr(model::Model, subdata::Dict{String, Any}, l_local::Int; t::Int = 1)
    l_global = subdata["branch_local_indices"][l_local]
    f = subdata["branch_full"]["f_bus"][l_global]
    tbus = subdata["branch_full"]["t_bus"][l_global]
    x = subdata["branch_full"]["br_x"][l_global]

    θf = model[:theta][f, t]
    θt = model[:theta][tbus, t]
    return (θf - θt) / x
end

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
                model = pdata["model"]
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
                λ += ρ * (θ[b, 1] - θ_avg[b])
                duals[zone_id][b] = λ
            end
        end
    end
end

function plot_admm_convergence(
    res_primal::Vector{Float64},
    res_dual::Vector{Float64},
    eps_primal::Vector{Float64},
    eps_dual::Vector{Float64},
    path::String,
)
    x_vals = 1:length(res_primal)
    p = Plot()
    p(;
        x = x_vals,
        y = res_primal,
        type = :scatter,
        mode = "lines",
        name = "‖rᵏ‖",
        line = (; color = "blue"),
    )
    p(;
        x = x_vals,
        y = res_dual,
        type = :scatter,
        mode = "lines",
        name = "‖sᵏ‖",
        line = (; color = "red"),
    )
    p(;
        x = x_vals,
        y = eps_primal,
        type = :scatter,
        mode = "lines",
        name = "εᵖʳⁱ",
        line = (; color = "purple", dash = "dash"),
    )
    p(;
        x = x_vals,
        y = eps_dual,
        type = :scatter,
        mode = "lines",
        name = "εᵈᵘᵃˡ",
        line = (; color = "orange", dash = "dash"),
    )

    p.layout.title = "Convergência do ADMM"
    p.layout.xaxis.title = "Iterações k"
    p.layout.yaxis.title = "Resíduo"
    p.layout.yaxis.type = "log"

    PlotlyLight.save(p, path)
    return println("Gráfico de convergência salvo em $path")
end

function plot_rho(rho_vec::Vector{Float64}, path::String)
    x_vals = 1:length(rho_vec)
    p = Plot()
    p(;
        x = x_vals,
        y = rho_vec,
        type = :scatter,
        mode = "lines+markers",
        name = "ρ",
        line = (; color = "green"),
    )

    p.layout.title = "Evolução de ρ (rho) no ADMM"
    p.layout.xaxis.title = "Iterações k"
    p.layout.yaxis.title = "ρ"
    p.layout.yaxis.type = "log"

    PlotlyLight.save(p, path)
    return println("Gráfico de rho salvo em $path")
end

function plot_theta_avg_heatmap(θ_avg_history::Dict{Int, Vector{Float64}}, path::String)
    buses = sort(collect(keys(θ_avg_history)))
    z_matrix = [θ_avg_history[b] for b in buses]

    p = Plot()
    p(;
        z = z_matrix,
        x = 1:length(z_matrix[1]),
        y = buses,
        type = :heatmap,
        colorscale = "Viridis",
        colorbar = (; title = "θ_avg"),
    )

    p.layout.title = "Evolução de θ_avg por Barra Compartilhada"
    p.layout.xaxis.title = "Iterações"
    p.layout.yaxis.title = "Ângulo de fase (rad) por Barra Compartilhada"

    PlotlyLight.save(p, path)
    return println("Heatmap salvo em $path")
end

function plot_theta_local_heatmaps(
    θ_local_history::Dict{Int, Dict{Int, Vector{Float64}}},
    base_path::String,
)
    for (zone_id, zone_data) in θ_local_history
        buses = sort(collect(keys(zone_data)))
        z_matrix = [zone_data[b] for b in buses]

        p = Plot()
        p(;
            z = z_matrix,
            x = 1:length(z_matrix[1]),
            y = buses,
            type = :heatmap,
            colorscale = "Viridis",
            colorbar = (; title = "θ (local)"),
        )

        p.layout.title = "Zona $zone_id – θ local por Barra Local"
        p.layout.xaxis.title = "Iterações"
        p.layout.yaxis.title = "Ângulo de fase (rad) por Barra Local"

        path = "$(base_path)_zone_$zone_id.html"
        PlotlyLight.save(p, path)
        println(" Heatmap local salvo em $path")
    end
end
