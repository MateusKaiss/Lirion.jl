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
        name = "rho",
        line = (; color = "green"),
    )

    p.layout.title = "Evolução de rho (rho) no ADMM"
    p.layout.xaxis.title = "Iterações k"
    p.layout.yaxis.title = "rho"
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


##########################
# WRAPPERS DE PLOT       #
##########################

function save_dc_bus_plots(
    res_pri,
    res_du,
    eps_pri,
    eps_du,
    rho_hist,
    θ_avg_history;
    stem::String = "admm",
)
    plot_admm_convergence(
        res_pri,
        res_du,
        eps_pri,
        eps_du,
        "plot_convergencia_$(stem)_dc.html",
    )
    plot_rho(rho_hist, "plot_rho_$(stem)_dc.html")
    return plot_theta_avg_heatmap(θ_avg_history, "plot_theta_avg_heatmap_$(stem)_dc.html")
end

function save_dc_branch_plots(
    res_pri,
    res_du,
    eps_pri,
    eps_du,
    rho_hist;
    stem::String = "admm_branch",
)
    plot_admm_convergence(
        res_pri,
        res_du,
        eps_pri,
        eps_du,
        "plot_convergencia_$(stem).html",
    )
    return plot_rho(rho_hist, "plot_rho_$(stem).html")
end

##########################
# WRAPPERS DE PLOT       #
##########################

function save_ac_bus_plots(
    res_pri,
    res_du,
    eps_pri,
    eps_du,
    rhoθ_hist,
    rhoV_hist,
    θ_avg_hist;
    stem::String = "admm_ac",
)
    plot_admm_convergence(
        res_pri,
        res_du,
        eps_pri,
        eps_du,
        "plot_convergencia_$(stem).html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoθ_hist),
                    y = rhoθ_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "rhoθ",
                );
                p.layout.title = "Evolução de rhoθ";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "rhoθ";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rho_theta_$(stem).html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoV_hist),
                    y = rhoV_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "rhoV",
                );
                p.layout.title = "Evolução de rhoV";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "rhoV";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rho_vmag_$(stem).html",
    )
    return plot_theta_avg_heatmap(θ_avg_hist, "plot_theta_avg_heatmap_$(stem).html")
end

function save_ac_branch_plots(
    res_pri,
    res_du,
    eps_pri,
    eps_du,
    rhoP_hist,
    rhoQ_hist;
    stem::String = "admm_ac_branch",
)
    plot_admm_convergence(
        res_pri,
        res_du,
        eps_pri,
        eps_du,
        "plot_convergencia_$(stem).html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoP_hist),
                    y = rhoP_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "rhoP",
                );
                p.layout.title = "Evolução de rhoP";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "rhoP";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rhoP_$(stem).html",
    )
    return PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoQ_hist),
                    y = rhoQ_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "rhoQ",
                );
                p.layout.title = "Evolução de rhoQ";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "rhoQ";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rhoQ_$(stem).html",
    )
end

# ---------- Plots (wrappers) ----------
function save_socp_v2rt_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
    rhoV2_hist, rhoR_hist, rhoT_hist; stem = "admm_socp_v2rt")
    plot_admm_convergence(r_hist, s_hist, eps_r_hist, eps_s_hist,
        "plot_convergencia_$(stem).html")
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(; x = 1:length(rhoV2_hist), y = rhoV2_hist, type = :scatter,
                    mode = "lines+markers", name = "rho(V²)");
                p.layout.title = "Evolução de rho(V²)";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "rho(V²)";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rho_v2_$(stem).html")
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(; x = 1:length(rhoR_hist), y = rhoR_hist, type = :scatter,
                    mode = "lines+markers", name = "rhoR");
                p.layout.title = "Evolução de rhoR";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "rhoR";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rhoR_$(stem).html")
    return PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(; x = 1:length(rhoT_hist), y = rhoT_hist, type = :scatter,
                    mode = "lines+markers", name = "rhoT");
                p.layout.title = "Evolução de rhoT";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "rhoT";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rhoT_$(stem).html")
end

function save_socp_branch_pq_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
    rhoP_hist, rhoQ_hist; stem = "admm_socp_branch")
    plot_admm_convergence(r_hist, s_hist, eps_r_hist, eps_s_hist,
        "plot_convergencia_$(stem).html")
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(; x = 1:length(rhoP_hist), y = rhoP_hist, type = :scatter,
                    mode = "lines+markers", name = "rhoP");
                p.layout.title = "Evolução de rhoP";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "rhoP";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rhoP_$(stem).html")
    return PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(; x = 1:length(rhoQ_hist), y = rhoQ_hist, type = :scatter,
                    mode = "lines+markers", name = "rhoQ");
                p.layout.title = "Evolução de rhoQ";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "rhoQ";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rhoQ_$(stem).html")
end

function save_socp_hybrid_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
    rhoV2_hist, rhoP_hist, rhoQ_hist, rhoR_hist, rhoT_hist;
    stem = "admm_socp_hybrid")
    plot_admm_convergence(r_hist, s_hist, eps_r_hist, eps_s_hist,
        "plot_convergencia_$(stem).html")
    for (lab, hist) in [("rho(V²)", rhoV2_hist), ("rhoP", rhoP_hist), ("rhoQ", rhoQ_hist),
        ("rhoR", rhoR_hist), ("rhoT", rhoT_hist)]
        fname = replace(lab, "/" => "") |> x -> lowercase(replace(x, "rho" => "rho"))
        PlotlyLight.save(
            (
                () -> (
                    p = Plot();
                    p(; x = 1:length(hist), y = hist, type = :scatter,
                        mode = "lines+markers", name = lab);
                    p.layout.title = "Evolução de $(lab)";
                    p.layout.xaxis.title = "Iterações";
                    p.layout.yaxis.title = lab;
                    p.layout.yaxis.type = "log";
                    p
                )
            )(), "plot_$(fname)_$(stem).html")
    end
end
