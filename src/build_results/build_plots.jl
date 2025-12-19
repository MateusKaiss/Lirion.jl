# ============================================================
# build_results/build_plots.jl
# ============================================================
# API ÚNICA: build_plots(hist; outdir=".", stem="admm")
# ============================================================

# (se você já dá `using PlotlyLight` no módulo, não precisa aqui)
# using PlotlyLight

# ----------------------------
# Plots básicos
# ----------------------------

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
    )
    p(;
        x = x_vals,
        y = res_dual,
        type = :scatter,
        mode = "lines",
        name = "‖sᵏ‖",
    )
    p(;
        x = x_vals,
        y = eps_primal,
        type = :scatter,
        mode = "lines",
        name = "epsᵖʳⁱ",
        line = (; dash = "dash"),
    )
    p(;
        x = x_vals,
        y = eps_dual,
        type = :scatter,
        mode = "lines",
        name = "epsᵈᵘᵃˡ",
        line = (; dash = "dash"),
    )

    p.layout.title = "Convergência do ADMM"
    p.layout.xaxis.title = "Iterações k"
    p.layout.yaxis.title = "Resíduo"
    p.layout.yaxis.type = "log"

    PlotlyLight.save(p, path)
    return nothing
end

function plot_rho(rho_vec::Vector{Float64}, path::String; label::String = "rho")
    x_vals = 1:length(rho_vec)
    p = Plot()
    p(;
        x = x_vals,
        y = rho_vec,
        type = :scatter,
        mode = "lines+markers",
        name = label,
    )

    p.layout.title = "Evolução de $(label) no ADMM"
    p.layout.xaxis.title = "Iterações k"
    p.layout.yaxis.title = label
    p.layout.yaxis.type = "log"

    PlotlyLight.save(p, path)
    return nothing
end

function plot_theta_avg_heatmap(theta_avg_history::Dict{Int, Vector{Float64}}, path::String)
    isempty(theta_avg_history) && return nothing

    buses = sort!(collect(keys(theta_avg_history)))
    z_matrix = [theta_avg_history[b] for b in buses]
    isempty(z_matrix) && return nothing
    isempty(z_matrix[1]) && return nothing

    p = Plot()
    p(;
        z = z_matrix,
        x = 1:length(z_matrix[1]),
        y = buses,
        type = :heatmap,
        colorscale = "Viridis",
        colorbar = (; title = "theta_avg"),
    )

    p.layout.title = "Evolução de theta_avg por Barra Compartilhada"
    p.layout.xaxis.title = "Iterações"
    p.layout.yaxis.title = "Ângulo de fase (rad) por Barra Compartilhada"

    PlotlyLight.save(p, path)
    return nothing
end

function plot_theta_local_heatmaps(
    theta_local_history::Dict{Int, Dict{Int, Vector{Float64}}},
    base_path::String,
)
    isempty(theta_local_history) && return nothing

    for (zone_id, zone_data) in theta_local_history
        isempty(zone_data) && continue
        buses = sort!(collect(keys(zone_data)))
        z_matrix = [zone_data[b] for b in buses]
        isempty(z_matrix) && continue
        isempty(z_matrix[1]) && continue

        p = Plot()
        p(;
            z = z_matrix,
            x = 1:length(z_matrix[1]),
            y = buses,
            type = :heatmap,
            colorscale = "Viridis",
            colorbar = (; title = "theta (local)"),
        )

        p.layout.title = "Zona $zone_id – theta local por Barra Local"
        p.layout.xaxis.title = "Iterações"
        p.layout.yaxis.title = "Ângulo de fase (rad) por Barra Local"

        path = "$(base_path)_zone_$zone_id.html"
        PlotlyLight.save(p, path)
    end
    return nothing
end

# ----------------------------
# Wrappers legados (opcional manter)
# OBS: eles salvam no diretório atual; prefira build_plots(...)
# ----------------------------

function save_dc_bus_plots(
    res_pri,
    res_du,
    eps_pri,
    eps_du,
    rho_hist,
    theta_avg_history;
    stem::String = "admm",
)
    plot_admm_convergence(
        res_pri,
        res_du,
        eps_pri,
        eps_du,
        "plot_convergencia_$(stem)_dc.html",
    )
    plot_rho(rho_hist, "plot_rho_$(stem)_dc.html"; label = "rho")
    plot_theta_avg_heatmap(theta_avg_history, "plot_theta_avg_heatmap_$(stem)_dc.html")
    return nothing
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
    plot_rho(rho_hist, "plot_rho_$(stem).html"; label = "rho")
    return nothing
end

# ----------------------------
# API única (recomendada)
# ----------------------------

"""
Gera plots a partir do ADMMHistory.

Convenções:
- hist.rhos pode conter chaves: :rho, :rho_theta, :rho_v, :rhoP, :rhoQ, :rhoV2, :rhoR, :rhoT...
- hist.heatmaps pode conter: :theta_avg, :theta_local
- hist.meta pode conter: :tag (String) ou :kind (Symbol) para nomear arquivos (opcional)
"""
function build_plots(hist::ADMMHistory; outdir::String = ".", stem::String = "admm")
    mkpath(outdir)

    # 0) stem automático via meta (se existir)
    if haskey(hist.meta, :tag) && hist.meta[:tag] isa AbstractString &&
       !isempty(hist.meta[:tag])
        stem = String(hist.meta[:tag])
    elseif haskey(hist.meta, :kind)
        stem = stem * "_" * string(hist.meta[:kind])
    end

    # 1) convergência
    if !isempty(hist.r_primal)
        plot_admm_convergence(
            hist.r_primal,
            hist.r_dual,
            hist.eps_primal,
            hist.eps_dual,
            joinpath(outdir, "plot_convergencia_$(stem).html"),
        )
    end

    # 2) rhos
    for (k, v) in hist.rhos
        isempty(v) && continue
        plot_rho(
            v,
            joinpath(outdir, "plot_$(String(k))_$(stem).html");
            label = String(k),
        )
    end

    # 3) heatmaps
    if haskey(hist.heatmaps, :theta_avg)
        plot_theta_avg_heatmap(
            hist.heatmaps[:theta_avg],
            joinpath(outdir, "plot_theta_avg_heatmap_$(stem).html"),
        )
    end

    if haskey(hist.heatmaps, :theta_local)
        plot_theta_local_heatmaps(
            hist.heatmaps[:theta_local],
            joinpath(outdir, "plot_theta_local_heatmaps_$(stem)"),
        )
    end

    return nothing
end
