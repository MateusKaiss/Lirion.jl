module PlotlyLightExt

using Lirion
using PlotlyLight

# ============================================================
# Primitives
# ============================================================

function _plot_convergence(
    res_primal::Vector{Float64},
    res_dual::Vector{Float64},
    eps_primal::Vector{Float64},
    eps_dual::Vector{Float64},
    path::String,
)
    x = 1:length(res_primal)
    p = Plot()
    p(; x, y = res_primal, type = :scatter, mode = "lines",       name = "‖rᵏ‖")
    p(; x, y = res_dual,   type = :scatter, mode = "lines",       name = "‖sᵏ‖")
    p(; x, y = eps_primal, type = :scatter, mode = "lines",       name = "epsᵖʳⁱ",  line = (; dash = "dash"))
    p(; x, y = eps_dual,   type = :scatter, mode = "lines",       name = "epsᵈᵘᵃˡ", line = (; dash = "dash"))
    p.layout.title         = "Convergência do ADMM"
    p.layout.xaxis.title   = "Iterações k"
    p.layout.yaxis.title   = "Resíduo"
    p.layout.yaxis.type    = "log"
    PlotlyLight.save(p, path)
    return nothing
end

function _plot_rho(rho_vec::Vector{Float64}, path::String; label::String = "rho")
    x = 1:length(rho_vec)
    p = Plot()
    p(; x, y = rho_vec, type = :scatter, mode = "lines+markers", name = label)
    p.layout.title       = "Evolução de $(label) no ADMM"
    p.layout.xaxis.title = "Iterações k"
    p.layout.yaxis.title = label
    p.layout.yaxis.type  = "log"
    PlotlyLight.save(p, path)
    return nothing
end

function _plot_theta_avg_heatmap(theta_avg_history::Dict{Int,Vector{Float64}}, path::String)
    isempty(theta_avg_history) && return nothing
    buses    = sort!(collect(keys(theta_avg_history)))
    z_matrix = [theta_avg_history[b] for b in buses]
    (isempty(z_matrix) || isempty(z_matrix[1])) && return nothing
    p = Plot()
    p(; z = z_matrix, x = 1:length(z_matrix[1]), y = buses,
       type = :heatmap, colorscale = "Viridis", colorbar = (; title = "theta_avg"))
    p.layout.title       = "Evolução de theta_avg por Barra Compartilhada"
    p.layout.xaxis.title = "Iterações"
    p.layout.yaxis.title = "Ângulo de fase (rad) por Barra Compartilhada"
    PlotlyLight.save(p, path)
    return nothing
end

function _plot_theta_local_heatmaps(
    theta_local_history::Dict{Int,Dict{Int,Vector{Float64}}},
    base_path::String,
)
    isempty(theta_local_history) && return nothing
    for (zone_id, zone_data) in theta_local_history
        isempty(zone_data) && continue
        buses    = sort!(collect(keys(zone_data)))
        z_matrix = [zone_data[b] for b in buses]
        (isempty(z_matrix) || isempty(z_matrix[1])) && continue
        p = Plot()
        p(; z = z_matrix, x = 1:length(z_matrix[1]), y = buses,
           type = :heatmap, colorscale = "Viridis", colorbar = (; title = "theta (local)"))
        p.layout.title       = "Zona $zone_id – theta local por Barra Local"
        p.layout.xaxis.title = "Iterações"
        p.layout.yaxis.title = "Ângulo de fase (rad) por Barra Local"
        PlotlyLight.save(p, "$(base_path)_zone_$zone_id.html")
    end
    return nothing
end

# ============================================================
# Legacy wrappers
# ============================================================

function Lirion.save_dc_bus_plots(
    res_pri::Vector{Float64},
    res_du::Vector{Float64},
    eps_pri::Vector{Float64},
    eps_du::Vector{Float64},
    rho_hist::Vector{Float64},
    theta_avg_history::Dict{Int,Vector{Float64}};
    stem::String = "admm",
)
    _plot_convergence(res_pri, res_du, eps_pri, eps_du, "plot_convergencia_$(stem)_dc.html")
    _plot_rho(rho_hist, "plot_rho_$(stem)_dc.html"; label = "rho")
    _plot_theta_avg_heatmap(theta_avg_history, "plot_theta_avg_heatmap_$(stem)_dc.html")
    return nothing
end

function Lirion.save_dc_branch_plots(
    res_pri::Vector{Float64},
    res_du::Vector{Float64},
    eps_pri::Vector{Float64},
    eps_du::Vector{Float64},
    rho_hist::Vector{Float64};
    stem::String = "admm_branch",
)
    _plot_convergence(res_pri, res_du, eps_pri, eps_du, "plot_convergencia_$(stem).html")
    _plot_rho(rho_hist, "plot_rho_$(stem).html"; label = "rho")
    return nothing
end

# ============================================================
# Main API
# ============================================================

function Lirion.build_plots(
    hist::Lirion.ADMMHistory;
    outdir::String = ".",
    stem::String   = "admm",
)
    mkpath(outdir)

    if haskey(hist.meta, :tag) && hist.meta[:tag] isa AbstractString && !isempty(hist.meta[:tag])
        stem = String(hist.meta[:tag])
    elseif haskey(hist.meta, :kind)
        stem = stem * "_" * string(hist.meta[:kind])
    end

    if !isempty(hist.r_primal)
        _plot_convergence(
            hist.r_primal, hist.r_dual, hist.eps_primal, hist.eps_dual,
            joinpath(outdir, "plot_convergencia_$(stem).html"),
        )
    end

    for (k, v) in hist.rhos
        isempty(v) && continue
        _plot_rho(v, joinpath(outdir, "plot_$(String(k))_$(stem).html"); label = String(k))
    end

    if haskey(hist.heatmaps, :theta_avg)
        _plot_theta_avg_heatmap(
            hist.heatmaps[:theta_avg],
            joinpath(outdir, "plot_theta_avg_heatmap_$(stem).html"),
        )
    end

    if haskey(hist.heatmaps, :theta_local)
        _plot_theta_local_heatmaps(
            hist.heatmaps[:theta_local],
            joinpath(outdir, "plot_theta_local_heatmaps_$(stem)"),
        )
    end

    return nothing
end

end # module
