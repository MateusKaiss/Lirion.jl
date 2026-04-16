# ============================================================
# build_plots.jl — stubs
# ============================================================

const _PLOTLYLIGHT_ERROR = """
Plotting requires PlotlyLight.jl. Load it first:

    using PlotlyLight
"""

function build_plots end

function build_plots(hist::ADMMHistory; outdir::String = ".", stem::String = "admm")
    return error(_PLOTLYLIGHT_ERROR)
end

function save_dc_bus_plots(args...; kwargs...)
    return error(_PLOTLYLIGHT_ERROR)
end

function save_dc_branch_plots(args...; kwargs...)
    return error(_PLOTLYLIGHT_ERROR)
end
