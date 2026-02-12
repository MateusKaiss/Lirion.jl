# ============================================================
# src/solve.jl
# ============================================================

# -----------------------------
# Model builders
# -----------------------------
dc() = DC
ac() = AC
socp() = SOCP

# -----------------------------
# Algorithm builders
# -----------------------------
centralized() = :centralized

@inline function _admm_cfg(method::Symbol; max_iter::Int = 350, kwargs...)
    return DistributedOPFModel(DC, method, max_iter; kwargs...)
end

admm_bus(; kwargs...) = _admm_cfg(:admm; kwargs...)
admm_branch(; kwargs...) = _admm_cfg(:admm_branch; kwargs...)
admm_hybrid(; kwargs...) = _admm_cfg(:admm_hybrid; kwargs...)

# -----------------------------
# Plot helpers
# -----------------------------
@inline function _input_basedir(path::AbstractString)::String
    p = abspath(String(path))
    return isdir(p) ? p : dirname(p)
end

@inline function _plots_outdir(input_path::AbstractString)::String
    return joinpath(_input_basedir(input_path), "plots_outputs")
end

@inline function _plots_stem(model_type)::String
    if model_type isa DistributedOPFModel
        return "admm_" * string(model_type.base_model) * "_" * string(model_type.method)
    else
        return "opf_" * lowercase(string(model_type))
    end
end

@inline function _maybe_build_plots!(
    write_plots::Bool,
    input_path::AbstractString,
    model_type,
    out::LirionRunResult,
)
    write_plots || return nothing

    diag = out.diagnostics
    if diag isa ADMMHistory
        outdir = _plots_outdir(input_path)
        stem = _plots_stem(model_type)
        mkpath(outdir)
        build_plots(diag; outdir = outdir, stem = stem)
    else
        @warn "Plots requested, but there is no ADMMHistory in this run (no plots generated)."
    end
    return nothing
end
