# -----------------------------
# Model builders (mantém Type)
# -----------------------------
dc() = DC
ac() = AC
socp() = SOCP

# -----------------------------
# Algorithm builders
# -----------------------------
centralized() = :centralized

@inline function admm(;
    coupling::Symbol = :hybrid,         # :state | :power | :hybrid
    sectioning::Symbol = :tie_line,     # por enquanto só :tie_line
    max_iter::Int = 350,
    kwargs...,
)
    coupling in (:state, :power, :hybrid) ||
        throw(ArgumentError("Invalid coupling=$coupling. Allowed: :state, :power, :hybrid"))

    sectioning === :tie_line ||
        throw(
            ArgumentError("Invalid sectioning=$sectioning. Currently supported: :tie_line"),
        )

    return ADMMAlgorithmSpec(coupling, sectioning, max_iter, (; kwargs...))
end

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
        # Ex: admm_ac_hybrid_tie_line
        return "admm_" *
               lowercase(string(model_type.base_model)) * "_" *
               string(model_type.coupling) * "_" *
               string(model_type.sectioning)
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

@inline function _assert_valid_path(path::AbstractString)::String
    p = abspath(String(path))
    (isfile(p) || isdir(p)) ||
        throw(ArgumentError("Invalid path: $p is not a file or directory."))
    return p
end
