# -----------------------------
# Model builders
# -----------------------------

dc() = DC
ac() = AC
acr() = ACR
ivr() = IVR
lpac() = LPAC
socwr() = SOCWR

# -----------------------------
# Algorithm builders
# -----------------------------

centralized() = CentralizedAlgorithm()

@inline _sectioning_name(s::Sectioning) = lowercase(string(nameof(typeof(s))))

@inline function _sectioning(x)::Sectioning
    x isa Sectioning && return x
    x === :tie_line && return TieLine()
    x === :dummy_bus && return DummyBus()
    x === :metis && return Metis()
    throw(ArgumentError("Invalid sectioning=$x. Allowed: :tie_line, :dummy_bus, :metis"))
end

@inline function admm(;
    coupling::Symbol = :hybrid,
    sectioning = :tie_line,
    max_iter::Int = 350,
    kwargs...,
)
    coupling in (:state, :power, :hybrid) ||
        throw(ArgumentError("Invalid coupling=$coupling. Allowed: :state, :power, :hybrid"))
    max_iter > 0 || throw(ArgumentError("max_iter must be > 0"))

    return ADMMAlgorithm(coupling, _sectioning(sectioning), max_iter, (; kwargs...))
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
        return "admm_" *
               lowercase(string(nameof(model_type.base_model))) * "_" *
               string(model_type.coupling) * "_" *
               _sectioning_name(model_type.sectioning)
    else
        return "opf_" * lowercase(string(nameof(model_type)))
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
