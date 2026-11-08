# ============================================================
# Model builders
# ============================================================

"""
    dc()

DC linear active-power approximation (LP). Hand the result to `solve`'s `model`
keyword: `solve(case; model = dc())`. The other eleven builders (`ac`, `acr`,
`ivr`, `lpac`, `socwr`, `socwrnl`, `qc`, `sdpwrm`, `socbf`, `bfa`, `nfa`) follow
the same pattern; the Theory pages cover what each one solves.
"""
dc() = DC

"AC optimal power flow in polar coordinates (NLP)."
ac() = AC

"AC optimal power flow in rectangular coordinates (NLP)."
acr() = ACR

"Current-voltage rectangular AC formulation (NLP)."
ivr() = IVR

"LP-based AC approximation (QCQP)."
lpac() = LPAC

"W-rectangular second-order-cone relaxation, Jabr (SOCP)."
socwr() = SOCWR

"Smooth-NLP form of the SOCWR relaxation (NLP)."
socwrnl() = SOCWRNL

"Quadratic-convex relaxation (QCQP)."
qc() = QC

"Semidefinite W-matrix relaxation (SDP)."
sdpwrm() = SDPWRM

"Branch-flow second-order-cone relaxation (SOCP)."
socbf() = SOCBF

"Linearized branch-flow approximation / LinDistFlow (LP)."
bfa() = BFA

"Network-flow approximation, active power only (LP)."
nfa() = NFA

# ============================================================
# Algorithm builders
# ============================================================

"""
    centralized()

The default algorithm: solve the whole system as one problem. Pass it to `solve`'s
`algorithm` keyword.
"""
centralized() = CentralizedAlgorithm()

@inline _sectioning_name(s::Sectioning) = lowercase(string(nameof(typeof(s))))

@inline function _sectioning(x)::Sectioning
    x isa Sectioning && return x
    x === :tie_line && return TieLine()
    x === :dummy_bus && return DummyBus()
    x === :split_bus && return SplitBus()
    x === :metis && return Metis()
    throw(
        ArgumentError(
            "Invalid sectioning=$x. Allowed: :tie_line, :dummy_bus, :split_bus, :metis",
        ),
    )
end

"""
    admm(; coupling=:hybrid, sectioning=:tie_line, max_iter=350, rho=..., kwargs...)

Build the ADMM distributed algorithm for `solve`'s `algorithm` keyword. `coupling` is
`:state`, `:power`, or `:hybrid`; `sectioning` is `:tie_line`, `:dummy_bus`, or
`:split_bus`. A single `rho` fills the per-coupling penalties; the ADMM settings page
lists the rest.
"""
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

"""
    aladin(; coupling=:state, sectioning=:tie_line, max_iter=100, rho=..., mu=..., kwargs...)

Build the ALADIN distributed algorithm for `solve`'s `algorithm` keyword. Takes the
local penalty `rho` and the coordinator penalty `mu`, plus `proximal_mode`,
`stopping_criterion`, and the globalization; the ALADIN settings page lists them.
"""
@inline function aladin(;
    coupling::Symbol = :state,
    sectioning = :tie_line,
    max_iter::Int = 100,
    kwargs...,
)
    coupling in (:state, :power, :hybrid) ||
        throw(ArgumentError("Invalid coupling=$coupling. Allowed: :state, :power, :hybrid"))
    max_iter > 0 || throw(ArgumentError("max_iter must be > 0"))

    if haskey(kwargs, :proximal_mode)
        pm = kwargs[:proximal_mode]
        pm in (:primal, :pflow) ||
            throw(ArgumentError(
                "Invalid proximal_mode=$pm. Allowed: :primal, :pflow",
            ))
        pm === :pflow && !(coupling in (:power, :hybrid)) &&
            throw(
                ArgumentError(
                    "proximal_mode=:pflow is supported with coupling=:power or " *
                    ":hybrid (got coupling=$(coupling)). Use one of those, or " *
                    "proximal_mode=:primal.",
                ),
            )
    end
    if haskey(kwargs, :sigma_pflow)
        kwargs[:sigma_pflow] > 0 ||
            throw(ArgumentError("sigma_pflow must be > 0"))
    end
    if haskey(kwargs, :sigma_qflow)
        kwargs[:sigma_qflow] > 0 ||
            throw(ArgumentError("sigma_qflow must be > 0"))
    end
    if haskey(kwargs, :stopping_criterion)
        sc = kwargs[:stopping_criterion]
        sc in (:boyd, :hfd) ||
            throw(ArgumentError(
                "Invalid stopping_criterion=$sc. Allowed: :boyd, :hfd",
            ))
    end

    return ALADINAlgorithm(coupling, _sectioning(sectioning), max_iter, (; kwargs...))
end

# ============================================================
# Plot helpers
# ============================================================

@inline function _input_basedir(path::AbstractString)::String
    p = abspath(String(path))
    return isdir(p) ? p : dirname(p)
end

@inline function _plots_outdir(input_path::AbstractString)::String
    return joinpath(_input_basedir(input_path), "plots_outputs")
end

@inline function _plots_stem(model_type)::String
    if model_type isa ADMMConfig
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
