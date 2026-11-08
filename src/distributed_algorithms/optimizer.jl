# ============================================================
# Dispatch
# ============================================================

@inline base_model_tag(::Type{AC}) = Val(:ac)
@inline base_model_tag(::Type{DC}) = Val(:dc)
@inline base_model_tag(::Type{SOCWR}) = Val(:socwr)
@inline base_model_tag(::Type{SOCWRNL}) = Val(:socwrnl)

# ============================================================
# Utilities
# ============================================================

is_gurobi(m::JuMP.Model)::Bool = occursin("Gurobi", MOI.get(m, MOI.SolverName()))

# ============================================================
# Default factories
# ============================================================

# ------------------------------------------------------------
# AC — Ipopt
# ------------------------------------------------------------

const _SILENT_FACTORY_AC = optimizer_with_attributes(
    Ipopt.Optimizer,
    # silent output
    "print_level" => 0,
    "sb" => "yes",
    "print_timing_statistics" => "no",
    "print_user_options" => "no",
    "file_print_level" => 0,
    # convergence criteria
    "tol" => 1e-7,
    "acceptable_tol" => 1e-5,
    "acceptable_iter" => 10,
    # raised from the Ipopt default (3000); propagates to SOCWRNL via _socwrnl_tight_factory
    "max_iter" => 20000,
    "dual_inf_tol" => 1e-4,
    "constr_viol_tol" => 1e-4,
    "compl_inf_tol" => 1e-4,
    "acceptable_dual_inf_tol" => 1e-2,
    "acceptable_constr_viol_tol" => 1e-3,
    "acceptable_compl_inf_tol" => 1e-3,
    # NLP & bounds
    "fixed_variable_treatment" => "make_parameter_nodual",
    "bound_relax_factor" => 1e-6,
    "check_derivatives_for_naninf" => "no",
    "jacobian_approximation" => "exact",
    "gradient_approximation" => "exact",
    "kappa_d" => 1e-5,
    # scaling
    "nlp_scaling_max_gradient" => 100.0,
    "nlp_scaling_min_value" => 1e-8,
    # initialization
    "bound_push" => 1e-8,
    "bound_frac" => 1e-6,
    "slack_bound_push" => 1e-8,
    "slack_bound_frac" => 1e-6,
    # linear solver
    "linear_solver" => "mumps",
    "linear_scaling_on_demand" => "yes",
    "fast_step_computation" => "no",
    "min_refinement_steps" => 1,
    "max_refinement_steps" => 8,
    "neg_curv_test_tol" => 1e-12,
    "neg_curv_test_reg" => "yes",
    "max_hessian_perturbation" => 1e20,
    "min_hessian_perturbation" => 1e-20,
    "mumps_mem_percent" => 500,
    "mumps_pivtol" => 1e-6,
    "mumps_pivtolmax" => 0.1,
    "mumps_pivot_order" => 7,
    "hessian_approximation" => "exact",
    "limited_memory_max_history" => 6,
)

# ------------------------------------------------------------
# DC — OSQP
# ------------------------------------------------------------

const _SILENT_FACTORY_DC = optimizer_with_attributes(
    OSQP.Optimizer,
    "eps_abs" => 1e-6,
    "eps_rel" => 1e-6,
    "max_iter" => 50_000,
    "verbose" => false,
    "adaptive_rho" => true,
    "linsys_solver" => "qdldl",
)

# ------------------------------------------------------------
# SOCWR — Clarabel
# ------------------------------------------------------------

const _SILENT_FACTORY_SOCWR = optimizer_with_attributes(
    Clarabel.Optimizer,
    "verbose" => false,
    # convergence
    "tol_feas" => 1e-4,
    "tol_gap_abs" => 1e-4,
    "tol_gap_rel" => 1e-4,
    "max_iter" => 1200,
    # wall-clock cap per zone solve (native Clarabel cannot be interrupted by the driver)
    "time_limit" => 60.0,
    # dynamic regularization
    "dynamic_regularization_enable" => true,
    "dynamic_regularization_eps" => 1e-10,
    "dynamic_regularization_delta" => 1e-4,
    # static regularization
    "static_regularization_enable" => true,
    # equilibration
    "equilibrate_enable" => true,
    "equilibrate_min_scaling" => 1e-4,
    "equilibrate_max_scaling" => 1e4,
    # linear solver
    "direct_solve_method" => :qdldl,
    # iterative refinement
    "iterative_refinement_reltol" => 1e-9,
    "iterative_refinement_abstol" => 1e-9,
    # step control
    "max_step_fraction" => 0.7,
    "linesearch_backtrack_step" => 0.5,
    # termination
    "tol_ktratio" => 1e-5,
    "min_terminate_step_length" => 1e-6,
    "min_switch_step_length" => 5e-2,
    # reduced tolerances (fallback acceptance)
    "reduced_tol_gap_abs" => 5e-4,
    "reduced_tol_gap_rel" => 5e-4,
    "reduced_tol_feas" => 5e-4,
    "reduced_tol_infeas_abs" => 5e-11,
    "reduced_tol_infeas_rel" => 5e-4,
    "reduced_tol_ktratio" => 1e-4,
)

# ------------------------------------------------------------
# SOCWRNL — Ipopt + MA57 (smooth-NLP SOCWR), tight duals
# ------------------------------------------------------------

# Opt-in zone scaling (LIRION_SOCWRNL_ZONE_SCALING): MA57 automatic scaling
# plus acceptable_* = 1e-9.
const _SOCWRNL_ZONE_SCALING_ENV = "LIRION_SOCWRNL_ZONE_SCALING"
@inline function _socwrnl_zone_scaling_enabled()::Bool
    v = lowercase(strip(get(ENV, _SOCWRNL_ZONE_SCALING_ENV, "")))
    return v == "1" || v == "true" || v == "yes" || v == "on"
end
# Opt-in zone robustness (LIRION_SOCWRNL_ZONE_ROBUST ∈ {pivtol, muadapt, both}):
#   pivtol  : "ma57_pivtol" => 1e-4, "ma57_pivtolmax" => 0.5
#   muadapt : "mu_strategy" => "adaptive"
#   both    : all three
const _SOCWRNL_ZONE_ROBUST_ENV = "LIRION_SOCWRNL_ZONE_ROBUST"
function _socwrnl_zone_robust_mode()::Symbol
    v = lowercase(strip(get(ENV, _SOCWRNL_ZONE_ROBUST_ENV, "")))
    isempty(v) && return :off
    v == "pivtol" && return :pivtol
    v == "muadapt" && return :muadapt
    v == "both" && return :both
    error("$(_SOCWRNL_ZONE_ROBUST_ENV)=$(repr(v)) is not one of pivtol | muadapt | both (or unset).")
end
const _SOCWRNL_ROBUST_PIVTOL = 1e-4
const _SOCWRNL_ROBUST_PIVTOLMAX = 0.5

const _SOCWRNL_ACCEPTABLE_KEYS = (
    "acceptable_tol", "acceptable_dual_inf_tol",
    "acceptable_constr_viol_tol", "acceptable_compl_inf_tol",
)

# Opt-in acceptable-point fallback:
#   LIRION_SOCWRNL_ZONE_ACCEPT     = <n>  acceptable_iter (>=1 enables)
#   LIRION_SOCWRNL_ZONE_ACCEPT_TOL = <x>  acceptable_* value (default 1e-8)
const _SOCWRNL_ZONE_ACCEPT_ENV = "LIRION_SOCWRNL_ZONE_ACCEPT"
const _SOCWRNL_ZONE_ACCEPT_TOL_ENV = "LIRION_SOCWRNL_ZONE_ACCEPT_TOL"
const _SOCWRNL_ACCEPT_TOL_DEFAULT = 1e-8
# Returns (acceptable_iter::Int, acceptable_tol::Float64); (0, NaN) when OFF.
function _socwrnl_zone_accept_spec()
    v = lowercase(strip(get(ENV, _SOCWRNL_ZONE_ACCEPT_ENV, "")))
    (isempty(v) || v == "0" || v == "off" || v == "false" || v == "no") &&
        return (0, NaN)
    n = tryparse(Int, v)
    (n === nothing || n < 1) && error(
        "$(_SOCWRNL_ZONE_ACCEPT_ENV)=$(repr(v)) must be a positive integer " *
        "(acceptable_iter) or unset/0/off.")
    tv = strip(get(ENV, _SOCWRNL_ZONE_ACCEPT_TOL_ENV, ""))
    tol = isempty(tv) ? _SOCWRNL_ACCEPT_TOL_DEFAULT : tryparse(Float64, tv)
    (tol === nothing || !(tol > 0)) && error(
        "$(_SOCWRNL_ZONE_ACCEPT_TOL_ENV)=$(repr(tv)) must be a positive float.")
    return (n, tol)
end

function _socwrnl_tight_factory()
    hsl = get(ENV, "LIRION_HSLLIB", "")
    use_ma57 = !isempty(hsl) && isfile(hsl)
    if !use_ma57 && get(ENV, "LIRION_STRICT_SOLVERS", "0") == "1"
        error(
            "STRICT SOLVERS: SOCWRNL requires Ipopt+MA57 but LIRION_HSLLIB is " *
            "unset/invalid ($(isempty(hsl) ? "unset" : hsl)); refusing to " *
            "silently fall back to MUMPS for a measurability run.",
        )
    end
    use_ma57 ||
        @warn "SOCWRNL: MA57 unavailable (LIRION_HSLLIB unset/invalid) — \
               falling back to MUMPS; duals may be less precise."

    drop = Set([
        "tol", "dual_inf_tol", "constr_viol_tol", "compl_inf_tol",
        "acceptable_iter", "bound_relax_factor", "linear_solver", "hsllib",
    ])
    zone_scaling = _socwrnl_zone_scaling_enabled()
    accept_iter, accept_tol = _socwrnl_zone_accept_spec()
    # either flag overrides the inherited acceptable_* (the tighter wins if both are set)
    override_acceptable = zone_scaling || accept_iter > 0
    override_acceptable && union!(drop, _SOCWRNL_ACCEPTABLE_KEYS)
    keep = Pair[]
    for p in _SILENT_FACTORY_AC.params
        nm = p.first isa MOI.RawOptimizerAttribute ? p.first.name : ""
        nm in drop && continue
        push!(keep, p)
    end

    # Any-valued pairs: acceptable_iter must stay an Int
    extra = Pair{MOI.RawOptimizerAttribute, Any}[
        MOI.RawOptimizerAttribute("tol") => 1e-9,
        MOI.RawOptimizerAttribute("dual_inf_tol") => 1e-9,
        MOI.RawOptimizerAttribute("constr_viol_tol") => 1e-9,
        MOI.RawOptimizerAttribute("compl_inf_tol") => 1e-9,
        MOI.RawOptimizerAttribute("bound_relax_factor") => 0.0,
        MOI.RawOptimizerAttribute("acceptable_iter") => (accept_iter > 0 ? accept_iter : 0),
    ]
    if override_acceptable
        # scaling intends 1e-9; accept intends accept_tol; both ⇒ the tighter.
        acc_tol = zone_scaling ? (accept_iter > 0 ? min(1e-9, accept_tol) : 1e-9) :
                                 accept_tol
        for nm in _SOCWRNL_ACCEPTABLE_KEYS
            push!(extra, MOI.RawOptimizerAttribute(nm) => acc_tol)
        end
    end
    if zone_scaling
        if use_ma57
            push!(extra, MOI.RawOptimizerAttribute("ma57_automatic_scaling") => "yes")
        elseif get(ENV, "LIRION_STRICT_SOLVERS", "0") == "1"
            error(
                "STRICT SOLVERS: $(_SOCWRNL_ZONE_SCALING_ENV) requires Ipopt+MA57 " *
                "(ma57_automatic_scaling is an MA57 option) but LIRION_HSLLIB is " *
                "unset/invalid ($(isempty(hsl) ? "unset" : hsl)).",
            )
        else
            @warn "SOCWRNL: $(_SOCWRNL_ZONE_SCALING_ENV) requested but MA57 is " *
                  "unavailable — applying only the acceptable_* = 1e-9 guard (MUMPS)."
        end
    end
    robust = _socwrnl_zone_robust_mode()
    if robust === :pivtol || robust === :both
        if use_ma57
            push!(extra, MOI.RawOptimizerAttribute("ma57_pivtol") => _SOCWRNL_ROBUST_PIVTOL)
            push!(extra, MOI.RawOptimizerAttribute("ma57_pivtolmax") => _SOCWRNL_ROBUST_PIVTOLMAX)
        elseif get(ENV, "LIRION_STRICT_SOLVERS", "0") == "1"
            error(
                "STRICT SOLVERS: $(_SOCWRNL_ZONE_ROBUST_ENV)=$(robust) requires Ipopt+MA57 " *
                "(ma57_pivtol is an MA57 option) but LIRION_HSLLIB is unset/invalid.",
            )
        else
            @warn "SOCWRNL: $(_SOCWRNL_ZONE_ROBUST_ENV)=$(robust) requested but MA57 is " *
                  "unavailable — ma57_pivtol/ma57_pivtolmax not applied (MUMPS)."
        end
    end
    if robust === :muadapt || robust === :both
        push!(extra, MOI.RawOptimizerAttribute("mu_strategy") => "adaptive")
    end
    if use_ma57
        push!(extra, MOI.RawOptimizerAttribute("linear_solver") => "ma57")
        push!(extra, MOI.RawOptimizerAttribute("hsllib") => hsl)
    else
        push!(extra, MOI.RawOptimizerAttribute("linear_solver") => "mumps")
    end

    return optimizer_with_attributes(Ipopt.Optimizer, vcat(keep, extra)...)
end

# ============================================================
# User-supplied optimizer
# ============================================================

function _silent_factory(opt::MOI.OptimizerWithAttributes)
    params = filter(p -> !(p.first isa MOI.Silent), opt.params)
    return optimizer_with_attributes(
        opt.optimizer_constructor,
        params...,
        MOI.Silent() => true,
    )
end

function _silent_factory(opt::Type{<:MOI.AbstractOptimizer})
    @info "Custom optimizer $opt passed; only MOI.Silent() applied. \
           Lirion's default tuning was skipped."
    return optimizer_with_attributes(opt, MOI.Silent() => true)
end

function _silent_factory(opt)
    @warn "Unknown optimizer specification $opt; output may not be suppressed."
    return opt
end

# ============================================================
# Entry point
# ============================================================

resolve_optimizer(::Val{:ac}, user = nothing) =
    user === nothing ? _SILENT_FACTORY_AC : _silent_factory(user)
resolve_optimizer(::Val{:dc}, user = nothing) =
    user === nothing ? _SILENT_FACTORY_DC : _silent_factory(user)
resolve_optimizer(::Val{:socwr}, user = nothing) =
    user === nothing ? _SILENT_FACTORY_SOCWR : _silent_factory(user)
function resolve_optimizer(::Val{:socwrnl}, user = nothing)
    user === nothing && return _socwrnl_tight_factory()
    # a user-supplied optimizer bypasses the tight tuning: error under STRICT, warn otherwise
    msg =
        "SOCWRNL: a custom optimizer was supplied; Lirion's measurability tuning " *
        "(dual_inf_tol = 1e-9, bound_relax_factor = 0.0, MA57) will be DISCARDED. " *
        "The SOCWRNL certificate budget tol·‖κ‖ depends on dual_inf_tol = 1e-9; a " *
        "looser dual tolerance makes the Tier-2 measurability result invalid."
    if get(ENV, "LIRION_STRICT_SOLVERS", "0") == "1"
        error("STRICT SOLVERS: " * msg)
    end
    @warn msg
    return _silent_factory(user)
end

@inline function resolve_optimizer(
    config::AbstractDistributedConfig,
    user = nothing,
)
    return resolve_optimizer(base_model_tag(config.base_model), user)
end

# ============================================================
# Feasibility check
# ============================================================

# Counters for reduced-accuracy solves and non-feasible dual statuses (instrumentation only).
const _REDUCED_ACCURACY_EVENTS = Threads.Atomic{Int}(0)
const _DUAL_STATUS_EVENTS = Threads.Atomic{Int}(0)

function _require_optimal!(model::JuMP.Model; iter::Int, zone_id::Int)
    ts = JuMP.termination_status(model)
    ps = JuMP.primal_status(model)
    ds = JuMP.dual_status(model)

    acceptable =
        ts == MOI.OPTIMAL || ts == MOI.LOCALLY_SOLVED ||
        ts == MOI.ALMOST_OPTIMAL || ts == MOI.ALMOST_LOCALLY_SOLVED

    # reduced-accuracy statuses are accepted, but never silently
    if ts == MOI.ALMOST_OPTIMAL || ts == MOI.ALMOST_LOCALLY_SOLVED
        Threads.atomic_add!(_REDUCED_ACCURACY_EVENTS, 1)
        @warn "Local zone solve accepted with a REDUCED-accuracy status: iter=$iter " *
              "zone=$zone_id termination_status=$ts primal_status=$ps dual_status=$ds. " *
              "Primal/dual values may only satisfy the solver's acceptable/reduced " *
              "tolerances (not the nominal ones)." maxlog = 5
    end
    # ALADIN consumes the zone duals; report a non-feasible dual status
    if acceptable && !(ds == MOI.FEASIBLE_POINT || ds == MOI.NEARLY_FEASIBLE_POINT)
        Threads.atomic_add!(_DUAL_STATUS_EVENTS, 1)
        @warn "Local zone solve iter=$iter zone=$zone_id: termination_status=$ts but " *
              "dual_status=$ds — zone multipliers are not a (nearly) feasible dual " *
              "point; ALADIN active sets / Hessians / κ̄ built from them are unreliable." maxlog = 5
    end

    if !acceptable
        raw = try
            JuMP.raw_status(model)
        catch
            "raw_status unavailable"
        end
        error("""
        iter=$iter zone=$zone_id
          termination_status = $ts
          primal_status       = $ps
          dual_status         = $ds
          raw_status          = $raw
        """)
    end

    return nothing
end
