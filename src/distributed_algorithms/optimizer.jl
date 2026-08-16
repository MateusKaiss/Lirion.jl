# ============================================================
# Dispatch
# ============================================================

@inline base_model_tag(::Type{AC}) = Val(:ac)
@inline base_model_tag(::Type{DC}) = Val(:dc)
@inline base_model_tag(::Type{SOCWR}) = Val(:socwr)

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
    # wall-clock cap per zone solve: a pathological large-zone SOC solve (e.g. case588
    # k4 socwr/tie zones) can grind the full 1200 iters over tens of minutes; the
    # distributed driver's ALADIN attempt-timeout cannot interrupt native Clarabel, so
    # the whole run hangs. 60s lets a stuck solve return so the caller can proceed.
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

@inline function resolve_optimizer(
    config::AbstractDistributedConfig,
    user = nothing,
)
    return resolve_optimizer(base_model_tag(config.base_model), user)
end

# ============================================================
# Feasibility check
# ============================================================

function _require_optimal!(model::JuMP.Model; iter::Int, zone_id::Int)
    ts = JuMP.termination_status(model)
    ps = JuMP.primal_status(model)
    ds = JuMP.dual_status(model)

    acceptable =
        ts == MOI.OPTIMAL || ts == MOI.LOCALLY_SOLVED ||
        ts == MOI.ALMOST_OPTIMAL || ts == MOI.ALMOST_LOCALLY_SOLVED

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
