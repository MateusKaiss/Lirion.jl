@inline base_model_tag(::Type{DC}) = Val(:dc)
@inline base_model_tag(::Type{AC}) = Val(:ac)
@inline base_model_tag(::Type{ACR}) = Val(:acr)
@inline base_model_tag(::Type{IVR}) = Val(:ivr)
@inline base_model_tag(::Type{LPAC}) = Val(:lpac)
@inline base_model_tag(::Type{SOCWR}) = Val(:socwr)
@inline base_model_tag(::Type{SOCBF}) = Val(:socbf)

# ============================================================
# Silent optimizer factories — built once at module load
# Core solvers (always available): Ipopt, OSQP, Clarabel
# Optional solvers (extensions):   Gurobi → GurobiExt.jl
#                                   SCS    → SCSExt.jl
#                                   HiGHS  → HiGHSExt.jl
#                                   GLPK   → GLPKExt.jl
# ============================================================

const _IPOPT_SILENT_FACTORY = optimizer_with_attributes(
    Ipopt.Optimizer,
    "print_level" => 0,
    "sb" => "yes",
    "print_timing_statistics" => "no",
    "print_user_options" => "no",
    "file_print_level" => 0,
)
const _CLARABEL_SILENT_FACTORY = optimizer_with_attributes(
    Clarabel.Optimizer,
    "verbose" => false,
    # convergence
    "tol_feas" => 1e-4,
    "tol_gap_abs" => 1e-4,
    "tol_gap_rel" => 1e-4,
    "max_iter" => 1200,
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

_silent_factory(::Type{Ipopt.Optimizer}) = _IPOPT_SILENT_FACTORY
_silent_factory(::Type{Clarabel.Optimizer}) = _CLARABEL_SILENT_FACTORY

function _silent_factory(opt)
    @warn "No silent configuration registered for optimizer $opt; output may not be suppressed."
    return opt
end

# ============================================================
# DC
# ============================================================

function resolve_optimizer(::Val{:dc}, inputs_optimizer = nothing)
    inputs_optimizer !== nothing && return _silent_factory(inputs_optimizer)

    return optimizer_with_attributes(
        OSQP.Optimizer,
        "eps_abs" => 1e-6,
        "eps_rel" => 1e-6,
        "max_iter" => 50_000,
        "verbose" => false,
        "adaptive_rho" => true,
        "linsys_solver" => "qdldl",
    )
end

# ============================================================
# AC — full Ipopt config built once at module load
# ============================================================

const _IPOPT_AC_FACTORY = optimizer_with_attributes(
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

function resolve_optimizer(::Val{:ac}, inputs_optimizer = nothing)
    if inputs_optimizer !== nothing && inputs_optimizer !== Ipopt.Optimizer
        return _silent_factory(inputs_optimizer)
    end
    return _IPOPT_AC_FACTORY
end

# ============================================================
# ACR — shares Ipopt AC factory (same non-convex NLP, different form)
# ============================================================

function resolve_optimizer(::Val{:acr}, inputs_optimizer = nothing)
    if inputs_optimizer !== nothing && inputs_optimizer !== Ipopt.Optimizer
        return _silent_factory(inputs_optimizer)
    end
    return _IPOPT_AC_FACTORY
end

function resolve_optimizer(::Val{:ivr}, inputs_optimizer = nothing)
    if inputs_optimizer !== nothing && inputs_optimizer !== Ipopt.Optimizer
        return _silent_factory(inputs_optimizer)
    end
    return _IPOPT_AC_FACTORY
end

function resolve_optimizer(::Val{:lpac}, inputs_optimizer = nothing)
    if inputs_optimizer !== nothing && inputs_optimizer !== Ipopt.Optimizer
        return _silent_factory(inputs_optimizer)
    end
    return _IPOPT_AC_FACTORY
end

# ============================================================
# SOCWR
# ============================================================

is_gurobi(m::JuMP.Model)::Bool = occursin("Gurobi", MOI.get(m, MOI.SolverName()))

function _gurobi_socwr_factory()
    return error("Gurobi not loaded. Add `using Gurobi` to use Gurobi for SOCWR.")
end
function _scs_socwr_factory()
    return error("SCS not loaded. Add `using SCS` to use SCS for SOCWR.")
end

function resolve_optimizer(::Val{:socwr}, inputs_optimizer = nothing)
    inputs_optimizer !== nothing && return _silent_factory(inputs_optimizer)

    for s in (:clarabel, :gurobi, :scs)
        s === :clarabel && return _CLARABEL_SILENT_FACTORY
        s === :gurobi && try
            return _gurobi_socwr_factory()
        catch e
            @warn "Gurobi unavailable for SOCWR ($(sprint(showerror, e))); trying next solver."
            continue
        end
        s === :scs && try
            return _scs_socwr_factory()
        catch e
            @warn "SCS unavailable for SOCWR ($(sprint(showerror, e))); trying next solver."
            continue
        end
    end

    return _CLARABEL_SILENT_FACTORY
end

# ============================================================
# SOCBF — same conic solver as SOCWR
# ============================================================

function resolve_optimizer(::Val{:socbf}, inputs_optimizer = nothing)
    inputs_optimizer !== nothing && return _silent_factory(inputs_optimizer)
    return _CLARABEL_SILENT_FACTORY
end

# ============================================================
# Entry point
# ============================================================

@inline function resolve_optimizer(config::DistributedOPFModel, inputs_optimizer = nothing)
    return resolve_optimizer(base_model_tag(config.base_model), inputs_optimizer)
end

# ============================================================
# Feasibility check
# ============================================================

function _require_optimal!(model::JuMP.Model; iter::Int, zone_id::Int)
    ts = JuMP.termination_status(model)
    ps = JuMP.primal_status(model)
    ds = JuMP.dual_status(model)

    acceptable = ts == MOI.OPTIMAL || ts == MOI.LOCALLY_SOLVED

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
