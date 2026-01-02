# ============================================================
# Optional-field helpers
# ============================================================

@inline _adaptive_on(config) = _getopt(config, :adaptive_rho, true)

@inline function _getopt(config, name::Symbol, default)
    if hasproperty(config, name)
        v = getproperty(config, name)
        return v === nothing ? default : v
    end
    return default
end

# ============================================================
# Optimizer resolution
# ============================================================

@inline function _silent_factory(opt)
    if opt === Ipopt.Optimizer
        return optimizer_with_attributes(
            Ipopt.Optimizer,
            "print_level" => 0,
            "sb" => "yes",
            "print_timing_statistics" => "no",
            "print_user_options" => "no",
            "file_print_level" => 0,
        )
    elseif opt === Gurobi.Optimizer
        return optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0)
    elseif opt === SCS.Optimizer
        return optimizer_with_attributes(SCS.Optimizer, "verbose" => 0)
    elseif opt === Clarabel.Optimizer
        return optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
    elseif isdefined(Main, :HiGHS) && opt === HiGHS.Optimizer
        return optimizer_with_attributes(
            HiGHS.Optimizer,
            "log_to_console" => false,
            "output_flag" => false,
        )
    elseif isdefined(Main, :GLPK) && opt === GLPK.Optimizer
        return optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => 0)
    elseif opt === COSMO.Optimizer
        return optimizer_with_attributes(COSMO.Optimizer, "verbose" => false)
    else
        return opt
    end
end

# ----------------------------
# DC
# ----------------------------
function resolve_optimizer(::Val{:dc}, config)
    opt = _getopt(config, :optimizer, nothing)
    if opt !== nothing
        if opt == COSMO.Optimizer
            return optimizer_with_attributes(
                COSMO.Optimizer,
                "eps_abs" => 1e-6,
                "eps_rel" => 1e-6,
                "max_iter" => 10_000,
                "verbose" => false,
                "decompose" => true,
            )
        end
        return _silent_factory(opt)
    end

    eps_abs = _getopt(config, :eps_primal, 1e-6)
    eps_rel = _getopt(config, :eps_dual, 1e-6)

    return optimizer_with_attributes(
        OSQP.Optimizer,
        "eps_abs" => eps_abs,
        "eps_rel" => eps_rel,
        "max_iter" => 50_000,
        "verbose" => false,
        "adaptive_rho" => true,
        "linsys_solver" => "qdldl",
    )
end

# ----------------------------
# AC
# ----------------------------
function resolve_optimizer(::Val{:ac}, config)
    opt = _getopt(config, :optimizer, nothing)
    if opt !== nothing
        if opt == Ipopt.Optimizer
            return optimizer_with_attributes(
                Ipopt.Optimizer,
                # saída silenciosa
                "print_level" => 0,
                "sb" => "yes",
                "print_timing_statistics" => "no",
                "print_user_options" => "no",
                "file_print_level" => 0,

                # critérios
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

                # init
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
        end
        return _silent_factory(opt)
    end

    return optimizer_with_attributes(
        Ipopt.Optimizer,
        # saída silenciosa
        "print_level" => 0,
        "sb" => "yes",
        "print_timing_statistics" => "no",
        "print_user_options" => "no",
        "file_print_level" => 0,

        # critérios
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

        # init
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
end

# ----------------------------
# SOCP
# ----------------------------
is_gurobi(m::JuMP.Model)::Bool = occursin("Gurobi", MOI.get(m, MOI.SolverName()))

function resolve_optimizer(::Val{:socp}, config)
    opt = _getopt(config, :optimizer, nothing)
    opt !== nothing && return _silent_factory(opt)

    prefer = _getopt(config, :prefer_socp_solver, :auto)

    order =
        prefer === :gurobi ? (:gurobi, :clarabel, :scs) :
        prefer === :clarabel ? (:clarabel, :gurobi, :scs) :
        prefer === :scs ? (:scs, :gurobi, :clarabel) :
        (:gurobi, :clarabel, :scs)

    for s in order
        try
            if s === :gurobi
                return _silent_factory(Gurobi.Optimizer)
            elseif s === :clarabel
                return _silent_factory(Clarabel.Optimizer)
            else
                return _silent_factory(SCS.Optimizer)
            end
        catch
        end
    end
    return _silent_factory(SCS.Optimizer)
end

@inline base_model_tag(::Type{DC}) = Val(:dc)
@inline base_model_tag(::Type{AC}) = Val(:ac)
@inline base_model_tag(::Type{SOCP}) = Val(:socp)

@inline function resolve_optimizer(config)
    return resolve_optimizer(base_model_tag(config.base_model), config)
end

# ============================================================
# Pretty time + misc
# ============================================================

@inline function _pretty_time(sec::Real)
    if sec < 60
        return "$(round(sec; digits = 3)) s"
    elseif sec < 3600
        return "$(round(sec / 60; digits = 2)) min"
    else
        return "$(round(sec / 3600; digits = 2)) h"
    end
end

@inline _wrap_angle_diff(Δ::Float64) = atan(sin(Δ), cos(Δ))

function compute_real_generation_cost(model, params)
    T = get(model.ext, :time_periods, 1)
    active_generators = findall(params.generator_status .== 1)

    total_cost = 0.0
    for t in 1:T, g in active_generators
        if haskey(params.generator_costs, t) && haskey(params.generator_costs[t], g)
            pgen = value(model[:p_generation][g, t]) * params.sbase
            coefs = params.generator_costs[t][g]
            degree = length(coefs)
            @inbounds for i in 1:degree
                total_cost += coefs[i] * pgen^(degree - i)
            end
        end
    end
    return total_cost
end

# ============================================================
# Finalization + summary
# ============================================================

@inline function _admm_totals(partitions)
    parts = partitions.parts
    total_real = 0.0
    total_expr = 0.0
    @inbounds @simd for k in eachindex(parts)
        p = parts[k]
        total_real += p.real_generation_cost
        total_expr += p.generation_cost_model_expr
    end
    return total_real, total_expr
end

@inline _elapsed_s(t0_ns::Integer) = (time_ns() - t0_ns) / 1e9

@inline function _finalize_admm_result(partitions, t0_ns::Integer, iter_done::Integer)
    total_real_cost, total_expr_cost = _admm_totals(partitions)
    elapsed_s = _elapsed_s(t0_ns)
    return (;
        total_real_cost,
        total_expr_cost,
        partitions,
        iter_done = Int(iter_done),
        elapsed_s,
    )
end

function _print_admm_summary!(
    tag::AbstractString,
    result;
    rho_lines = (),
)
    println("\n$tag Resumo:")
    println("→ Total (value(:generation_cost))     = $(result.total_expr_cost)")
    println("→ Total (reavaliado via p_generation) = $(result.total_real_cost)")
    for line in rho_lines
        println(line)
    end
    println("Iterações realizadas: $(result.iter_done)")
    println("Tempo ADMM: $(_pretty_time(result.elapsed_s))")
    return nothing
end

# ------------------------------------------------------------
# Diagnostics / history
# ------------------------------------------------------------

@inline function _record_iter!(
    hist::ADMMHistory,
    r::Float64,
    s::Float64,
    eps_r::Float64,
    eps_s::Float64;
    rhos::Dict{Symbol, <:Real} = Dict{Symbol, Real}(),
)
    push!(hist.r_primal, r)
    push!(hist.r_dual, s)
    push!(hist.eps_primal, eps_r)
    push!(hist.eps_dual, eps_s)
    for (k, v) in rhos
        _push_rho!(hist, k, v)
    end
    return nothing
end
