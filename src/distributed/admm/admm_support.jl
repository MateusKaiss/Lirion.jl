@inline _adaptive_on(config) =
    hasfield(typeof(config), :adaptive_rho) ? getfield(config, :adaptive_rho) : true

function _set_model_silent!(m::JuMP.Model)
    _safe(k, v) = (
        try
            set_optimizer_attribute(m, k, v)
        catch
        end
    )

    try
        set_silent(m)
    catch
    end
    try
        set_optimizer_attribute(m, JuMP.MOI.Silent(), true)
    catch
    end

    solver = ""
    try
        backend = JuMP.backend(m)
        solver = lowercase(JuMP.MOI.get(backend, JuMP.MOI.SolverName()))
    catch
    end

    if occursin("ipopt", solver)
        _safe("print_level", 0)
        _safe("sb", "yes")
    elseif occursin("gurobi", solver)
        _safe("OutputFlag", 0)
    elseif occursin("cplex", solver)
        _safe("CPXPARAM_ScreenOutput", 0)
    elseif occursin("highs", solver)
        _safe("log_to_console", false)
        _safe("output_flag", false)
    elseif occursin("scs", solver)
        _safe("verbose", false)
    elseif occursin("clarabel", solver)
        _safe("verbose", false)
    elseif occursin("glpk", solver)
        _safe("msg_lev", 0)
    else
        _safe("verbose", false)
        _safe("output_flag", 0)
        _safe("print_level", 0)
    end
    return nothing
end

function resolve_optimizer_dc(config)
    opt = try
        getfield(config, :optimizer)
    catch
        nothing
    end

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
        else
            return opt
        end
    end

    eps_abs = try
        getfield(config, Symbol("eps_primal"))
    catch
        1e-6
    end
    eps_rel = try
        getfield(config, Symbol("eps_dual"))
    catch
        1e-6
    end

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

###########################
# HELPERS GENÉRICOS AC    #
###########################
function resolve_optimizer_ac(config)
    if hasproperty(config, :optimizer) && getproperty(config, :optimizer) !== nothing
        return getproperty(config, :optimizer)
    else
        return optimizer_with_attributes(
            Ipopt.Optimizer,

            # ---------------------------
            # Saída silenciosa
            # ---------------------------
            "print_level" => 0,
            "sb" => "yes",
            "print_timing_statistics" => "no",
            "print_user_options" => "no",
            "file_print_level" => 0,

            # ---------------------------
            # Critérios de parada
            # ---------------------------
            "tol" => 1e-7,
            "acceptable_tol" => 1e-5,
            "acceptable_iter" => 10,
            "dual_inf_tol" => 1e-4,
            "constr_viol_tol" => 1e-4,
            "compl_inf_tol" => 1e-4,
            "acceptable_dual_inf_tol" => 1e-2,
            "acceptable_constr_viol_tol" => 1e-3,
            "acceptable_compl_inf_tol" => 1e-3,

            # ---------------------------
            # NLP & Bounds
            # ---------------------------
            "fixed_variable_treatment" => "make_parameter_nodual",
            "bound_relax_factor" => 1e-6,
            "check_derivatives_for_naninf" => "no",
            "jacobian_approximation" => "exact",
            "gradient_approximation" => "exact",
            "kappa_d" => 1e-5,

            # ---------------------------
            # Escalonamento
            # ---------------------------
            "nlp_scaling_max_gradient" => 100.0,
            "nlp_scaling_min_value" => 1e-8,

            # ---------------------------
            # Inicialização
            # ---------------------------
            "bound_push" => 1e-8,
            "bound_frac" => 1e-6,
            "slack_bound_push" => 1e-8,
            "slack_bound_frac" => 1e-6,

            # ---------------------------
            # Linear Solver & Refinement
            # ---------------------------
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
end


"Retorna true se o modelo `m` está usando Gurobi."
is_gurobi(m::JuMP.Model)::Bool = occursin("Gurobi", MOI.get(m, MOI.SolverName()))


# ---- Otimizador: tenta Gurobi -> Clarabel -> SCS (com override) ----
function resolve_socp_optimizer(config)
    # 1) Prioridade absoluta ao que veio no config (se existir e não-nothing)
    if Base.hasfield(typeof(config), :optimizer)
        opt = getfield(config, :optimizer)
        if opt !== nothing
            return opt              # ex.: Clarabel.Optimizer, Gurobi.Optimizer, SCS.Optimizer
        end
    end

    # 2) Ordem preferida opcional: :auto | :gurobi | :clarabel | :scs
    prefer = (
        Base.hasfield(typeof(config), :prefer_socp_solver) ?
        getfield(config, :prefer_socp_solver) : :auto
    )

    order =
        prefer === :gurobi ? (:gurobi, :clarabel, :scs) :
        prefer === :clarabel ? (:clarabel, :gurobi, :scs) :
        prefer === :scs ? (:scs, :gurobi, :clarabel) :
        (:gurobi, :clarabel, :scs)   # :auto

    # 3) Retorna o primeiro construtor disponível (factory), sem instanciar aqui
    for s in order
        try
            if s === :gurobi
                return Gurobi.Optimizer
            elseif s === :clarabel
                return Clarabel.Optimizer
            else
                return SCS.Optimizer
            end
        catch
            # pacote ausente? tenta o próximo
        end
    end
    return SCS.Optimizer
end

# (mantido, mas sem forçar nada além do que você já usa)
safe_set_optimizer_attr!(m, key, val) = (
    try
        set_optimizer_attribute(m, key, val)
    catch
        nothing
    end
)

"Configura tolerâncias/accel do SCS (mantido; use se quiser)."
function set_socp_defaults!(m; eps_abs = 1e-6, eps_rel = 1e-6, max_iters = 200_000,
    acceleration_lookback = 10, acceleration_interval = 10,
    alpha = 1.8, normalize = true, scale = 1.0)
    safe_set_optimizer_attr!(m, "eps_abs", eps_abs)
    safe_set_optimizer_attr!(m, "eps_rel", eps_rel)
    safe_set_optimizer_attr!(m, "max_iters", max_iters)
    safe_set_optimizer_attr!(m, "acceleration_lookback", acceleration_lookback)
    safe_set_optimizer_attr!(m, "acceleration_interval", acceleration_interval)
    safe_set_optimizer_attr!(m, "alpha", alpha)
    safe_set_optimizer_attr!(m, "normalize", normalize)
    safe_set_optimizer_attr!(m, "scale", scale)
    return nothing
end

@inline function _pretty_time(sec::Real)
    if sec < 60
        return "$(round(sec; digits = 3)) s"
    elseif sec < 3600
        return "$(round(sec/60; digits = 2)) min"
    else
        return "$(round(sec/3600; digits = 2)) h"
    end
end

# ---------- rho adaptativo (paridade com DC) ----------
function adaptive_rho!(rho::Base.RefValue{Float64}, r::Float64, s::Float64;
    μ_incr = 2.0, μ_decr = 2.0, τ_incr = 1.05, τ_decr = 0.98)
    if r > μ_incr * s
        rho[] *= τ_incr
    elseif s > μ_decr * r
        rho[] *= τ_decr
    end
end
##########################
# rho ADAPTATIVO           #
##########################

#############################
# rho ADAPTATIVO AC           #
#############################

function adaptive_rho_ac!(
    rho::Base.RefValue{Float64},
    r::Float64,
    s::Float64;
    μ_incr = 2.0,
    μ_decr = 2.0,
    τ_incr = 1.05,
    τ_decr = 0.98,
)
    if r > μ_incr * s
        rho[] *= τ_incr
    elseif s > μ_decr * r
        rho[] *= τ_decr
    end
end

function adaptive_rho_dc!(
    rho::Base.RefValue{Float64},
    r::Float64,
    s::Float64;
    μ_incr = 3.0,
    μ_decr = 3.0,
    τ_incr = 1.05,
    τ_decr = 0.98,
)
    if r > μ_incr * s
        rho[] *= τ_incr
    elseif s > μ_decr * r
        rho[] *= τ_decr
    end
end

