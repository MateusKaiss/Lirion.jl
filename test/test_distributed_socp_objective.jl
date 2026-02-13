# ============================================
# ADMM-SOCP vs SOCP (objective) — tests (API solve)
# Corrigido p/ NOVA API: admm(coupling=:state, sectioning=:tie_line, ...)
# ============================================

const ADMM_MAX_IT = 5000

# -------------------------------------------------
# Clarabel robust preset (string kept for compatibility)
# -------------------------------------------------
const CLARABEL_ROBUST = """
optimizer_with_attributes(Clarabel.Optimizer,
    "tol_feas" => 1e-4,
    "tol_gap_abs" => 1e-4,
    "tol_gap_rel" => 1e-4,
    "max_iter" => 1200,

    "dynamic_regularization_enable" => true,
    "dynamic_regularization_eps"   => 1e-10,
    "dynamic_regularization_delta" => 1e-4,

    "static_regularization_enable" => true,

    "equilibrate_enable" => true,
    "equilibrate_min_scaling" => 1e-4,
    "equilibrate_max_scaling" => 1e4,

    "direct_solve_method" => :qdldl,

    "iterative_refinement_reltol" => 1e-9,
    "iterative_refinement_abstol" => 1e-9,

    "max_step_fraction" => 0.7,
    "linesearch_backtrack_step" => 0.5,

    "tol_ktratio" => 1e-5,
    "min_terminate_step_length" => 1e-6,
    "min_switch_step_length" => 5e-2,

    "reduced_tol_gap_abs" => 5e-4,
    "reduced_tol_gap_rel" => 5e-4,
    "reduced_tol_feas"    => 5e-4,
    "reduced_tol_infeas_abs" => 5e-11,
    "reduced_tol_infeas_rel" => 5e-4,
    "reduced_tol_ktratio"    => 1e-4
)"""

# -------- helper local (caso o runtests não tenha definido) --------
# ENV["LIRION_SOCP_PREFER"] pode ser "auto" | "gurobi" | "cosmo" | "scs"
if !@isdefined(pick_socp_optimizer)
    const _SOCP_CHAIN_MAP_LOCAL = Dict(
        :auto => (Gurobi.Optimizer, COSMO.Optimizer, SCS.Optimizer),
        :gurobi => (Gurobi.Optimizer, COSMO.Optimizer, SCS.Optimizer),
        :cosmo => (COSMO.Optimizer, Gurobi.Optimizer, SCS.Optimizer),
        :scs => (SCS.Optimizer, Gurobi.Optimizer, COSMO.Optimizer),
    )
    function pick_socp_optimizer(;
        prefer::Symbol = Symbol(lowercase(get(ENV, "LIRION_SOCP_PREFER", "auto"))),
    )
        chain = get(_SOCP_CHAIN_MAP_LOCAL, prefer, _SOCP_CHAIN_MAP_LOCAL[:auto])
        for factory in chain
            try
                _ = factory()
                return factory
            catch
            end
        end
        return SCS.Optimizer
    end
end

# Defaults: usa os do runtests se existirem; senão, cai para o picker local
const _CENTRAL_DEFAULT =
    @isdefined(CENTRAL_SOCP_SOLVER) ? CENTRAL_SOCP_SOLVER :
    (
        @isdefined(DEFAULT_SOCP_SOLVER) ? DEFAULT_SOCP_SOLVER :
        pick_socp_optimizer(; prefer = :gurobi)
    )

const _DISTRIB_DEFAULT =
    @isdefined(DISTRIBUTED_SOCP_SOLVER) ? DISTRIBUTED_SOCP_SOLVER :
    (
        @isdefined(DEFAULT_SOCP_SOLVER) ? DEFAULT_SOCP_SOLVER :
        pick_socp_optimizer(; prefer = :auto)
    )

# --- lista de casos particionados (2/3/4 zonas) ---
const distributed_socp_zone_cases = [
    # "pglib_opf_case73_ieee_rts__api.m",
    # "pglib_opf_case73_ieee_rts__sad.m",
    ## "pglib_opf_case57_ieee__api.m",
    # "pglib_opf_case57_ieee__sad.m",

    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
    "pglib_opf_case73_ieee_rts.m",
    "pglib_opf_case200_activ.m",

    # "pglib_opf_case118_ieee.m",
    # "pglib_opf_case57_ieee.m",
    # "pglib_opf_case500_goc.m",
]

function track_failure_socp(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] =
            Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

# ------------------------------------------------------------
# Objective extraction (current return formats)
# ------------------------------------------------------------
@inline function _objective_centralized_socp(out::Lirion.LirionRunResult)
    d = out.result
    d isa AbstractDict || error("Centralized out.result must be Dict, got $(typeof(d))")
    haskey(d, "objective") ||
        error("Centralized result has no key 'objective'. Keys=$(collect(keys(d)))")
    return Float64(d["objective"])
end

@inline function _objective_distributed_socp(out::Lirion.LirionRunResult)
    r = out.result
    r isa NamedTuple || error("Distributed out.result must be NamedTuple, got $(typeof(r))")
    hasproperty(r, :total_expr_cost) ||
        error(
            "Distributed result has no field :total_expr_cost. Fields=$(propertynames(r))",
        )
    return Float64(getproperty(r, :total_expr_cost))
end

# ------------------------------------------------------------
# Helper: call solve with/without optimizer kwarg
# ------------------------------------------------------------
@inline function _solve_maybe_optimizer(
    file::String;
    model,
    algorithm,
    optimizer,
    outputs = false,
    plots = false,
)
    if optimizer === nothing
        return Lirion.solve(
            file;
            model = model,
            algorithm = algorithm,
            outputs = outputs,
            plots = plots,
        )
    else
        return Lirion.solve(
            file;
            model = model,
            algorithm = algorithm,
            optimizer = optimizer,
            outputs = outputs,
            plots = plots,
        )
    end
end

# ------------------------------------------------------------
# Small helper: parse "on/off/nothing" into Bool (default=true)
# ------------------------------------------------------------
@inline function _parse_adaptive(admm_adaptive::Union{Nothing, String, Bool})
    admm_adaptive === nothing && return true
    admm_adaptive isa Bool && return admm_adaptive
    s = lowercase(String(admm_adaptive))
    if s == "on"
        return true
    elseif s == "off"
        return false
    else
        error(
            "admm_adaptive must be \"on\", \"off\", Bool, or nothing; got: $admm_adaptive",
        )
    end
end

function compare_distributed_socp_vs_centralized(
    file::String;
    label::String,
    central_optimizer = _CENTRAL_DEFAULT,
    distributed_optimizer = _DISTRIB_DEFAULT,
    rtol = 0.05,
    atol = 1e-6,

    # (mantidos por compat / log) — na NOVA API vamos usar um rho “base”
    rho_v::Float64 = 1.0e6,
    rho_p::Float64 = 1.0e6,
    rho_q::Float64 = 1.0e6, admm_max_iter::Int = ADMM_MAX_IT,
    admm_eps_abs::Float64 = 1e-4,
    admm_eps_rel::Float64 = 1e-4,
    admm_adaptive::Union{Nothing, Bool, String} = true,

    # MATCH do seu run “bom”:
    admm_mu_incr::Float64 = 10.0,
    admm_mu_decr::Float64 = 10.0,
    admm_tau_incr::Float64 = 1.05,
    admm_tau_decr::Float64 = 0.98,
)
    case_name = basename(file)

    # -----------------------------------------------------------
    # Centralized SOCP
    # -----------------------------------------------------------
    cent_obj = try
        out_cent = Lirion.solve(
            file;
            model = Lirion.socp(),
            algorithm = Lirion.centralized(),
            optimizer = central_optimizer,
            outputs = false,
            plots = false,
        )
        _objective_centralized_socp(out_cent)
    catch e
        @warn "Centralized SOCP failed on $case_name: $e"
        track_failure_socp(label, "Centralized", case_name)
        return
    end

    # -----------------------------------------------------------
    # Distributed ADMM-SOCP (NOVA API)
    # -----------------------------------------------------------
    dist_obj = try
        # Para SOCP (bus/state), usamos rho “base”.
        # Mantemos a intenção antiga (rho_v/p/q) como “sugestão” pro rho:
        rho_base = max(rho_v, rho_p, rho_q)

        alg = Lirion.admm(;
            coupling = :state,          # substitui admm_bus
            sectioning = :tie_line,
            max_iter = Int(admm_max_iter),
            rho = Float64(rho_base),
            eps_abs = Float64(admm_eps_abs),
            eps_rel = Float64(admm_eps_rel),
            adaptive_rho = _parse_adaptive(admm_adaptive),
            mu_incr = Float64(admm_mu_incr),
            mu_decr = Float64(admm_mu_decr),
            tau_incr = Float64(admm_tau_incr),
            tau_decr = Float64(admm_tau_decr),
        )

        out_dist = _solve_maybe_optimizer(
            file;
            model = Lirion.socp(),
            algorithm = alg,
            optimizer = distributed_optimizer,
            outputs = false,
            plots = false,
        )

        _objective_distributed_socp(out_dist)
    catch e
        @warn "Distributed ADMM-SOCP failed on $case_name: $e"
        track_failure_socp(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized SOCP = $(round(cent_obj, digits=6)) | distributed-ADMM-SOCP (total_expr_cost) = $(round(dist_obj, digits=6)) | max_iter=$(admm_max_iter) | rho=$(max(rho_v,rho_p,rho_q)) | mu=$(admm_mu_incr) | tau_incr=$(admm_tau_incr) | tau_decr=$(admm_tau_decr)"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_socp_objective_tests(;
    case_list = distributed_socp_zone_cases,
    label = "ADMM-SOCP vs SOCP",
    central_optimizer = _CENTRAL_DEFAULT,
    distributed_optimizer = nothing,
    rtol = 0.05,
    atol = 1e-6,
    admm_eps_abs = 1e-4,
    admm_eps_rel = 1e-4,
)
    if isempty(case_list)
        @info "Empty SOCP distributed case list. No tests will be run."
        return true
    end

    test_dir = joinpath(@__DIR__, "distributed")
    selected_files = [
        joinpath(test_dir, fname)
        for fname in case_list
        if isfile(joinpath(test_dir, fname))
    ]

    n = length(selected_files)
    for (i, file) in enumerate(selected_files)
        case_name = basename(file)

        rho_v_case = 1.0e6
        rho_p_case = 1.0e6
        rho_q_case = 1.0e6

        admm_max_iter_case = 10_000
        admm_eps_abs_case = admm_eps_abs
        admm_eps_rel_case = admm_eps_rel

        @info "[$label] Rodando caso $i/$n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rhoV=$(rho_v_case), rhoP=$(rho_p_case), rhoQ=$(rho_q_case), " *
              "eps_abs=$(admm_eps_abs_case), eps_rel=$(admm_eps_rel_case), mu=10, tau_incr=1.05, tau_decr=0.98, optimizer=$(distributed_optimizer === nothing ? "nothing(default)" : string(distributed_optimizer)))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_socp_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer,
                distributed_optimizer = distributed_optimizer,
                rtol = rtol,
                atol = atol,
                rho_v = rho_v_case,
                rho_p = rho_p_case,
                rho_q = rho_q_case,
                admm_max_iter = admm_max_iter_case,
                admm_eps_abs = admm_eps_abs_case,
                admm_eps_rel = admm_eps_rel_case,
                admm_adaptive = true,
                admm_mu_incr = 10.0,
                admm_mu_decr = 10.0,
                admm_tau_incr = 1.05,
                admm_tau_decr = 0.98,
            )
        end
    end

    return true
end
