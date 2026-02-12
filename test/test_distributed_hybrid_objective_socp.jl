# ==============================================
# ADMM-Hybrid-SOCP vs SOCP (objective) – tests (API solve)
# ==============================================

using Test
using Lirion
using Glob
using SCS
using Gurobi
using COSMO
using Clarabel

# -------------------------------------------------
# (Opcional) pacote de atributos robustos do Clarabel
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
if !@isdefined(pick_socp_optimizer)
    const _SOCP_CHAIN_MAP_LOCAL = Dict(
        :auto   => (Gurobi.Optimizer, COSMO.Optimizer, SCS.Optimizer),
        :gurobi => (Gurobi.Optimizer, COSMO.Optimizer, SCS.Optimizer),
        :cosmo  => (COSMO.Optimizer,  Gurobi.Optimizer, SCS.Optimizer),
        :scs    => (SCS.Optimizer,    Gurobi.Optimizer, COSMO.Optimizer),
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

# -------------------------------------------------
# Cases
# -------------------------------------------------
const distributed_hybrid_socp_zone_cases = [
    # "pglib_opf_case73_ieee_rts__api.m",
    # "pglib_opf_case73_ieee_rts__sad.m",
    ## "pglib_opf_case57_ieee__api.m",
    # "pglib_opf_case57_ieee__sad.m",

    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
    "pglib_opf_case73_ieee_rts.m",
    # "pglib_opf_case57_ieee.m",
    # "pglib_opf_case118_ieee.m",
    # "pglib_opf_case200_activ.m",
    # "pglib_opf_case500_goc.m",
]

function track_failure_hybrid_socp(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] =
            Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

# ------------------------------------------------------------
# Objective extraction (same pattern as branch-AC)
# ------------------------------------------------------------
@inline function _objective_centralized_socp(out::Lirion.LirionRunResult)
    d = out.result
    d isa AbstractDict || error("Centralized out.result must be Dict, got $(typeof(d))")
    haskey(d, "objective") || error("Centralized result has no key 'objective'. Keys=$(collect(keys(d)))")
    return d["objective"]
end

@inline function _objective_distributed_hybrid_socp(out::Lirion.LirionRunResult)
    r = out.result
    r isa NamedTuple || error("Distributed out.result must be NamedTuple, got $(typeof(r))")
    # ADMM distributed results expose these (same as your branch tests)
    hasproperty(r, :total_expr_cost) || error("Distributed result has no field :total_expr_cost. Fields=$(propertynames(r))")
    return getproperty(r, :total_expr_cost)
end

function compare_distributed_hybrid_socp_vs_centralized(
    file::String;
    label::String,
    central_optimizer = _CENTRAL_DEFAULT,
    distributed_optimizer = _DISTRIB_DEFAULT,
    rtol = 0.01,
    atol = 1e-6,
    # ADMM defaults
    rho_v::Float64 = 1.0e3,
    rho_flow_p::Float64 = 3.0e2,
    rho_flow_q::Float64 = 3.0e2,
    admm_max_iter::Int = 350,
    admm_eps_abs::Float64 = 1e-4,
    admm_eps_rel::Float64 = 1e-4,
    admm_adaptive::Union{Nothing, Bool} = nothing,
    admm_mu_incr::Float64 = 3.0,
    admm_mu_decr::Float64 = 3.0,
    admm_tau_incr::Float64 = 1.05,
    admm_tau_decr::Float64 = 0.98,
)
    case_name = basename(file)

    # ===========================================================
    # Caso especial: 2000-bus / 5 zonas
    # ===========================================================
    if case_name == "pglib_opf_case2000_5zones_goc.m"
        @info "⚙️ [HYBRID-SOCP] Caso especial $case_name com Clarabel robusto (central) e Clarabel com rho=1e8 (distribuído)."

        cent_obj = try
            out_cent = Lirion.solve(
                file;
                model = Lirion.socp(),
                algorithm = Lirion.centralized(),
                optimizer = eval(Meta.parse(CLARABEL_ROBUST)),
                outputs = false,
                plots = false,
            )
            _objective_centralized_socp(out_cent)
        catch e
            @warn "Centralized SOCP failed on $case_name: $e"
            track_failure_hybrid_socp(label, "Centralized", case_name)
            return
        end

        dist_obj = try
            alg = Lirion.admm_hybrid(
                max_iter = 700,
                rho_v = 1e8,
                rho_flow_p = 1e8,
                rho_flow_q = 1e8,
                eps_abs = 1e-4,
                eps_rel = 1e-4,
                adaptive_rho = true,
                mu_incr = 3.0,
                mu_decr = 3.0,
                tau_incr = 1.1,
                tau_decr = 0.9,
            )

            out_dist = Lirion.solve(
                file;
                model = Lirion.socp(),
                algorithm = alg,
                optimizer = Clarabel.Optimizer,
                outputs = false,
                plots = false,
            )
            _objective_distributed_hybrid_socp(out_dist)
        catch e
            @warn "Distributed ADMM-Hybrid-SOCP failed on $case_name: $e"
            track_failure_hybrid_socp(label, "Distributed", case_name)
            return
        end

        @info "File $case_name | centralized SOCP = $(round(cent_obj, digits=6)) | distributed-hybrid-SOCP (total_expr_cost) = $(round(dist_obj, digits=6))"
        @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
        return
    end

    # ===========================================================
    # Casos gerais
    # ===========================================================
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
        track_failure_hybrid_socp(label, "Centralized", case_name)
        return
    end

    dist_obj = try
        # Hybrid-SOCP needs separate rho_v / rho_flow_p / rho_flow_q (do NOT use --admm-rho here)
        alg = Lirion.admm_hybrid(
            max_iter = admm_max_iter,
            rho_v = Float64(rho_v),
            rho_flow_p = Float64(rho_flow_p),
            rho_flow_q = Float64(rho_flow_q),
            eps_abs = Float64(admm_eps_abs),
            eps_rel = Float64(admm_eps_rel),
            adaptive_rho = isnothing(admm_adaptive) ? true : Bool(admm_adaptive),
            mu_incr = Float64(admm_mu_incr),
            mu_decr = Float64(admm_mu_decr),
            tau_incr = Float64(admm_tau_incr),
            tau_decr = Float64(admm_tau_decr),
        )

        out_dist = Lirion.solve(
            file;
            model = Lirion.socp(),
            algorithm = alg,
            optimizer = distributed_optimizer,
            outputs = false,
            plots = false,
        )
        _objective_distributed_hybrid_socp(out_dist)
    catch e
        @warn "Distributed ADMM-Hybrid-SOCP failed on $case_name: $e"
        track_failure_hybrid_socp(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized SOCP = $(round(cent_obj, digits=6)) | distributed-hybrid-SOCP (total_expr_cost) = $(round(dist_obj, digits=6)) | max_iter=$(admm_max_iter), rhoV=$(rho_v), rhoP=$(rho_flow_p), rhoQ=$(rho_flow_q)"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_hybrid_socp_objective_tests(;
    case_list = distributed_hybrid_socp_zone_cases,
    label = "ADMM-Hybrid-SOCP vs SOCP",
    central_optimizer = _CENTRAL_DEFAULT,
    distributed_optimizer = _DISTRIB_DEFAULT,
    rtol = 0.01,
    atol = 1e-6,
    admm_eps_abs = 1e-4,
    admm_eps_rel = 1e-4,
)
    if isempty(case_list)
        @info "Empty hybrid-SOCP distributed case list. No tests will be run."
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

        # defaults por case
        rho_v_case = 1.0e4
        rho_flow_p_case = 3.0e2
        rho_flow_q_case = 3.0e2
        admm_max_iter_case = 5000
        admm_eps_abs_case = admm_eps_abs
        admm_eps_rel_case = admm_eps_rel

        @info "[$label] Rodando caso $i/$n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rhoV=$(rho_v_case), rhoP=$(rho_flow_p_case), rhoQ=$(rho_flow_q_case), " *
              "eps_pr=$(admm_eps_abs_case), eps_dual=$(admm_eps_rel_case))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_hybrid_socp_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer,
                distributed_optimizer = distributed_optimizer,
                rtol = rtol,
                atol = atol,
                rho_v = rho_v_case,
                rho_flow_p = rho_flow_p_case,
                rho_flow_q = rho_flow_q_case,
                admm_max_iter = admm_max_iter_case,
                admm_eps_abs = admm_eps_abs_case,
                admm_eps_rel = admm_eps_rel_case,
            )
        end
    end
    return true
end