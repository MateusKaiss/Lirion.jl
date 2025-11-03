# ==============================================
# ADMM-Hybrid-SOCP vs SOCP (objetivo) – testes
# ==============================================
using Test
using Lirion
using Glob
using SCS
using Gurobi
using COSMO
using Clarabel

# -------------------------------------------------
# Pacote de atributos robustos do Clarabel
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

# Defaults
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

const distributed_hybrid_socp_zone_cases = [
    # 14 buses
    "pglib_opf_case14_2zones_ieee_socp.m",
    "pglib_opf_case14_3zones_ieee_socp.m",
    # 24 buses
    "pglib_opf_case24_2zones_ieee_rts.m",
    "pglib_opf_case24_3zones_ieee_rts.m",
    # 30 buses
    "pglib_opf_case30_2zones_as.m",
    "pglib_opf_case30_3zones_as.m",
    # # 57 buses
    # # "pglib_opf_case57_2zones_ieee.m", problem
    # # "pglib_opf_case57_4zones_ieee.m", problem
    # # # 118 buses
    "pglib_opf_case118_2zones_ieee.m",
    "pglib_opf_case118_3zones_ieee.m",
    "pglib_opf_case118_5zones_ieee.m",

    # # 500 buses
    "pglib_opf_case500_3zones_goc.m",
    # "pglib_opf_case500_4zones_goc.m", problema
    # # 793
    # # "pglib_opf_case793_11zones_goc.m",
    # # 2000
    # "pglib_opf_case2000_5zones_goc.m",
    # # 3970
    # "pglib_opf_case3970_25zones_goc.m",
    # # 10000
    # "pglib_opf_case10000_50zones_goc.m",
    # # 20758
    # "pglib_opf_case20758_42zones_epigrids.m",
    # # 30000
    # "pglib_opf_case30000_20zones_goc.m",
    # # 78484
    # "pglib_opf_case78484_157zones_epigrids.m",
]
function track_failure_hybrid_socp(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] = Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

function compare_distributed_hybrid_socp_vs_centralized(
    file::String;
    label::String,
    central_optimizer = _CENTRAL_DEFAULT,
    distributed_optimizer = _DISTRIB_DEFAULT,
    rtol = 0.05,
    atol = 1e-6,
    rho_v::Float64 = 1.0e3,
    rho_flow_p::Float64 = 3.0e2,
    rho_flow_q::Float64 = 3.0e2,
    admm_max_iter::Int = 350,
    admm_eps_primal::Float64 = 1e-4,
    admm_eps_dual::Float64 = 1e-4,
)
    case_name = basename(file)

    # ===========================================================
    # Caso especial: 2000-bus / 5 zonas
    # ===========================================================
    if case_name == "pglib_opf_case2000_5zones_goc.m"
        @info "⚙️ [HYBRID-SOCP] Caso especial $case_name com Clarabel robusto (central) e Clarabel com ρ=1e8 (distribuído)."

        # --- Centralizado Clarabel robusto ---
        cent_obj = try
            args_cent = String[
                file,
                "--model", "socp",
                "--optimizer", CLARABEL_ROBUST,
            ]
            result_dict, _, _ = Lirion.main(args_cent)
            result_dict["objective"]
        catch e
            @warn "Centralized SOCP failed on $case_name: $e"
            track_failure_hybrid_socp(label, "Centralized", case_name)
            return
        end

        # --- Distribuído ADMM-Hybrid-SOCP com Clarabel e args específicos ---
        dist_obj = try
            args_dist = String[
                file,
                "--model", "distributed-admm-hybrid-socp",
                "--optimizer", "Clarabel.Optimizer",
                "--admm-max-iter", "700",
                "--admm-rho-v", "1e8",
                "--admm-rho-flow-p", "1e8",
                "--admm-rho-flow-q", "1e8",
                "--admm-eps-primal", "1e-4",
                "--admm-eps-dual", "1e-4",
                "--admm-adaptive", "on",
                "--admm-mu-incr", "3.0",
                "--admm-mu-decr", "3.0",
                "--admm-tau-incr", "1.1",
                "--admm-tau-decr", "0.9",
            ]
            obj, _ = Lirion.main(args_dist)
            obj
        catch e
            @warn "Distributed ADMM-Hybrid-SOCP failed on $case_name: $e"
            track_failure_hybrid_socp(label, "Distributed", case_name)
            return
        end

        @info "File $case_name | centralized SOCP = $(round(cent_obj, digits=6)) | distributed-hybrid-SOCP = $(round(dist_obj, digits=6))"
        @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
        return
    end

    # ===========================================================
    # Casos gerais (fallback)
    # ===========================================================
    cent_obj = try
        args_cent = String[
            file,
            "--model", "socp",
            "--optimizer", string(central_optimizer),
        ]
        result_dict, _, _ = Lirion.main(args_cent)
        result_dict["objective"]
    catch e
        @warn "Centralized SOCP failed on $case_name: $e"
        track_failure_hybrid_socp(label, "Centralized", case_name)
        return
    end

    dist_obj = try
        args_dist = String[
            file,
            "--model", "distributed-admm-hybrid-socp",
            "--optimizer", string(distributed_optimizer),
            "--admm-rho-v", string(rho_v),
            "--admm-rho-flow-p", string(rho_flow_p),
            "--admm-rho-flow-q", string(rho_flow_q),
            "--admm-max-iter", string(admm_max_iter),
            "--admm-eps-primal", string(admm_eps_primal),
            "--admm-eps-dual", string(admm_eps_dual),
        ]
        obj, _ = Lirion.main(args_dist)
        obj
    catch e
        @warn "Distributed ADMM-Hybrid-SOCP failed on $case_name: $e"
        track_failure_hybrid_socp(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized SOCP = $(round(cent_obj, digits=6)) | distributed-hybrid-SOCP = $(round(dist_obj, digits=6)) | max_iter=$(admm_max_iter), ρV²=$(rho_v), ρP=$(rho_flow_p), ρQ=$(rho_flow_q)"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_hybrid_socp_objective_tests(;
    case_list = distributed_hybrid_socp_zone_cases,
    label = "ADMM-Hybrid-SOCP vs SOCP",
    central_optimizer = _CENTRAL_DEFAULT,
    distributed_optimizer = _DISTRIB_DEFAULT,
    rtol = 0.05, atol = 1e-6,
    admm_eps_primal = 1e-4, admm_eps_dual = 1e-4,
)
    if isempty(case_list)
        @info "Empty hybrid-SOCP distributed case list. No tests will be run."
        return true
    end

    test_dir = joinpath(@__DIR__, "matpower_cases")
    selected_files = [
        joinpath(test_dir, fname)
        for fname in case_list
        if isfile(joinpath(test_dir, fname))
    ]

    n = length(selected_files)
    for (i, file) in enumerate(selected_files)
        case_name = basename(file)

        rho_v_case = 1.0e3
        rho_flow_p_case = 3.0e2
        rho_flow_q_case = 3.0e2
        admm_max_iter_case = 350
        admm_eps_primal_case = admm_eps_primal
        admm_eps_dual_case = admm_eps_dual

        if occursin("_3zones_", case_name) || occursin("_4zones_", case_name)
            admm_max_iter_case = max(admm_max_iter_case, 800)
        end
        # if occursin("case118_", case_name)
        #     admm_max_iter_case = max(admm_max_iter_case, 1500)
        # end
        # if occursin("case500_", case_name)
        #     admm_max_iter_case = max(admm_max_iter_case, 2000)
        # end

        @info "[$label] Rodando caso $i de $n: $case_name " *
              "(max_iter=$(admm_max_iter_case), ρV²=$(rho_v_case), ρP=$(rho_flow_p_case), ρQ=$(rho_flow_q_case), " *
              "eps_pr=$(admm_eps_primal_case), eps_dual=$(admm_eps_dual_case))"

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
                admm_eps_primal = admm_eps_primal_case,
                admm_eps_dual = admm_eps_dual_case,
            )
        end
    end
    return true
end