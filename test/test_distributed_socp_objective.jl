# ============================================
# ADMM-SOCP vs SOCP (objetivo) – testes
# ============================================
using Test
using Lirion
using Glob
using SCS
using Gurobi
using COSMO
using Clarabel

const ADMM_MAX_IT = 5000

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
        prefer::Symbol = Symbol(
            lowercase(get(ENV, "LIRION_SOCP_PREFER", "auto")),
        ),
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
    "pglib_opf_case118_ieee.m",
    "pglib_opf_case200_activ.m",

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
function compare_distributed_socp_vs_centralized(
    file::String;
    label::String,
    central_optimizer = _CENTRAL_DEFAULT,
    distributed_optimizer = _DISTRIB_DEFAULT,
    rtol = 0.05,
    atol = 1e-6,
    rho_v::Float64 = 1.0e3,
    rho_flow_p::Float64 = 3.0e2,
    rho_flow_q::Float64 = 3.0e2,
    admm_eps_abs::Float64 = 1e-4,
    admm_eps_rel::Float64 = 1e-4,
)
    case_name = basename(file)

    # =====================================================================
    # Caso especial: 2000-bus (usa Clarabel tunado no centralizado)
    # =====================================================================
    if case_name == "pglib_opf_case2000_5zones_goc.m"
        @info "⚙️ Usando Clarabel robusto no caso especial $case_name"

        # --- Centralizado SOCP (Clarabel robusto) ---
        cent_obj = try
            args_cent = String[
                file,
                "--model", "socp",
                "--optimizer",
                """optimizer_with_attributes(Clarabel.Optimizer,
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
                )""",
            ]
            result_dict, _, _ = Lirion.main(args_cent)
            result_dict["objective"]
        catch e
            @warn "Centralized SOCP failed on $case_name: $e"
            track_failure_socp(label, "Centralized", case_name)
            return
        end
ADMM_MAX_IT
        # --- Distribuído ADMM-SOCP (Clarabel padrão) ---
        dist_obj = try
            args_dist = String[
                file,
                "--model", "distributed-admm-socp",
                "--optimizer", "Clarabel.Optimizer",
                "--admm-max-iter", "700",
                "--admm-rho-v", "1e6",
                "--admm-rho-flow-p", "1e6",
                "--admm-rho-flow-q", "1e6",
                "--admm-eps-abs", "1e-4",
                "--admm-eps-rel", "1e-4",
                "--admm-adaptive", "on",
                "--admm-mu-incr", "5.0",
                "--admm-mu-decr", "5.0",
                "--admm-tau-incr", "1.1",
                "--admm-tau-decr", "0.9",
            ]
            obj, _ = Lirion.main(args_dist)
            obj
        catch e
            @warn "Distributed ADMM-SOCP failed on $case_name: $e"
            track_failure_socp(label, "Distributed", case_name)
            return
        end

        @info "File $case_name | centralized SOCP = $(round(cent_obj, digits=6)) | distributed-ADMM-SOCP = $(round(dist_obj, digits=6))"
        @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
        return
    end

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
        track_failure_socp(label, "Centralized", case_name)
        return
    end

    dist_obj = try
        args_dist = String[
            file,
            "--model", "distributed-admm-socp",
            "--optimizer", string(distributed_optimizer),
            "--admm-rho-v", string(rho_v),
            "--admm-rho-flow-p", string(rho_flow_p),
            "--admm-rho-flow-q", string(rho_flow_q),
            "--admm-max-iter", string(ADMM_MAX_IT),
            "--admm-eps-abs", string(admm_eps_abs),
            "--admm-eps-rel", string(admm_eps_rel),
        ]
        obj, _ = Lirion.main(args_dist)
        obj
    catch e
        @warn "Distributed ADMM-SOCP failed on $case_name: $e"
        track_failure_socp(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized SOCP = $(round(cent_obj, digits=6)) | distributed-ADMM-SOCP = $(round(dist_obj, digits=6)) | ADMM_MAX_IT=$(ADMM_MAX_IT)"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_socp_objective_tests(;
    case_list = distributed_socp_zone_cases,
    label = "ADMM-SOCP vs SOCP",
    central_optimizer = _CENTRAL_DEFAULT,
    distributed_optimizer = _DISTRIB_DEFAULT,
    rtol = 0.05, atol = 1e-6,
    admm_eps_abs = 1e-4, admm_eps_rel = 1e-4,
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

        # penalidades constantes
        rho_v_case = 1.0e7
        rho_flow_p_case = 1.0e7
        rho_flow_q_case = 1.0e7

        admm_max_iter_case = 5000
        admm_eps_abs_case = admm_eps_abs
        admm_eps_rel_case = admm_eps_rel

        @info "[$label] Rodando caso $i de $n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rhoV²=$(rho_v_case), rhoP=$(rho_flow_p_case), rhoQ=$(rho_flow_q_case), " *
              "eps_pr=$(admm_eps_abs_case), eps_dual=$(admm_eps_rel_case))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_socp_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer,
                distributed_optimizer = distributed_optimizer,
                rtol = rtol,
                atol = atol,
                rho_v = rho_v_case,
                rho_flow_p = rho_flow_p_case,
                rho_flow_q = rho_flow_q_case,
                admm_eps_abs = admm_eps_abs_case,
                admm_eps_rel = admm_eps_rel_case,
            )
        end
    end
    return true
end
