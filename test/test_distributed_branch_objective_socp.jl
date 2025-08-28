# ==============================================
# ADMM-Branch-SOCP vs SOCP (objetivo) – testes
# ==============================================
using Test
using PowerOptimization
using Glob
using SCS

const distributed_branch_socp_zone_cases = [
    # 14 buses
    "pglib_opf_case14_2zones_ieee_socp.m",
    "pglib_opf_case14_3zones_ieee_socp.m",
    # 24 buses
    "pglib_opf_case24_2zones_ieee_rts.m",
    "pglib_opf_case24_3zones_ieee_rts.m",
    # 30 buses
    "pglib_opf_case30_2zones_as.m",
    "pglib_opf_case30_3zones_as.m",
    # 57/118/500 
    # "pglib_opf_case57_2zones_ieee.m",
    # "pglib_opf_case57_4zones_ieee.m",
    # "pglib_opf_case118_2zones_ieee.m",
    # "pglib_opf_case500_2zones_goc.m",
]

function track_failure_branch_socp(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] =
            Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

function compare_distributed_branch_socp_vs_centralized(
    file::String;
    label::String,
    central_optimizer = SCS.Optimizer,        # SOCP central
    distributed_optimizer = SCS.Optimizer,    # ADMM-Branch-SOCP
    rtol = 0.05,
    atol = 1e-6,
    # ADMM penalties
    rho_flow_p::Float64 = 3.0e2,
    rho_flow_q::Float64 = 3.0e2,
    # ADMM settings
    admm_max_iter::Int = 1000,
    admm_eps_primal::Float64 = 1e-5,
    admm_eps_dual::Float64 = 1e-5,
)
    case_name = basename(file)

    # --- SOCP centralizado ---
    cent_obj = try
        args_cent = String[
            file,
            "--model", "socp",
            "--optimizer", string(central_optimizer),
        ]
        result_dict, _, _ = PowerOptimization.main(args_cent)
        result_dict["objective"]
    catch e
        @warn "Centralized SOCP failed on $case_name: $e"
        track_failure_branch_socp(label, "Centralized", case_name)
        return
    end

    # --- ADMM-Branch-SOCP distribuído ---
    dist_obj = try
        args_dist = String[
            file,
            "--model", "distributed-admm-branch-socp",
            "--optimizer", string(distributed_optimizer),

            "--admm-rho-flow-p", string(rho_flow_p),
            "--admm-rho-flow-q", string(rho_flow_q),

            "--admm-max-iter",   string(admm_max_iter),
            "--admm-eps-primal", string(admm_eps_primal),
            "--admm-eps-dual",   string(admm_eps_dual),
        ]
        obj, _ = PowerOptimization.main(args_dist)
        obj
    catch e
        @warn "Distributed ADMM-Branch-SOCP failed on $case_name: $e"
        track_failure_branch_socp(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized SOCP = $(round(cent_obj, digits=6)) | distributed-branch-SOCP = $(round(dist_obj, digits=6))"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_branch_socp_objective_tests(;
    case_list = distributed_branch_socp_zone_cases,
    label = "ADMM-Branch-SOCP vs SOCP",
    central_optimizer = SCS.Optimizer,
    distributed_optimizer = SCS.Optimizer,
    rtol = 0.05, atol = 1e-6,
    admm_eps_primal = 1e-5, admm_eps_dual = 1e-5,
)
    if isempty(case_list)
        @info "Empty distributed (branch-SOCP) case list. No tests will be run."
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

        # hiperparâmetros default validados
        rho_flow_p_case = 1e3
        rho_flow_q_case = 1e3
        admm_max_iter_case = 350
        admm_eps_primal_case = admm_eps_primal
        admm_eps_dual_case   = admm_eps_dual

        # ajustes leves por tamanho/zonas (mantém ρP/ρQ)
        if occursin("_3zones_", case_name) || occursin("_4zones_", case_name)
            admm_max_iter_case = max(admm_max_iter_case, 1200)
        end
        if occursin("case118_", case_name)
            admm_max_iter_case = max(admm_max_iter_case, 1500)
        end
        if occursin("case500_", case_name)
            admm_max_iter_case = max(admm_max_iter_case, 2000)
        end

        @info "[$label] Rodando caso $i de $n: $case_name " *
              "(max_iter=$(admm_max_iter_case), ρP=$(rho_flow_p_case), ρQ=$(rho_flow_q_case), " *
              "eps_pr=$(admm_eps_primal_case), eps_dual=$(admm_eps_dual_case))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_branch_socp_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer,
                distributed_optimizer = distributed_optimizer,
                rtol = rtol,
                atol = atol,
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