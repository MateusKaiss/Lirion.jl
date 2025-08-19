# ==============================================
# ADMM-Branch-DC vs DC (objetivo) – testes
# ==============================================

# Casos particionados (2/3/4 zonas) — mesmos do teste por barramento
const distributed_branch_zone_cases = [
    # 14 buses
    "pglib_opf_case14_2zones_ieee.m",
    "pglib_opf_case14_3zones_ieee.m",
    # # 24 buses
    "pglib_opf_case24_2zones_ieee_rts.m",
    "pglib_opf_case24_3zones_ieee_rts.m",
    # # 30 buses
    "pglib_opf_case30_2zones_as.m",
    "pglib_opf_case30_3zones_as.m",
    # 57 buses
    # "pglib_opf_case57_2zones_ieee.m",
    # ########"pglib_opf_case57_4zones_ieee.m", #####problem
    # # # 118 buses
    # "pglib_opf_case118_2zones_ieee.m",
    # # # 500 buses
    # "pglib_opf_case500_2zones_goc.m",
]

function track_failure_branch(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] =
            Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

"""
Compara objetivo Centralizado DC vs Distribuído ADMM-Branch-DC no mesmo caso.
"""
function compare_distributed_branch_vs_centralized(
    file::String;
    label::String,
    central_optimizer = Ipopt.Optimizer,
    distributed_optimizer = Ipopt.Optimizer,
    rtol = 0.05,
    atol = 1e-6,
    admm_rho = 1_000.0,
    admm_max_iter = 350,
    admm_eps_primal = 1e-6,
    admm_eps_dual = 1e-6,
)
    case_name = basename(file)

    # --- Centralizado DC ---
    cent_obj = try
        args_cent = String[
            file,
            "--model", "dc",
            "--optimizer", string(central_optimizer),
        ]
        result_dict, _, _ = PowerOptimization.main(args_cent)
        result_dict["objective"]
    catch e
        @warn "Centralized DC failed on $case_name: $e"
        track_failure_branch(label, "Centralized", case_name)
        return
    end

    # --- Distribuído ADMM-Branch-DC ---
    dist_obj = try
        args_dist = String[
            file,
            "--model", "distributed-admm-branch-dc",
            "--optimizer", string(distributed_optimizer),
            "--admm-rho", string(admm_rho),
            "--admm-max-iter", string(admm_max_iter),
            "--admm-eps-primal", string(admm_eps_primal),
            "--admm-eps-dual", string(admm_eps_dual),
        ]
        obj, _ = PowerOptimization.main(args_dist)
        obj
    catch e
        @warn "Distributed ADMM-Branch-DC failed on $case_name: $e"
        track_failure_branch(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized = $(round(cent_obj, digits=6)) | distributed-branch = $(round(dist_obj, digits=6))"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

"""
Roda a suíte de testes para ADMM-Branch-DC vs DC nos casos de zonas.
"""
function run_distributed_branch_objective_tests(;
    case_list = distributed_branch_zone_cases,
    label = "ADMM-Branch-DC vs DC",
    central_optimizer = get_available_optimizer(Gurobi.Optimizer, Ipopt.Optimizer),
    distributed_optimizer = Ipopt.Optimizer,
    rtol = 0.05, atol = 1e-6,
    admm_eps_primal = 1e-6, admm_eps_dual = 1e-6,
)
    if isempty(case_list)
        @info "Empty distributed (branch) case list. No tests will be run."
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

        admm_max_iter_case, admm_rho_case =
            if occursin("_3zones_", case_name)
                (400, 1_000.0)
            elseif occursin("_2zones_", case_name)
                (350, 1_000.0)
            else
                (300, 1_000.0)
            end

        admm_eps_primal_case = admm_eps_primal
        admm_eps_dual_case = admm_eps_dual

        if occursin("case500_", case_name)
            admm_max_iter_case = 600
            admm_rho_case = 1_000_000.0
            admm_eps_primal_case = 1e-7
            admm_eps_dual_case = 1e-7
        end

        @info "[$label] Rodando caso $i de $n: $case_name (max_iter=$(admm_max_iter_case), rho=$(admm_rho_case), eps_pr=$(admm_eps_primal_case), eps_dual=$(admm_eps_dual_case))"
        @testset "$(label) - $(case_name)" begin
            compare_distributed_branch_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer,
                distributed_optimizer = distributed_optimizer,
                rtol = rtol,
                atol = atol,
                admm_rho = admm_rho_case,
                admm_max_iter = admm_max_iter_case,
                admm_eps_primal = admm_eps_primal_case,
                admm_eps_dual = admm_eps_dual_case,
            )
        end
    end
    return true
end
