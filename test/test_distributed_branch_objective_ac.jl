# ==============================================
# ADMM-Branch-AC vs AC (objetivo) – testes
# ==============================================
const distributed_branch_ac_zone_cases = [
    # 14 buses
    "pglib_opf_case14_2zones_ieee.m",
    "pglib_opf_case14_3zones_ieee.m",
    # 24 buses
    "pglib_opf_case24_2zones_ieee_rts.m",
    "pglib_opf_case24_3zones_ieee_rts.m",
    # 30 buses
    "pglib_opf_case30_2zones_as.m",
    "pglib_opf_case30_3zones_as.m",
    # 57 buses
    # "pglib_opf_case57_2zones_ieee.m",
    # "pglib_opf_case57_4zones_ieee.m", 

    ####
    # "pglib_opf_case118_2zones_ieee.m",
    # "pglib_opf_case500_2zones_goc.m",
]

function track_failure_branch_ac(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] =
            Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

function compare_distributed_branch_ac_vs_centralized(
    file::String;
    label::String,
    central_optimizer = Ipopt.Optimizer,
    distributed_optimizer = Ipopt.Optimizer,
    rtol = 0.05,
    atol = 1e-6,
    admm_rho_flow = 1.0e5,
    admm_rho_flow_p = nothing,
    admm_rho_flow_q = nothing,
    admm_max_iter = 400,
    admm_eps_primal = 1e-6,
    admm_eps_dual = 1e-6,
)
    case_name = basename(file)

    # --- Centralizado AC ---
    cent_obj = try
        args_cent = String[
            file,
            "--model", "ac",
            "--optimizer", string(central_optimizer),
        ]
        result_dict, _, _ = PowerOptimization.main(args_cent)
        result_dict["objective"]
    catch e
        @warn "Centralized AC failed on $case_name: $e"
        track_failure_branch_ac(label, "Centralized", case_name)
        return
    end

    # --- Distribuído ADMM-Branch-AC ---
    dist_obj = try
        args_dist = String[
            file,
            "--model", "distributed-admm-branch-ac",
            "--optimizer", string(distributed_optimizer),
            "--admm-max-iter", string(admm_max_iter),
            "--admm-eps-primal", string(admm_eps_primal),
            "--admm-eps-dual", string(admm_eps_dual),
        ]
        if admm_rho_flow_p === nothing && admm_rho_flow_q === nothing
            push!(args_dist, "--admm-rho-flow")
            push!(args_dist, string(admm_rho_flow))
        else
            if admm_rho_flow_p !== nothing
                push!(args_dist, "--admm-rho-flow-p")
                push!(args_dist, string(admm_rho_flow_p))
            end
            if admm_rho_flow_q !== nothing
                push!(args_dist, "--admm-rho-flow-q")
                push!(args_dist, string(admm_rho_flow_q))
            end
        end

        obj, _ = PowerOptimization.main(args_dist)
        obj
    catch e
        @warn "Distributed ADMM-Branch-AC failed on $case_name: $e"
        track_failure_branch_ac(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized AC = $(round(cent_obj, digits=6)) | distributed-branch-AC = $(round(dist_obj, digits=6))"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_branch_ac_objective_tests(;
    case_list = distributed_branch_ac_zone_cases,
    label = "ADMM-Branch-AC vs AC",
    central_optimizer = Ipopt.Optimizer,
    distributed_optimizer = Ipopt.Optimizer,
    rtol = 0.05, atol = 1e-6,
    admm_eps_primal = 1e-6, admm_eps_dual = 1e-6,
)
    if isempty(case_list)
        @info "Empty distributed (branch-AC) case list. No tests will be run."
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
                (500, 1.0e6)
            elseif occursin("_2zones_", case_name)
                (400, 1.0e5)
            else
                (450, 1.0e5)
            end

        admm_eps_primal_case = admm_eps_primal
        admm_eps_dual_case = admm_eps_dual

        #OVERRIDE
        # if occursin("case118_", case_name)
        #     admm_max_iter_case = max(admm_max_iter_case, 600)
        #     admm_rho_case = 3.0e5
        #     admm_eps_primal_case = min(admm_eps_primal_case, 1e-7)
        #     admm_eps_dual_case = min(admm_eps_dual_case, 1e-7)
        # end

        if occursin("_3zones_", case_name)
            admm_eps_primal_case = 1e-8
            admm_eps_dual_case = 1e-8
        end

        @info "[$label] Rodando caso $i de $n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rho_flow=$(admm_rho_case), " *
              "eps_pr=$(admm_eps_primal_case), eps_dual=$(admm_eps_dual_case))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_branch_ac_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer,
                distributed_optimizer = distributed_optimizer,
                rtol = rtol,
                atol = atol,
                admm_rho_flow = admm_rho_case,
                admm_max_iter = admm_max_iter_case,
                admm_eps_primal = admm_eps_primal_case,
                admm_eps_dual = admm_eps_dual_case,
            )
        end
    end
    return true
end