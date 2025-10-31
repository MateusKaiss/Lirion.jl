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
    # 57 buses (descomente quando quiser)
    # "pglib_opf_case57_2zones_ieee.m",
    # "pglib_opf_case57_4zones_ieee.m",
    # 118/500 buses
    # "pglib_opf_case118_2zones_ieee.m",
    # "pglib_opf_case118_3zones_ieee.m",
    "pglib_opf_case118_5zones_ieee.m",
    # # 500 buses
    # "pglib_opf_case500_3zones_goc.m",
    # "pglib_opf_case500_4zones_goc.m",
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
    admm_eps_primal = 1e-4,
    admm_eps_dual = 1e-4,
    admm_mu_incr = 2.0,
    admm_mu_decr = 2.0,
    admm_tau_incr = 1.05,
    admm_tau_decr = 0.98,
    admm_adaptive::Union{Nothing, String} = nothing,
    additional_flags::Vector{String} = String[],
)
    case_name = basename(file)

    # --- Centralizado AC ---
    cent_obj = try
        args_cent = String[
            file,
            "--model", "ac",
            "--optimizer", string(central_optimizer),
        ]
        result_dict, _, _ = Lirion.main(args_cent)
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
            "--admm-mu-incr", string(admm_mu_incr),
            "--admm-mu-decr", string(admm_mu_decr),
            "--admm-tau-incr", string(admm_tau_incr),
            "--admm-tau-decr", string(admm_tau_decr),
        ]
        if admm_rho_flow_p === nothing && admm_rho_flow_q === nothing
            push!(args_dist, "--admm-rho-flow", string(admm_rho_flow))
        else
            if admm_rho_flow_p !== nothing
                push!(args_dist, "--admm-rho-flow-p", string(admm_rho_flow_p))
            end
            if admm_rho_flow_q !== nothing
                push!(args_dist, "--admm-rho-flow-q", string(admm_rho_flow_q))
            end
        end
        if admm_adaptive !== nothing
            push!(args_dist, "--admm-adaptive", admm_adaptive)
        end

        append!(args_dist, additional_flags)

        obj, _ = Lirion.main(args_dist)
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
    admm_eps_primal = 1e-4, admm_eps_dual = 1e-4,
    admm_mu_incr = 2.0,
    admm_mu_decr = 2.0,
    admm_tau_incr = 1.05,
    admm_tau_decr = 0.98,
)
    if isempty(case_list)
        @info "Empty distributed (branch-AC) case list. No tests will be run."
        return true
    end

    test_dir = joinpath(@__DIR__, "matpower_cases")
    selected_files =
        [joinpath(test_dir, f) for f in case_list if isfile(joinpath(test_dir, f))]

    n = length(selected_files)
    for (i, file) in enumerate(selected_files)
        case_name = basename(file)

        admm_max_iter_case, admm_rho_case =
            if occursin("_3zones_", case_name)
                (400, 1.0e6)
            elseif occursin("_2zones_", case_name)
                (350, 1.0e5)
            else
                (300, 1.0e5)
            end
        admm_eps_primal_case = admm_eps_primal
        admm_eps_dual_case = admm_eps_dual
        μ_incr_case = admm_mu_incr
        μ_decr_case = admm_mu_decr
        τ_incr_case = admm_tau_incr
        τ_decr_case = admm_tau_decr
        adaptive_case = nothing
        extra_flags = String[]

        if occursin("_3zones_", case_name)
            admm_eps_primal_case = 1e-4
            admm_eps_dual_case = 1e-4
        end
        if case_name == "pglib_opf_case118_2zones_ieee.m"
            admm_max_iter_case = 500
            admm_eps_primal_case = 1e-4
            admm_eps_dual_case = 1e-4
            admm_rho_case = 1.0e5
            adaptive_case = "on"
            μ_incr_case, μ_decr_case = 2.0, 2.0
            τ_incr_case, τ_decr_case = 1.05, 0.95
        end
        if case_name in
           ("pglib_opf_case118_5zones_ieee.m", "pglib_opf_case118_3zones_ieee.m")
            admm_max_iter_case = 500
            admm_eps_primal_case = 1e-4
            admm_eps_dual_case = 1e-4
            admm_rho_case = 1.0e5
            adaptive_case = "on"
            μ_incr_case, μ_decr_case = 3.0, 3.0
            τ_incr_case, τ_decr_case = 1.1, 0.9
        end

        if case_name == "pglib_opf_case500_3zones_goc.m" ||
           case_name == "pglib_opf_case500_4zones_goc.m"
            admm_max_iter_case = 700
            admm_eps_primal_case = 1e-4
            admm_eps_dual_case = 1e-4
            adaptive_case = "on"
            μ_incr_case, μ_decr_case = 3.0, 3.0
            τ_incr_case, τ_decr_case = 1.3, 0.7
            extra_flags = String[
                "--admm-rho-theta", "1e5",
                "--admm-rho-v", "1e5",
            ]
        end
        
        if case_name == "pglib_opf_case2000_5zones_goc.m"
            admm_max_iter_case   = 700
            admm_eps_primal_case = 1e-4
            admm_eps_dual_case   = 1e-4
            admm_rho_case        = 1.0e4   # rho único (P e Q)
            adaptive_case        = "on"
            μ_incr_case, μ_decr_case = 3.0, 3.0
            τ_incr_case, τ_decr_case = 1.1, 0.9
        end

        if case_name == "pglib_opf_case3970_25zones_goc.m"
            admm_max_iter_case   = 700
            admm_eps_primal_case = 3e-4
            admm_eps_dual_case   = 3e-4
            admm_rho_case        = 1.0e4   # rho único (P e Q)
            adaptive_case        = "on"
            μ_incr_case, μ_decr_case = 3.0, 3.0
            τ_incr_case, τ_decr_case = 1.1, 0.9
        end

        if case_name == "pglib_opf_case10000_goc.m"
            admm_max_iter_case   = 800
            admm_eps_primal_case = 3e-4
            admm_eps_dual_case   = 3e-4
            admm_rho_case        = 1.0e3   # rho único (P e Q). 
            adaptive_case        = "on"
            μ_incr_case, μ_decr_case = 3.0, 3.0
            τ_incr_case, τ_decr_case = 1.05, 0.95
        end


        @info "[$label] Rodando caso $i de $n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rho_flow=$(admm_rho_case), " *
              "eps_pr=$(admm_eps_primal_case), eps_dual=$(admm_eps_dual_case), " *
              "μ=(incr=$(μ_incr_case),decr=$(μ_decr_case)), " *
              "τ=(incr=$(τ_incr_case),decr=$(τ_decr_case)), " *
              "adaptive=$(adaptive_case))"

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
                admm_mu_incr = μ_incr_case,
                admm_mu_decr = μ_decr_case,
                admm_tau_incr = τ_incr_case,
                admm_tau_decr = τ_decr_case,
                admm_adaptive = adaptive_case,
                additional_flags = extra_flags,
            )
        end
    end
    return true
end
