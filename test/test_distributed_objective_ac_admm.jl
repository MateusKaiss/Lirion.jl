# ==============================================
# ADMM-AC vs AC (objetivo) – testes (compat μ/τ)
# ==============================================

# Casos particionados (2/3/4 zonas)
const distributed_ac_zone_cases = [
    # 14 buses
    "pglib_opf_case14_2zones_ieee.m",
    "pglib_opf_case14_3zones_ieee.m",
    # 24 buses
    "pglib_opf_case24_2zones_ieee_rts.m",
    "pglib_opf_case24_3zones_ieee_rts.m",
    # 30 buses
    "pglib_opf_case30_2zones_as.m",
    "pglib_opf_case30_3zones_as.m",
    # 57/118/500 (exemplos comentados)
    # "pglib_opf_case57_2zones_ieee.m",
    # "pglib_opf_case57_4zones_ieee.m",
    # "pglib_opf_case118_2zones_ieee.m",
    # "pglib_opf_case500_2zones_goc.m",
]

function track_failure_ac(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] =
            Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end
function compare_distributed_ac_vs_centralized(
    file::String;
    label::String,
    central_optimizer = Ipopt.Optimizer,
    distributed_optimizer = Ipopt.Optimizer,
    rtol = 0.05,
    atol = 1e-6,
    # parâmetros padrões
    admm_rho = 1_000.0,
    admm_max_iter = 350,
    admm_eps_primal = 1e-6,
    admm_eps_dual = 1e-6,
    admm_mu_incr = 2.0,
    admm_mu_decr = 2.0,
    admm_tau_incr = 1.05,
    admm_tau_decr = 0.98,
    admm_rho_theta::Union{Nothing,Real} = nothing,
    admm_rho_v::Union{Nothing,Real} = nothing,
    admm_adaptive::Union{Nothing,String} = nothing, # "on" | "off" | nothing (não passa flag)
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
        track_failure_ac(label, "Centralized", case_name)
        return
    end

    # --- Distribuído ADMM-AC ---
    dist_obj = try
        args_dist = String[
            file,
            "--model", "distributed-admm-ac",
            "--optimizer", string(distributed_optimizer),
            "--admm-max-iter", string(admm_max_iter),
            "--admm-eps-primal", string(admm_eps_primal),
            "--admm-eps-dual", string(admm_eps_dual),
            "--admm-mu-incr", string(admm_mu_incr),
            "--admm-mu-decr", string(admm_mu_decr),
            "--admm-tau-incr", string(admm_tau_incr),
            "--admm-tau-decr", string(admm_tau_decr),
        ]
        if admm_rho_theta !== nothing
            push!(args_dist, "--admm-rho-theta"); push!(args_dist, string(admm_rho_theta))
        end
        if admm_rho_v !== nothing
            push!(args_dist, "--admm-rho-v");     push!(args_dist, string(admm_rho_v))
        end
        if admm_rho_theta === nothing && admm_rho_v === nothing
            push!(args_dist, "--admm-rho");       push!(args_dist, string(admm_rho))
        end
        if admm_adaptive !== nothing
            push!(args_dist, "--admm-adaptive");  push!(args_dist, admm_adaptive)
        end

        obj, _ = PowerOptimization.main(args_dist)
        obj
    catch e
        @warn "Distributed ADMM-AC failed on $case_name: $e"
        track_failure_ac(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized AC = $(round(cent_obj, digits=6)) | distributed ADMM-AC = $(round(dist_obj, digits=6))"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_ac_objective_tests(;
    case_list = distributed_ac_zone_cases,
    label = "ADMM-AC vs AC",
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
        @info "Empty distributed (AC) case list. No tests will be run."
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
                (600, 1_000.0)
            elseif occursin("_2zones_", case_name)
                (600, 1_000.0)
            else
                (600, 1_000.0)
            end

        admm_eps_primal_case = admm_eps_primal
        admm_eps_dual_case   = admm_eps_dual

        admm_rho_theta_case = nothing
        admm_rho_v_case     = nothing
        admm_adaptive_case  = nothing

        if occursin("case118_", case_name)
            admm_max_iter_case = max(admm_max_iter_case, 500)
        end
        if occursin("case500_", case_name)
            admm_max_iter_case = max(admm_max_iter_case, 800)
            admm_rho_case = 1_000_000.0
            admm_eps_primal_case = min(admm_eps_primal_case, 1e-4)
            admm_eps_dual_case   = min(admm_eps_dual_case,   1e-4)
        end

        if case_name == "pglib_opf_case118_2zones_ieee.m"
            admm_max_iter_case   = 500
            admm_eps_primal_case = 1e-4
            admm_eps_dual_case   = 1e-4
            admm_rho_theta_case  = 1e5
            admm_rho_v_case      = 1e5
            admm_adaptive_case   = "on"
            admm_mu_incr = 5.0
            admm_mu_decr = 5.0
            admm_tau_incr = 1.25
            admm_tau_decr = 0.75
        end

        @info "[$label] Rodando caso $i de $n: $case_name (max_iter=$(admm_max_iter_case), rho=$(admm_rho_case), ρθ=$(admm_rho_theta_case), ρV=$(admm_rho_v_case), eps_pr=$(admm_eps_primal_case), eps_dual=$(admm_eps_dual_case), μ=(incr=$(admm_mu_incr),decr=$(admm_mu_decr)), τ=(incr=$(admm_tau_incr),decr=$(admm_tau_decr)), adaptive=$(admm_adaptive_case))"
        @testset "$(label) - $(case_name)" begin
            compare_distributed_ac_vs_centralized(
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
                admm_mu_incr = admm_mu_incr,
                admm_mu_decr = admm_mu_decr,
                admm_tau_incr = admm_tau_incr,
                admm_tau_decr = admm_tau_decr,
                admm_rho_theta = admm_rho_theta_case,
                admm_rho_v = admm_rho_v_case,
                admm_adaptive = admm_adaptive_case,
            )
        end
    end
    return true
end

