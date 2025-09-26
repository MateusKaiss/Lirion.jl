# =============================
# Casos particionados (2/3 zonas)
# =============================
const distributed_zone_cases = [
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
    "pglib_opf_case57_2zones_ieee.m",
    "pglib_opf_case57_4zones_ieee.m",
    # 118 buses
    "pglib_opf_case118_2zones_ieee.m",
    "pglib_opf_case118_3zones_ieee.m",
    "pglib_opf_case118_5zones_ieee.m",
    # 500 buses
    "pglib_opf_case500_3zones_goc.m",
    "pglib_opf_case500_4zones_goc.m",
    # 793
    # "pglib_opf_case793_11zones_goc.m",
    # 2000
    "pglib_opf_case2000_11zonesbus_goc.m",
]

function track_failure_dist(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] =
            Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

function compare_distributed_vs_centralized(
    file::String;
    label::String,
    central_optimizer = get_available_optimizer(OSQP.Optimizer, Ipopt.Optimizer),
    distributed_optimizer::Union{Nothing, Function} = nothing,
    rtol = 0.05,
    atol = 1e-6,
    admm_rho = 1_000.0,
    admm_max_iter = 350,
    admm_eps_primal = 1e-6,
    admm_eps_dual = 1e-6,
    # α=1.0 => over-relaxation inativo
    admm_alpha::Float64 = 1.0,
    admm_mu_incr::Float64 = 3.0,
    admm_mu_decr::Float64 = 3.0,
    admm_tau_incr::Float64 = 1.05,
    admm_tau_decr::Float64 = 0.98,
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
        track_failure_dist(label, "Centralized", case_name)
        return
    end

    # --- Distribuído ADMM-DC (bus) ---
    dist_obj = try
        args_dist = String[
            file,
            "--model", "distributed-admm-dc",
            "--admm-rho", string(admm_rho),
            "--admm-max-iter", string(admm_max_iter),
            "--admm-eps-primal", string(admm_eps_primal),
            "--admm-eps-dual", string(admm_eps_dual),
            "--admm-alpha", string(admm_alpha),
            "--admm-mu-incr", string(admm_mu_incr),
            "--admm-mu-decr", string(admm_mu_decr),
            "--admm-tau-incr", string(admm_tau_incr),
            "--admm-tau-decr", string(admm_tau_decr),
        ]
        # Só passamos --optimizer se o usuário quiser sobrescrever o default (OSQP)
        if !isnothing(distributed_optimizer)
            push!(args_dist, "--optimizer")
            push!(args_dist, string(distributed_optimizer))
        end
        obj, _ = PowerOptimization.main(args_dist)
        obj
    catch e
        @warn "Distributed ADMM-DC failed on $case_name: $e"
        track_failure_dist(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized = $(round(cent_obj, digits=6)) | distributed = $(round(dist_obj, digits=6))"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_objective_tests(;
    case_list = distributed_zone_cases,
    label = "ADMM-DC (bus) vs DC",
    central_optimizer = get_available_optimizer(OSQP.Optimizer, Ipopt.Optimizer),
    distributed_optimizer::Union{Nothing, Function} = nothing,
    rtol = 0.05, atol = 1e-6,
    admm_eps_primal = 1e-4, admm_eps_dual = 1e-4,
    admm_alpha::Float64 = 1.0,  # over-relaxation inativo
)
    if isempty(case_list)
        @info "Empty distributed case list. No tests will be run."
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

        # defaults por zonas
        admm_max_iter_case, admm_rho_case =
            if occursin("_3zones_", case_name)
                (700, 1_000.0)
            elseif occursin("_2zones_", case_name)
                (300, 1_000.0)
            else
                (350, 1_000.0)
            end

        admm_eps_primal_case = admm_eps_primal
        admm_eps_dual_case   = admm_eps_dual

        # parâmetros adaptativos: defaults gerais
        admm_mu_incr_case = 3.0
        admm_mu_decr_case = 3.0
        admm_tau_incr_case = 1.05
        admm_tau_decr_case = 0.98

        # ---- OVERRIDES pontuais ----
        if occursin("case118_", case_name)
            admm_max_iter_case = 700
            admm_rho_case = 1_000_000.0
            admm_eps_primal_case = 1e-4
            admm_eps_dual_case   = 1e-4
        end
        if occursin("case2000_11zonesbus", case_name)
            admm_max_iter_case   = 900
            admm_rho_case        = 10_000.0
            admm_eps_primal_case = 1e-4
            admm_eps_dual_case   = 1e-4
            admm_mu_incr_case    = 3.0
            admm_mu_decr_case    = 3.0
            admm_tau_incr_case   = 1.2
            admm_tau_decr_case   = 0.8
        end

        @info "[$label] Rodando caso $i de $n: $case_name (max_iter=$(admm_max_iter_case), rho=$(admm_rho_case), eps_pr=$(admm_eps_primal_case), eps_dual=$(admm_eps_dual_case), alpha=$(admm_alpha), μ_incr=$(admm_mu_incr_case), μ_decr=$(admm_mu_decr_case), τ_incr=$(admm_tau_incr_case), τ_decr=$(admm_tau_decr_case))"
        @testset "$(label) - $(case_name)" begin
            compare_distributed_vs_centralized(
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
                admm_alpha = admm_alpha,
                admm_mu_incr = admm_mu_incr_case,
                admm_mu_decr = admm_mu_decr_case,
                admm_tau_incr = admm_tau_incr_case,
                admm_tau_decr = admm_tau_decr_case,
            )
        end
    end
    return true
end
