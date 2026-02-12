# ==============================================
# ADMM-Branch-DC vs DC (objetivo) – testes
# ==============================================

# Casos particionados (2/3/4 zonas)
const distributed_branch_zone_cases = [
    # "pglib_opf_case73_ieee_rts__api.m",
    # "pglib_opf_case73_ieee_rts__sad.m",
    ## "pglib_opf_case57_ieee__api.m",
    # "pglib_opf_case57_ieee__sad.m",

    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
    "pglib_opf_case73_ieee_rts.m",
    "pglib_opf_case57_ieee.m",
    "pglib_opf_case118_ieee.m",
    # "pglib_opf_case200_activ.m",
    # "pglib_opf_case500_goc.m",

]

function track_failure_branch(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] =
            Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

function compare_distributed_branch_vs_centralized(
    file::String;
    label::String,
    central_optimizer = get_available_optimizer(OSQP.Optimizer, Ipopt.Optimizer),
    distributed_optimizer = nothing,
    rtol = 0.01,
    atol = 1e-6,
    admm_rho = 1_000.0,
    admm_max_iter = 5000,
    admm_eps_abs = 1e-4,
    admm_eps_rel = 1e-4,
    admm_alpha::Float64 = 1.0,
    admm_mu_incr::Float64 = 3.0,
    admm_mu_decr::Float64 = 3.0,
    admm_tau_incr::Float64 = 1.05,
    admm_tau_decr::Float64 = 0.98,
    admm_warm_start::Bool = false,
)
    case_name = basename(file)

    # ----------------------------
    # Centralizado DC -> objective
    # ----------------------------
    cent_obj = try
        args_cent = String[
            file,
            "--model", "dc",
            "--optimizer", string(central_optimizer),
        ]
        result_dict, _, _ = Lirion.main(args_cent)
        Float64(result_dict["objective"])
    catch e
        @warn "Centralized DC failed on $case_name: $e"
        track_failure_branch(label, "Centralized", case_name)
        return
    end

    # --------------------------------------
    # Distribuído ADMM-Branch-DC -> NamedTuple
    # --------------------------------------
    dist_obj = try
        args_dist = String[
            file,
            "--model", "distributed-admm-branch-dc",
            "--admm-rho", string(admm_rho),
            "--admm-max-iter", string(admm_max_iter),
            "--admm-eps-abs", string(admm_eps_abs),
            "--admm-eps-rel", string(admm_eps_rel),
            "--admm-mu-incr", string(admm_mu_incr),
            "--admm-mu-decr", string(admm_mu_decr),
            "--admm-tau-incr", string(admm_tau_incr),
            "--admm-tau-decr", string(admm_tau_decr),
        ]

        if admm_alpha != 1.0
            push!(args_dist, "--admm-alpha", string(admm_alpha))
        end
        if admm_warm_start
            push!(args_dist, "--admm-warm-start")
        end
        if distributed_optimizer !== nothing
            push!(args_dist, "--optimizer", string(distributed_optimizer))
        end

        res = Lirion.main(args_dist)

        # aqui é o ponto: res é NamedTuple, pega o custo certo
        Float64(res.total_real_cost)
        # (se você quiser o outro: Float64(res.total_expr_cost))
    catch e
        @warn "Distributed ADMM-Branch-DC failed on $case_name: $e"
        track_failure_branch(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized = $(round(cent_obj, digits=6)) | distributed-branch = $(round(dist_obj, digits=6))"

    abs_err = abs(dist_obj - cent_obj)
    rel_err = abs_err / max(abs(cent_obj), 1.0)
    @info "Errors: abs=$(abs_err) | rel=$(rel_err)"

    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end


"""
Roda a suíte de testes para ADMM-Branch-DC vs DC nos casos de zonas.
"""
function run_distributed_branch_objective_tests(;
    case_list = distributed_branch_zone_cases,
    label = "ADMM-Branch-DC vs DC",
    central_optimizer = get_available_optimizer(OSQP.Optimizer, Ipopt.Optimizer),
    distributed_optimizer = nothing,
    rtol = 0.01, atol = 1e-6,
    admm_eps_abs = 1e-4, admm_eps_rel = 1e-4,
    admm_alpha::Float64 = 1.0,
    admm_mu_incr::Float64 = 3.0,
    admm_mu_decr::Float64 = 3.0,
    admm_tau_incr::Float64 = 1.05,
    admm_tau_decr::Float64 = 0.98,
)
    if isempty(case_list)
        @info "Empty distributed (branch) case list. No tests will be run."
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
        admm_max_iter_case = 5000
        admm_rho_case = 100000.0
        admm_eps_abs_case = admm_eps_abs
        admm_eps_rel_case = admm_eps_rel
        admm_alpha_case = admm_alpha
        admm_mu_incr_case = admm_mu_incr
        admm_mu_decr_case = admm_mu_decr
        admm_tau_incr_case = admm_tau_incr
        admm_tau_decr_case = admm_tau_decr
        distributed_optimizer_case = distributed_optimizer
        admm_warm_start_case = false

        # --------------------------------------------
        # Override por caso (ex: case197_snem)
        # --------------------------------------------
        if case_name == "pglib_opf_case500_goc.m"
            admm_rho_case = 5000.0
            central_optimizer_case = Ipopt.Optimizer
            distributed_optimizer_case = Ipopt.Optimizer
        end




        @info "[$label] Rodando caso $i de $n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rho=$(admm_rho_case), " *
              "eps_pr=$(admm_eps_abs_case), eps_dual=$(admm_eps_rel_case), " *
              "alpha=$(admm_alpha_case), mu_incr=$(admm_mu_incr_case), mu_decr=$(admm_mu_decr_case), " *
              "tau_incr=$(admm_tau_incr_case), tau_decr=$(admm_tau_decr_case), " *
              "warm=$(admm_warm_start_case), optimizer=$(distributed_optimizer_case === nothing ? "default" : string(distributed_optimizer_case)))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_branch_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer,
                distributed_optimizer = distributed_optimizer_case,
                rtol = rtol,
                atol = atol,
                admm_rho = admm_rho_case,
                admm_max_iter = admm_max_iter_case,
                admm_eps_abs = admm_eps_abs_case,
                admm_eps_rel = admm_eps_rel_case,
                admm_alpha = admm_alpha_case,
                admm_mu_incr = admm_mu_incr_case,
                admm_mu_decr = admm_mu_decr_case,
                admm_tau_incr = admm_tau_incr_case,
                admm_tau_decr = admm_tau_decr_case,
                admm_warm_start = admm_warm_start_case,
            )
        end
    end
    return true
end
