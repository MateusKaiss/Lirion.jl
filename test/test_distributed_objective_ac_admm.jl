# ==============================================
# ADMM-AC vs AC (objetivo) – testes (compat mu/tau)
# ==============================================

# Casos particionados (2/3/4 zonas)
const distributed_ac_zone_cases = [
    # "pglib_opf_case73_ieee_rts__api.m",
    # "pglib_opf_case73_ieee_rts__sad.m",
    ## "pglib_opf_case57_ieee__api.m",
    # "pglib_opf_case57_ieee__sad.m",

    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",

    "pglib_opf_case73_ieee_rts.m",
    "pglib_opf_case118_ieee.m",
    # "pglib_opf_case200_activ.m",
    
    # "pglib_opf_case57_ieee.m",
    # "pglib_opf_case500_goc.m",
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
    rtol = 0.02,
    atol = 1e-6,

    admm_rho = 100.0,
    admm_max_iter = 5000,
    admm_eps_abs = 1e-4,
    admm_eps_rel = 1e-4,
    admm_mu_incr = 2.0,
    admm_mu_decr = 2.0,
    admm_tau_incr = 1.05,
    admm_tau_decr = 0.98,
    admm_rho_theta::Union{Nothing, Real} = nothing,
    admm_rho_v::Union{Nothing, Real} = nothing,
    admm_adaptive::Union{Nothing, String} = nothing, # "on" | "off" | nothing
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
            "--admm-eps-abs", string(admm_eps_abs),
            "--admm-eps-rel", string(admm_eps_rel),

            "--admm-mu-incr", string(admm_mu_incr),
            "--admm-mu-decr", string(admm_mu_decr),
            "--admm-tau-incr", string(admm_tau_incr),
            "--admm-tau-decr", string(admm_tau_decr),
        ]

        if admm_rho_theta !== nothing
            push!(args_dist, "--admm-rho-theta"); push!(args_dist, string(admm_rho_theta))
        end
        if admm_rho_v !== nothing
            push!(args_dist, "--admm-rho-v"); push!(args_dist, string(admm_rho_v))
        end
        if admm_rho_theta === nothing && admm_rho_v === nothing
            push!(args_dist, "--admm-rho"); push!(args_dist, string(admm_rho))
        end
        if admm_adaptive !== nothing
            push!(args_dist, "--admm-adaptive"); push!(args_dist, admm_adaptive)
        end

        obj, _ = Lirion.main(args_dist)
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
    rtol = 0.01, atol = 1e-6,
    admm_eps_abs = 1e-4, admm_eps_rel = 1e-4,
    admm_mu_incr = 2.0,
    admm_mu_decr = 2.0,
    admm_tau_incr = 1.05,
    admm_tau_decr = 0.98,
)
    if isempty(case_list)
        @info "Empty distributed (AC) case list. No tests will be run."
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

        # ----------------------------
        # Defaults (sempre definidos)
        # ----------------------------
        admm_max_iter_case = 1000
        admm_rho_case = 1.0

        admm_eps_abs_case = admm_eps_abs
        admm_eps_rel_case = admm_eps_rel

        admm_rho_theta_case = nothing
        admm_rho_v_case = nothing
        admm_adaptive_case = nothing

        # defaults mu/tau por caso = os padrões da função (não quebra se não cair em if)
        admm_mu_incr_case = admm_mu_incr
        admm_mu_decr_case = admm_mu_decr
        admm_tau_incr_case = admm_tau_incr
        admm_tau_decr_case = admm_tau_decr

        # defaults de solver por caso = os padrões da função (idem)
        central_optimizer_case = central_optimizer
        distributed_optimizer_case = distributed_optimizer

        # ----------------------------
        # Overrides por caso específico
        # ----------------------------
        if case_name == "pglib_opf_case200_activ.m" ||
           case_name == "pglib_opf_case500_goc.m"

            admm_rho_case = 10000000.0
            admm_max_iter_case = 100000
            admm_eps_abs_case = 1e-4
            admm_eps_rel_case = 1e-4
            admm_adaptive_case = "on"

            admm_mu_incr_case = 10.0
            admm_mu_decr_case = 10.0
            admm_tau_incr_case = 1.2
            admm_tau_decr_case = 0.8

            central_optimizer_case = Ipopt.Optimizer
            distributed_optimizer_case = Ipopt.Optimizer

        elseif case_name == "pglib_opf_case73_ieee_rts.m" ||
               case_name == "pglib_opf_case57_ieee.m" ||
               case_name == "pglib_opf_case118_ieee.m"

            admm_rho_case = 10000000.0
            admm_max_iter_case = 100000
            admm_eps_abs_case = 1e-4
            admm_eps_rel_case = 1e-4
            admm_adaptive_case = "on"

            admm_mu_incr_case = 10.0
            admm_mu_decr_case = 10.0
            admm_tau_incr_case = 1.05
            admm_tau_decr_case = 0.98

            central_optimizer_case = Ipopt.Optimizer
            distributed_optimizer_case = Ipopt.Optimizer
        end

        @info "[$label] Rodando caso $i de $n: $case_name (max_iter=$(admm_max_iter_case), rho=$(admm_rho_case), rho_theta=$(admm_rho_theta_case), rhoV=$(admm_rho_v_case), eps_pr=$(admm_eps_abs_case), eps_dual=$(admm_eps_rel_case), mu=(incr=$(admm_mu_incr_case),decr=$(admm_mu_decr_case)), tau=(incr=$(admm_tau_incr_case),decr=$(admm_tau_decr_case)), adaptive=$(admm_adaptive_case))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_ac_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer_case,
                distributed_optimizer = distributed_optimizer_case,
                rtol = rtol,
                atol = atol,

                admm_rho = admm_rho_case,
                admm_max_iter = admm_max_iter_case,
                admm_eps_abs = admm_eps_abs_case,
                admm_eps_rel = admm_eps_rel_case,

                admm_mu_incr = admm_mu_incr_case,
                admm_mu_decr = admm_mu_decr_case,
                admm_tau_incr = admm_tau_incr_case,
                admm_tau_decr = admm_tau_decr_case,

                admm_rho_theta = admm_rho_theta_case,
                admm_rho_v = admm_rho_v_case,
                admm_adaptive = admm_adaptive_case,
            )
        end
    end

    return true
end
