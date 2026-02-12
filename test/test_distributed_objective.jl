# =============================
# Casos particionados (2/3 zonas)
# =============================
const distributed_zone_cases = [

    # "pglib_opf_case73_ieee_rts__api.m",
    # "pglib_opf_case73_ieee_rts__sad.m",
    ## "pglib_opf_case57_ieee__api.m",
    # "pglib_opf_case57_ieee__sad.m",

    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",

    "pglib_opf_case57_ieee.m",
    "pglib_opf_case73_ieee_rts.m",
    "pglib_opf_case118_ieee.m",
    
    # "pglib_opf_case500_goc.m",

    # "pglib_opf_case200_activ.m",

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
    admm_eps_abs = 1e-6,
    admm_eps_rel = 1e-6,
    admm_alpha::Float64 = 1.0,
    admm_mu_incr::Float64 = 3.0,
    admm_mu_decr::Float64 = 3.0,
    admm_tau_incr::Float64 = 1.05,
    admm_tau_decr::Float64 = 0.98,
    admm_accel::Bool = false,                       # se true, inclui --admm-accel
    admm_adaptive::Union{Nothing, String} = nothing,  # "on" | "off" | nothing
)
    case_name = basename(file)

    # --- Centralizado DC ---
    cent_obj = try
        args_cent = String[
            file,
            "--model", "dc",
            "--optimizer", string(central_optimizer),
        ]
        result_dict, _, _ = Lirion.main(args_cent)
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
            "--admm-eps-abs", string(admm_eps_abs),
            "--admm-eps-rel", string(admm_eps_rel),
            "--admm-alpha", string(admm_alpha),
            "--admm-mu-incr", string(admm_mu_incr),
            "--admm-mu-decr", string(admm_mu_decr),
            "--admm-tau-incr", string(admm_tau_incr),
            "--admm-tau-decr", string(admm_tau_decr),
        ]
        if admm_accel
            push!(args_dist, "--admm-accel")
        end
        if !isnothing(admm_adaptive)
            push!(args_dist, "--admm-adaptive")
            push!(args_dist, admm_adaptive)
        end
        # Só passamos --optimizer se o usuário quiser sobrescrever o default (OSQP)
        if !isnothing(distributed_optimizer)
            push!(args_dist, "--optimizer")
            push!(args_dist, string(distributed_optimizer))
        end

        obj, _ = Lirion.main(args_dist)
        obj
    catch e
        @warn "Distributed ADMM-DC failed on $case_name: $e"
        track_failure_dist(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized = $(round(cent_obj, digits=6)) | distributed = $(round(dist_obj, digits=6))"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

# -----------------------------------------------------------
# run_distributed_objective_tests: overrides p/ 30000 + novos params
# -----------------------------------------------------------
function run_distributed_objective_tests(;
    case_list = distributed_zone_cases,
    label = "ADMM-DC (bus) vs DC",
    central_optimizer = get_available_optimizer(OSQP.Optimizer, Ipopt.Optimizer),
    distributed_optimizer::Union{Nothing, Function} = nothing,
    rtol = 0.05, atol = 1e-6,
    admm_eps_abs = 1e-4, admm_eps_rel = 1e-4,
    admm_alpha::Float64 = 1.0,  # over-relaxation inativo
    admm_accel::Bool = false,
    admm_adaptive::Union{Nothing, String} = nothing,
)
    if isempty(case_list)
        @info "Empty distributed case list. No tests will be run."
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

        # defaults por zonas
        admm_max_iter_case, admm_rho_case =  (5000, 1000.0)


        admm_eps_abs_case = admm_eps_abs
        admm_eps_rel_case = admm_eps_rel
        admm_alpha_case = admm_alpha

        # parâmetros adaptativos: defaults gerais
        admm_mu_incr_case = 3.0
        admm_mu_decr_case = 3.0
        admm_tau_incr_case = 1.05
        admm_tau_decr_case = 0.98
        admm_accel_case = admm_accel
        admm_adaptive_case = admm_adaptive  # nothing | "on" | "off"

        if case_name == "pglib_opf_case73_ieee_rts.m" ||
        case_name == "pglib_opf_case57_ieee.m" ||
        case_name == "pglib_opf_case118_ieee.m" ||
        case_name == "pglib_opf_case200_activ.m"

            admm_rho_case = 100000.0
            admm_max_iter_case = 100000
            admm_eps_abs_case = 1e-4
            admm_eps_rel_case = 1e-4
            admm_adaptive_case = "off"

            central_optimizer_case = Ipopt.Optimizer
            distributed_optimizer_case = Ipopt.Optimizer
        end
        central_optimizer_case = Ipopt.Optimizer
        distributed_optimizer_case = Ipopt.Optimizer

        @info "[$label] Rodando caso $i de $n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rho=$(admm_rho_case), " *
              "eps_pr=$(admm_eps_abs_case), eps_dual=$(admm_eps_rel_case), " *
              "alpha=$(admm_alpha_case), mu_incr=$(admm_mu_incr_case), mu_decr=$(admm_mu_decr_case), " *
              "tau_incr=$(admm_tau_incr_case), tau_decr=$(admm_tau_decr_case), " *
              "accel=$(admm_accel_case), adaptive=$(admm_adaptive_case), " *
              "CENTRAL_SOLVER=$(typeof(central_optimizer_case)))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer_case,
                distributed_optimizer = distributed_optimizer,
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
                admm_accel = admm_accel_case,
                admm_adaptive = admm_adaptive_case,
            )
        end
    end
    return true
end
