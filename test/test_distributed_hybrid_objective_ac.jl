const distributed_hybrid_ac_zone_cases = [
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

function track_failure_hybrid_ac(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] =
            Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

function compare_distributed_hybrid_ac_vs_centralized(
    file::String;
    label::String,
    central_optimizer = Ipopt.Optimizer,
    distributed_optimizer = Ipopt.Optimizer,
    rtol = 0.01,
    atol = 1e-6,
    rho_base::Union{Nothing, Float64} = 1.0e3,
    rho_theta::Union{Nothing, Float64} = nothing,
    rho_v::Union{Nothing, Float64} = nothing,
    rho_flow_p::Union{Nothing, Float64} = 1.0e5,
    rho_flow_q::Union{Nothing, Float64} = 1.0e5,
    admm_max_iter::Int = 500,
    admm_eps_abs::Float64 = 1e-7,
    admm_eps_rel::Float64 = 1e-7,
    admm_mu_incr::Float64 = 2.0,
    admm_mu_decr::Float64 = 2.0,
    admm_tau_incr::Float64 = 1.05,
    admm_tau_decr::Float64 = 0.98,
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
        track_failure_hybrid_ac(label, "Centralized", case_name)
        return
    end

    # --- Distribuído ADMM-Hybrid-AC ---
    dist_obj = try
        args_dist = String[
            file,
            "--model", "distributed-admm-hybrid-ac",
            "--optimizer", string(distributed_optimizer),
            "--admm-max-iter", string(admm_max_iter),
            "--admm-eps-abs", string(admm_eps_abs),
            "--admm-eps-rel", string(admm_eps_rel),
            "--admm-mu-incr", string(admm_mu_incr),
            "--admm-mu-decr", string(admm_mu_decr),
            "--admm-tau-incr", string(admm_tau_incr),
            "--admm-tau-decr", string(admm_tau_decr),
        ]

        if rho_base !== nothing
            push!(args_dist, "--admm-rho")
            push!(args_dist, string(rho_base))
        end
        if rho_theta !== nothing
            push!(args_dist, "--admm-rho-theta")
            push!(args_dist, string(rho_theta))
        end
        if rho_v !== nothing
            push!(args_dist, "--admm-rho-v")
            push!(args_dist, string(rho_v))
        end
        if rho_flow_p !== nothing
            push!(args_dist, "--admm-rho-flow-p")
            push!(args_dist, string(rho_flow_p))
        end
        if rho_flow_q !== nothing
            push!(args_dist, "--admm-rho-flow-q")
            push!(args_dist, string(rho_flow_q))
        end
        if admm_adaptive !== nothing
            push!(args_dist, "--admm-adaptive")
            push!(args_dist, admm_adaptive)
        end

        obj, _ = Lirion.main(args_dist)
        obj
    catch e
        @warn "Distributed ADMM-Hybrid-AC failed on $case_name: $e"
        track_failure_hybrid_ac(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized AC = $(round(cent_obj, digits=6)) | " *
          "distributed-hybrid-AC = $(round(dist_obj, digits=6))"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_hybrid_ac_objective_tests(;
    case_list = distributed_hybrid_ac_zone_cases,
    label = "ADMM-Hybrid-AC vs AC",
    central_optimizer = Ipopt.Optimizer,
    distributed_optimizer = Ipopt.Optimizer,
    rtol = 0.01, atol = 1e-6,
    admm_eps_abs = 1e-4, admm_eps_rel = 1e-4,
    admm_mu_incr::Float64 = 2.0,
    admm_mu_decr::Float64 = 2.0,
    admm_tau_incr::Float64 = 1.05,
    admm_tau_decr::Float64 = 0.98,
)
    if isempty(case_list)
        @info "Empty hybrid-AC case list. No tests will be run."
        return true
    end

    test_dir = joinpath(@__DIR__, "distributed")

    selected_files = [
        joinpath(test_dir, fname)
        for fname in case_list
        if isfile(joinpath(test_dir, fname))
    ]

    @info "hybrid-ac: test_dir = $test_dir"
    @info "hybrid-ac: selected_files = $(length(selected_files))"
    isempty(selected_files) && error("Nenhum arquivo encontrado. Confira nomes e pasta.")

    n = length(selected_files)
    for (i, file) in enumerate(selected_files)
        case_name = basename(file)

        # defaults por padrão (2/3 zonas)
        admm_max_iter_case, rho_bus_base, rho_flow_base =
            if occursin("_3zones_", case_name)
                (400, 1.0e3, 1.0e5)
            elseif occursin("_2zones_", case_name)
                (350, 1.0e3, 1.0e5)
            else
                (300, 1.0e3, 1.0e5)
            end

        admm_eps_abs_case = admm_eps_abs
        admm_eps_rel_case = admm_eps_rel
        admm_mu_incr_case = admm_mu_incr
        admm_mu_decr_case = admm_mu_decr
        admm_tau_incr_case = admm_tau_incr
        admm_tau_decr_case = admm_tau_decr
        admm_adaptive_case = nothing
        if case_name == "pglib_opf_case118_2zones_ieee.m" ||
           case_name == "pglib_opf_case118_3zones_ieee.m" ||
           case_name == "pglib_opf_case118_5zones_ieee.m"
            admm_max_iter_case = 500
            admm_eps_abs_case = 1e-4
            admm_eps_rel_case = 1e-4
            rho_bus_base = 1.0e5     # --admm-rho-theta / --admm-rho-v
            rho_flow_base = 1.0e5     # --admm-rho-flow-p / --admm-rho-flow-q
            admm_adaptive_case = "on"      # --admm-adaptive on

            admm_mu_incr_case = 2.0
            admm_mu_decr_case = 2.0
            admm_tau_incr_case = 1.05
            admm_tau_decr_case = 0.95
        end

        if case_name == "pglib_opf_case500_3zones_ieee.m" ||
           case_name == "pglib_opf_case500_4zones_ieee.m"
            admm_max_iter_case = 500
            admm_eps_abs_case = 1e-4
            admm_eps_rel_case = 1e-4
            rho_bus_base = 1e4     # --admm-rho-theta / --admm-rho-v
            rho_flow_base = 1e4     # --admm-rho-flow-p / --admm-rho-flow-q
            admm_adaptive_case = "on"      # --admm-adaptive on
            admm_mu_incr_case = 3.0
            admm_mu_decr_case = 3.0
            admm_tau_incr_case = 1.1
            admm_tau_decr_case = 0.9
        end

        if case_name == "pglib_opf_case2000_5zones_goc.m"
            admm_max_iter_case = 500
            admm_eps_abs_case = 1e-4
            admm_eps_rel_case = 1e-4
            rho_bus_base = 1.0e4   # --admm-rho-theta / --admm-rho-v
            rho_flow_base = 1.0e4   # --admm-rho-flow-p / --admm-rho-flow-q
            admm_adaptive_case = "on"
            admm_mu_incr_case, admm_mu_decr_case = 2.0, 2.0
            admm_tau_incr_case, admm_tau_decr_case = 1.1, 0.9
        end

        if case_name == "pglib_opf_case3970_25zones_goc.m"
            admm_max_iter_case = 500
            admm_eps_abs_case = 3e-4
            admm_eps_rel_case = 3e-4
            rho_bus_base = 1.0e4   # --admm-rho-theta / --admm-rho-v
            rho_flow_base = 1.0e4   # --admm-rho-flow-p / --admm-rho-flow-q
            admm_adaptive_case = "on"
            admm_mu_incr_case, admm_mu_decr_case = 2.0, 2.0
            admm_tau_incr_case, admm_tau_decr_case = 1.3, 0.7
        end

        if case_name == "pglib_opf_case10000_10zones_ac_goc.m"
            admm_max_iter_case = 800
            admm_eps_abs_case = 3e-4
            admm_eps_rel_case = 3e-4
            rho_bus_base = 1.0e3   # --admm-rho-theta / --admm-rho-v
            rho_flow_base = 1.0e3   # --admm-rho-flow-p / --admm-rho-flow-q
            admm_adaptive_case = "on"
            admm_mu_incr_case, admm_mu_decr_case = 3.0, 3.0
            admm_tau_incr_case, admm_tau_decr_case = 1.05, 0.95
        end

        @info "[$label] Rodando caso $i de $n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rho_theta/rhoV=$(rho_bus_base), rho_flow=$(rho_flow_base), " *
              "eps_pr=$(admm_eps_abs_case), eps_dual=$(admm_eps_rel_case), " *
              "mu=(incr=$(admm_mu_incr_case),decr=$(admm_mu_decr_case)), " *
              "tau=(incr=$(admm_tau_incr_case),decr=$(admm_tau_decr_case)), " *
              "adaptive=$(admm_adaptive_case))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_hybrid_ac_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer,
                distributed_optimizer = distributed_optimizer,
                rtol = rtol,
                atol = atol,
                rho_base = nothing,                 # usar rhos separados
                rho_theta = rho_bus_base,
                rho_v = rho_bus_base,
                rho_flow_p = rho_flow_base,
                rho_flow_q = rho_flow_base,
                admm_max_iter = admm_max_iter_case,
                admm_eps_abs = admm_eps_abs_case,
                admm_eps_rel = admm_eps_rel_case,
                admm_mu_incr = admm_mu_incr_case,
                admm_mu_decr = admm_mu_decr_case,
                admm_tau_incr = admm_tau_incr_case,
                admm_tau_decr = admm_tau_decr_case,
                admm_adaptive = admm_adaptive_case,
            )
        end
    end
    return true
end
