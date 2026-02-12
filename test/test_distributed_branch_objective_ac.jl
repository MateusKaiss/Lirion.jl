# ==============================================
# ADMM-Branch-AC vs AC (objetivo) – testes
# ==============================================
const distributed_branch_ac_zone_cases = [
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
    rtol = 0.01,
    atol = 1e-6,
    admm_rho_flow = 1.0e5,
    admm_rho_flow_p = nothing,
    admm_rho_flow_q = nothing,
    admm_max_iter = 5000,
    admm_eps_abs = 1e-4,
    admm_eps_rel = 1e-4,
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
            "--admm-eps-abs", string(admm_eps_abs),
            "--admm-eps-rel", string(admm_eps_rel),
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
    rtol = 0.01, atol = 1e-6,
    admm_eps_abs = 1e-4, admm_eps_rel = 1e-4,
    admm_mu_incr = 3.0,
    admm_mu_decr = 3.0,
    admm_tau_incr = 1.05,
    admm_tau_decr = 0.98,
)
    if isempty(case_list)
        @info "Empty distributed (branch-AC) case list. No tests will be run."
        return true
    end

    test_dir = joinpath(@__DIR__, "distributed")
    selected_files =
        [joinpath(test_dir, f) for f in case_list if isfile(joinpath(test_dir, f))]

    n = length(selected_files)
    for (i, file) in enumerate(selected_files)
        case_name = basename(file)

        admm_max_iter_case, admm_rho_case = (5000, 10000.0)

        admm_eps_abs_case = admm_eps_abs
        admm_eps_rel_case = admm_eps_rel
        mu_incr_case = admm_mu_incr
        mu_decr_case = admm_mu_decr
        tau_incr_case = admm_tau_incr
        tau_decr_case = admm_tau_decr
        adaptive_case = nothing
        extra_flags = String[]

        @info "[$label] Rodando caso $i de $n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rho_flow=$(admm_rho_case), " *
              "eps_pr=$(admm_eps_abs_case), eps_dual=$(admm_eps_rel_case), " *
              "mu=(incr=$(mu_incr_case),decr=$(mu_decr_case)), " *
              "tau=(incr=$(tau_incr_case),decr=$(tau_decr_case)), " *
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
                admm_eps_abs = admm_eps_abs_case,
                admm_eps_rel = admm_eps_rel_case,
                admm_mu_incr = mu_incr_case,
                admm_mu_decr = mu_decr_case,
                admm_tau_incr = tau_incr_case,
                admm_tau_decr = tau_decr_case,
                admm_adaptive = adaptive_case,
                additional_flags = extra_flags,
            )
        end
    end
    return true
end
