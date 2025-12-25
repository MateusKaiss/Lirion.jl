# Casos particionados (2/3/4 zonas) — iguais aos do branch
const distributed_hybrid_zone_cases = [
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
    # # 118 buses
    # "pglib_opf_case118_2zones_ieee.m",
    # "pglib_opf_case118_3zones_ieee.m",
    # "pglib_opf_case118_5zones_ieee.m",
    # # 500 buses
    # "pglib_opf_case500_3zones_goc.m",
    # "pglib_opf_case500_4zones_goc.m",
    # # 793
    # "pglib_opf_case793_11zones_goc.m",
    # # 2000
    # "pglib_opf_case2000_5zones_goc.m",
    # "pglib_opf_case3970_25zones_goc.m",
    # "pglib_opf_case10000_50zones_goc.m",
    # "pglib_opf_case20758_42zones_epigrids.m",
    # "pglib_opf_case30000_20zones_goc.m",
    # "pglib_opf_case78484_157zones_epigrids.m",
]

function track_failure_hybrid(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] =
            Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

function compare_distributed_hybrid_vs_centralized(
    file::String;
    label::String,
    central_optimizer = get_available_optimizer(OSQP.Optimizer, Ipopt.Optimizer),
    distributed_optimizer::Union{Nothing, Function} = nothing,
    rtol = 0.05,
    atol = 1e-6,
    admm_rho = 1_000.0,
    admm_max_iter = 350,
    admm_eps_abs = 1e-4,
    admm_eps_rel = 1e-4,
    admm_alpha::Float64 = 1.0,
    admm_warm_start::Bool = false,
    extra_flags::Vector{String} = String[],
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
        track_failure_hybrid(label, "Centralized", case_name)
        return
    end

    # --- Distribuído ADMM-Hybrid-DC ---
    dist_obj = try
        args_dist = String[
            file,
            "--model", "distributed-admm-hybrid-dc",
            "--admm-rho", string(admm_rho),
            "--admm-max-iter", string(admm_max_iter),
            "--admm-eps-abs", string(admm_eps_abs),
            "--admm-eps-rel", string(admm_eps_rel),
            "--admm-alpha", string(admm_alpha),
        ]
        if !isempty(extra_flags)
            append!(args_dist, extra_flags)
        end
        if admm_warm_start
            push!(args_dist, "--admm-warm-start")
        end
        if !isnothing(distributed_optimizer)
            push!(args_dist, "--optimizer")
            push!(args_dist, string(distributed_optimizer))
        end
        obj, _ = Lirion.main(args_dist)
        obj
    catch e
        @warn "Distributed ADMM-Hybrid-DC failed on $case_name: $e"
        track_failure_hybrid(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized = $(round(cent_obj, digits=6)) | distributed-hybrid = $(round(dist_obj, digits=6))"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_hybrid_objective_tests(;
    case_list = distributed_hybrid_zone_cases,
    label = "ADMM-Hybrid-DC vs DC",
    central_optimizer = get_available_optimizer(OSQP.Optimizer, Ipopt.Optimizer),
    distributed_optimizer::Union{Nothing, Function} = nothing,
    rtol = 0.05, atol = 1e-6,
    admm_eps_abs = 1e-4, admm_eps_rel = 1e-4,
    admm_alpha::Float64 = 1.0,
)
    if isempty(case_list)
        @info "Empty hybrid case list. No tests will be run."
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
                (300, 1_000.0)
            elseif occursin("_5zones_", case_name)
                (700, 1_000.0)
            else
                (700, 1_000.0)
            end

        admm_eps_abs_case = admm_eps_abs
        admm_eps_rel_case = admm_eps_rel
        admm_alpha_case = admm_alpha

        extra_flags_case = String[]
        admm_warm_start_case = false

        if occursin("case2000_5zones_", case_name)
            admm_max_iter_case = 600
            admm_alpha_case = 1.2
            admm_rho_case = 1.0
            extra_flags_case = [
                "--admm-rho-theta", "700",
                "--admm-rho-flow", "300",
            ]
        end

        if occursin("case3970_25zones_", case_name)
            admm_max_iter_case = 600
            admm_alpha_case = 1.2
            admm_rho_case = 1.0
            extra_flags_case = [
                "--admm-rho-theta", "700",
                "--admm-rho-flow", "300",
            ]
        end

        if occursin("case10000_50zones_", case_name)
            admm_max_iter_case = 600
            admm_alpha_case = 1.2
            admm_rho_case = 1.0
            extra_flags_case = [
                "--admm-rho-theta", "700",
                "--admm-rho-flow", "300",
            ]
        end

        if occursin("case20758_42zones_", case_name)
            admm_max_iter_case = 600
            admm_alpha_case = 1.1
            admm_rho_case = 1.0
            admm_warm_start_case = true
            extra_flags_case = [
                "--admm-rho-theta", "700",
                "--admm-rho-flow", "300",
            ]
        end

        if occursin("case30000_20zones_", case_name)
            admm_max_iter_case = 1500
            admm_alpha_case = 1.0
            admm_rho_case = 1.0
            admm_eps_abs_case = 3e-4
            admm_eps_rel_case = 3e-4
            admm_warm_start_case = false
            extra_flags_case = String[
                "--admm-adaptive", "on",
                "--admm-rho-theta", "1000.0",
                "--admm-rho-flow", "1000.0",
                "--admm-mu-incr", "3.0",
                "--admm-mu-decr", "3.0",
                "--admm-tau-incr", "1.050",
                "--admm-tau-decr", "0.92",
                "--optimizer", "OSQP.Optimizer",
            ]
        end

        if occursin("case78484_157zones_", case_name)
            admm_max_iter_case = 2000
            admm_alpha_case = 1.0
            admm_rho_case = 1.0
            admm_eps_abs_case = 5e-4
            admm_eps_rel_case = 5e-4
            admm_warm_start_case = false
            extra_flags_case = String[
                "--admm-mu-incr", "3.0",
                "--admm-mu-decr", "3.0",
                "--admm-tau-incr", "1.05",
                "--admm-tau-decr", "0.94",
                "--admm-rho-theta", "50000",
                "--admm-rho-flow", "50000",
            ]
        end

        @info "[$label] Rodando caso $i de $n: $case_name (max_iter=$(admm_max_iter_case), rho=$(admm_rho_case), eps_pr=$(admm_eps_abs_case), eps_dual=$(admm_eps_rel_case), alpha=$(admm_alpha_case), warm_start=$(admm_warm_start_case), extra_flags=$(join(extra_flags_case, ' ')))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_hybrid_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer,
                distributed_optimizer = distributed_optimizer,
                rtol = rtol,
                atol = atol,
                admm_rho = admm_rho_case,
                admm_max_iter = admm_max_iter_case,
                admm_eps_abs = admm_eps_abs_case,
                admm_eps_rel = admm_eps_rel_case,
                admm_alpha = admm_alpha_case,
                admm_warm_start = admm_warm_start_case,
                extra_flags = extra_flags_case,
            )
        end
    end
    return true
end
