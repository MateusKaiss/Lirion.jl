const distributed_hybrid_ac_zone_cases = [
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
    # 118/500 buses (descomente quando quiser)
    # "pglib_opf_case118_2zones_ieee.m",
    # "pglib_opf_case500_2zones_goc.m",
]

function track_failure_hybrid_ac(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] = Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

function compare_distributed_hybrid_ac_vs_centralized(
    file::String;
    label::String,
    central_optimizer = Ipopt.Optimizer,
    distributed_optimizer = Ipopt.Optimizer,
    rtol = 0.05,
    atol = 1e-6,
    rho_base::Union{Nothing,Float64} = 1.0e3,
    rho_theta::Union{Nothing,Float64} = nothing,
    rho_v::Union{Nothing,Float64} = nothing,
    rho_flow_p::Union{Nothing,Float64} = 1.0e5,
    rho_flow_q::Union{Nothing,Float64} = 1.0e5,
    admm_max_iter::Int = 500,
    admm_eps_primal::Float64 = 1e-7,
    admm_eps_dual::Float64 = 1e-7,
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
            "--admm-eps-primal", string(admm_eps_primal),
            "--admm-eps-dual", string(admm_eps_dual),
        ]

        if rho_base !== nothing
            push!(args_dist, "--admm-rho")
            push!(args_dist, string(rho_base))
        end
        if rho_theta !== nothing
            push!(args_dist, "--admm-rho-theta"); push!(args_dist, string(rho_theta))
        end
        if rho_v !== nothing
            push!(args_dist, "--admm-rho-v");     push!(args_dist, string(rho_v))
        end
        if rho_flow_p !== nothing
            push!(args_dist, "--admm-rho-flow-p"); push!(args_dist, string(rho_flow_p))
        end
        if rho_flow_q !== nothing
            push!(args_dist, "--admm-rho-flow-q"); push!(args_dist, string(rho_flow_q))
        end

        obj, _ = PowerOptimization.main(args_dist)
        obj
    catch e
        @warn "Distributed ADMM-Hybrid-AC failed on $case_name: $e"
        track_failure_hybrid_ac(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized AC = $(round(cent_obj, digits=6)) | distributed-hybrid-AC = $(round(dist_obj, digits=6))"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_hybrid_ac_objective_tests(;
    case_list = distributed_hybrid_ac_zone_cases,
    label = "ADMM-Hybrid-AC vs AC",
    central_optimizer = Ipopt.Optimizer,
    distributed_optimizer = Ipopt.Optimizer,
    rtol = 0.05, atol = 1e-6,
    admm_eps_primal = 1e-8, admm_eps_dual = 1e-8,
)
    if isempty(case_list)
        @info "Empty hybrid-AC case list. No tests will be run."
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

        admm_max_iter_case, rho_bus_base, rho_flow_base =
            if occursin("_3zones_", case_name)
                (600, 1.0e3, 1.0e5)
            elseif occursin("_2zones_", case_name)
                (500, 1.0e3, 1.0e5)
            else
                (500, 1.0e3, 1.0e5)
            end

        admm_eps_primal_case = admm_eps_primal
        admm_eps_dual_case   = admm_eps_dual

        # ---- overrides para casos maiores (descomente quando ativar) ----
        # if occursin("case118_", case_name)
        #     admm_max_iter_case = max(admm_max_iter_case, 700)
        #     rho_bus_base  = 3.0e3
        #     rho_flow_base = 3.0e5
        #     admm_eps_primal_case = min(admm_eps_primal_case, 1e-8)
        #     admm_eps_dual_case   = min(admm_eps_dual_case, 1e-8)
        # end
        # if occursin("case500_", case_name)
        #     admm_max_iter_case = 800
        #     rho_bus_base  = 5.0e3
        #     rho_flow_base = 5.0e5
        #     admm_eps_primal_case = 1e-8
        #     admm_eps_dual_case   = 1e-8
        # end

        @info "[$label] Rodando caso $i de $n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rho_bus=$(rho_bus_base), rho_flow=$(rho_flow_base), " *
              "eps_pr=$(admm_eps_primal_case), eps_dual=$(admm_eps_dual_case))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_hybrid_ac_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer,
                distributed_optimizer = distributed_optimizer,
                rtol = rtol,
                atol = atol,
                rho_base = nothing, 
                rho_theta = rho_bus_base,
                rho_v     = rho_bus_base,
                rho_flow_p = rho_flow_base,
                rho_flow_q = rho_flow_base,
                admm_max_iter = admm_max_iter_case,
                admm_eps_primal = admm_eps_primal_case,
                admm_eps_dual = admm_eps_dual_case,
            )
        end
    end
    return true
end