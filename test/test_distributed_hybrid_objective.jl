# ==============================================
# ADMM-Hybrid-DC vs DC (objetivo) – testes
# ==============================================

# Casos particionados (2/3/4 zonas) — iguais aos do branch
const distributed_hybrid_zone_cases = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
    "pglib_opf_case73_ieee_rts.m",
    "pglib_opf_case57_ieee.m",
    # "pglib_opf_case118_ieee.m",
    "pglib_opf_case200_activ.m",
    # "pglib_opf_case500_goc.m",
]

function track_failure_hybrid(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] =
            Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

"""
Compara objetivo Centralizado DC vs Distribuído ADMM-Hybrid-DC no mesmo caso.
- Central: result_dict["objective"]
- Dist: res.total_real_cost (NamedTuple)
"""
function compare_distributed_hybrid_vs_centralized(
    file::String;
    label::String,
    central_optimizer = get_available_optimizer(OSQP.Optimizer, Ipopt.Optimizer),
    distributed_optimizer = nothing,
    rtol = 0.01,
    atol = 1e-6,
    admm_rho = 1_000.0,
    admm_max_iter = 350,
    admm_eps_abs = 5e-4,
    admm_eps_rel = 5e-4,
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
        Float64(result_dict["objective"])
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
        if distributed_optimizer !== nothing
            push!(args_dist, "--optimizer", string(distributed_optimizer))
        end

        res = Lirion.main(args_dist)
        Float64(res.total_real_cost)   # <- ponto principal
        # se quiser comparar o outro: Float64(res.total_expr_cost)
    catch e
        @warn "Distributed ADMM-Hybrid-DC failed on $case_name: $e"
        track_failure_hybrid(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized = $(round(cent_obj, digits=6)) | distributed-hybrid = $(round(dist_obj, digits=6))"

    abs_err = abs(dist_obj - cent_obj)
    rel_err = abs_err / max(abs(cent_obj), 1.0)
    @info "Errors: abs=$(abs_err) | rel=$(rel_err)"

    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

"""
Roda a suíte de testes para ADMM-Hybrid-DC vs DC nos casos de zonas.
"""
function run_distributed_hybrid_objective_tests(;
    case_list = distributed_hybrid_zone_cases,
    label = "ADMM-Hybrid-DC vs DC",
    central_optimizer = get_available_optimizer(OSQP.Optimizer, Ipopt.Optimizer),
    distributed_optimizer = nothing,
    rtol = 0.01, atol = 1e-6,
    admm_eps_abs = 2e-4, admm_eps_rel = 2e-4,
    admm_alpha::Float64 = 1.0,
)
    if isempty(case_list)
        @info "Empty hybrid case list. No tests will be run."
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
        admm_rho_case = 10000.0

        admm_eps_abs_case = admm_eps_abs
        admm_eps_rel_case = admm_eps_rel
        admm_alpha_case = admm_alpha

        extra_flags_case = String[]
        admm_warm_start_case = false

        # otimizadores por caso (não pisa nas variáveis globais)
        central_optimizer_case = central_optimizer
        distributed_optimizer_case = distributed_optimizer

        # --------------------------------------------
        # Override por caso (ex: pglib_opf_case500_goc.m)
        # --------------------------------------------
        if case_name == "pglib_opf_case500_goc.m"
            admm_rho_case = 5000.0
            central_optimizer_case = Ipopt.Optimizer
            distributed_optimizer_case = Ipopt.Optimizer
        end

        @info "[$label] Rodando caso $i de $n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rho=$(admm_rho_case), " *
              "eps_pr=$(admm_eps_abs_case), eps_dual=$(admm_eps_rel_case), " *
              "alpha=$(admm_alpha_case), warm_start=$(admm_warm_start_case), " *
              "optimizer=$(distributed_optimizer_case === nothing ? "default" : string(distributed_optimizer_case)), " *
              "extra_flags=$(join(extra_flags_case, ' ')))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_hybrid_vs_centralized(
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
                admm_alpha = admm_alpha_case,
                admm_warm_start = admm_warm_start_case,
                extra_flags = extra_flags_case,
            )
        end
    end

    return true
end