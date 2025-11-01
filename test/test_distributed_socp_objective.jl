# ============================================
# ADMM-SOCP vs SOCP (objetivo) – testes
# ============================================
using Test
using Lirion
using Glob
using SCS
using Gurobi
using COSMO

const ADMM_MAX_IT = 500

# -------- helper local (caso o runtests não tenha definido) --------
# ENV["LIRION_SOCP_PREFER"] pode ser "auto" | "gurobi" | "cosmo" | "scs"
if !@isdefined(pick_socp_optimizer)
    const _SOCP_CHAIN_MAP_LOCAL = Dict(
        :auto   => (Gurobi.Optimizer, COSMO.Optimizer, SCS.Optimizer),
        :gurobi => (Gurobi.Optimizer, COSMO.Optimizer, SCS.Optimizer),
        :cosmo  => (COSMO.Optimizer,  Gurobi.Optimizer, SCS.Optimizer),
        :scs    => (SCS.Optimizer,    Gurobi.Optimizer, COSMO.Optimizer),
    )
    function pick_socp_optimizer(; prefer::Symbol = Symbol(
        lowercase(get(ENV, "LIRION_SOCP_PREFER", "auto"))
    ))
        chain = get(_SOCP_CHAIN_MAP_LOCAL, prefer, _SOCP_CHAIN_MAP_LOCAL[:auto])
        for factory in chain
            try
                _ = factory()
                return factory
            catch
            end
        end
        return SCS.Optimizer
    end
end

# Defaults: usa os do runtests se existirem; senão, cai para o picker local
const _CENTRAL_DEFAULT = @isdefined(CENTRAL_SOCP_SOLVER)     ? CENTRAL_SOCP_SOLVER     :
                         (@isdefined(DEFAULT_SOCP_SOLVER)    ? DEFAULT_SOCP_SOLVER     :
                          pick_socp_optimizer(; prefer = :gurobi))

const _DISTRIB_DEFAULT = @isdefined(DISTRIBUTED_SOCP_SOLVER) ? DISTRIBUTED_SOCP_SOLVER :
                         (@isdefined(DEFAULT_SOCP_SOLVER)    ? DEFAULT_SOCP_SOLVER     :
                          pick_socp_optimizer(; prefer = :auto))

# --- lista de casos particionados (2/3/4 zonas) ---
const distributed_socp_zone_cases = [
    # 14 buses
    "pglib_opf_case14_2zones_ieee_socp.m",
    "pglib_opf_case14_3zones_ieee_socp.m",
    # 24 buses
    "pglib_opf_case24_2zones_ieee_rts.m",
    "pglib_opf_case24_3zones_ieee_rts.m",
    # 30 buses
    "pglib_opf_case30_2zones_as.m",
    "pglib_opf_case30_3zones_as.m",
    # 57/118/500 (descomentando conforme seus arquivos)
    # "pglib_opf_case57_2zones_ieee.m",
    # "pglib_opf_case57_4zones_ieee.m",
    # "pglib_opf_case118_2zones_ieee.m",
    # "pglib_opf_case500_2zones_goc.m",
]

function track_failure_socp(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] =
            Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

function compare_distributed_socp_vs_centralized(
    file::String;
    label::String,
    central_optimizer = _CENTRAL_DEFAULT,     # agora com fallback (Gurobi→COSMO→SCS)
    distributed_optimizer = _DISTRIB_DEFAULT, # idem para o distribuído
    rtol = 0.05,
    atol = 1e-6,
    # ADMM penalties
    rho_v::Float64 = 1.0e3,
    rho_flow_p::Float64 = 3.0e2,
    rho_flow_q::Float64 = 3.0e2,
    admm_eps_primal::Float64 = 1e-4,
    admm_eps_dual::Float64 = 1e-4,
)
    case_name = basename(file)

    # --- Centralizado SOCP ---
    cent_obj = try
        args_cent = String[
            file,
            "--model", "socp",
            "--optimizer", string(central_optimizer),
        ]
        result_dict, _, _ = Lirion.main(args_cent)
        result_dict["objective"]
    catch e
        @warn "Centralized SOCP failed on $case_name: $e"
        track_failure_socp(label, "Centralized", case_name)
        return
    end

    # --- Distribuído ADMM-SOCP (travado em 300) ---
    dist_obj = try
        args_dist = String[
            file,
            "--model", "distributed-admm-socp",
            "--optimizer", string(distributed_optimizer),
            "--admm-rho-v", string(rho_v),
            "--admm-rho-flow-p", string(rho_flow_p),
            "--admm-rho-flow-q", string(rho_flow_q),
            "--admm-max-iter", string(ADMM_MAX_IT),           # <- FIXO 300
            "--admm-eps-primal", string(admm_eps_primal),
            "--admm-eps-dual", string(admm_eps_dual),
        ]
        obj, _ = Lirion.main(args_dist)
        obj
    catch e
        @warn "Distributed ADMM-SOCP failed on $case_name: $e"
        track_failure_socp(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized SOCP = $(round(cent_obj, digits=6)) | distributed-ADMM-SOCP = $(round(dist_obj, digits=6)) | ADMM_MAX_IT=$(ADMM_MAX_IT)"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_socp_objective_tests(;
    case_list = distributed_socp_zone_cases,
    label = "ADMM-SOCP vs SOCP",
    central_optimizer = _CENTRAL_DEFAULT,
    distributed_optimizer = _DISTRIB_DEFAULT,
    rtol = 0.05, atol = 1e-6,
    admm_eps_primal = 1e-4, admm_eps_dual = 1e-4,
)
    if isempty(case_list)
        @info "Empty SOCP distributed case list. No tests will be run."
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

        # penalidades constantes
        rho_v_case = 1.0e3
        rho_flow_p_case = 3.0e2
        rho_flow_q_case = 3.0e2

        admm_max_iter_case = ADMM_MAX_IT
        admm_eps_primal_case = admm_eps_primal
        admm_eps_dual_case = admm_eps_dual

        @info "[$label] Rodando caso $i de $n: $case_name " *
              "(max_iter=$(admm_max_iter_case), ρV²=$(rho_v_case), ρP=$(rho_flow_p_case), ρQ=$(rho_flow_q_case), " *
              "eps_pr=$(admm_eps_primal_case), eps_dual=$(admm_eps_dual_case))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_socp_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer,
                distributed_optimizer = distributed_optimizer,
                rtol = rtol,
                atol = atol,
                rho_v = rho_v_case,
                rho_flow_p = rho_flow_p_case,
                rho_flow_q = rho_flow_q_case,
                admm_eps_primal = admm_eps_primal_case,
                admm_eps_dual = admm_eps_dual_case,
            )
        end
    end
    return true
end