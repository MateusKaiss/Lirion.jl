# ==============================================
# ADMM-Branch-DC vs DC (objective) – tests (API solve)
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
    # "pglib_opf_case57_ieee.m",
    # "pglib_opf_case118_ieee.m",
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

# ------------------------------------------------------------
# Objective extraction (current return formats)
# ------------------------------------------------------------
@inline function _objective_centralized_dc(out::Lirion.LirionRunResult)
    d = out.result
    d isa AbstractDict || error("Centralized out.result must be Dict, got $(typeof(d))")
    haskey(d, "objective") || error("Centralized result has no key 'objective'. Keys=$(collect(keys(d)))")
    return Float64(d["objective"])
end

@inline function _objective_distributed_branch_dc(out::Lirion.LirionRunResult)
    r = out.result
    r isa NamedTuple || error("Distributed out.result must be NamedTuple, got $(typeof(r))")

    # You used total_real_cost before; keep it as default.
    hasproperty(r, :total_real_cost) ||
        error("Distributed result has no field :total_real_cost. Fields=$(propertynames(r))")

    return Float64(getproperty(r, :total_real_cost))

    # If you ever want to compare expr cost instead:
    # return Float64(getproperty(r, :total_expr_cost))
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
    # Centralized DC -> objective
    # ----------------------------
    cent_obj = try
        out_cent = Lirion.solve(
            file;
            model = Lirion.dc(),
            algorithm = Lirion.centralized(),
            optimizer = central_optimizer,
            outputs = false,
            plots = false,
        )
        _objective_centralized_dc(out_cent)
    catch e
        @warn "Centralized DC failed on $case_name: $e"
        track_failure_branch(label, "Centralized", case_name)
        return
    end

    # --------------------------------------
    # Distributed ADMM-Branch-DC -> NamedTuple
    # --------------------------------------
    dist_obj = try
        alg = Lirion.admm_branch(
            max_iter = admm_max_iter,
            rho = Float64(admm_rho),
            eps_abs = Float64(admm_eps_abs),
            eps_rel = Float64(admm_eps_rel),
            # NOTE: your config uses alpha_overrelax (not alpha)
            alpha_overrelax = (admm_alpha == 1.0 ? nothing : Float64(admm_alpha)),
            warm_start = Bool(admm_warm_start),
            adaptive_rho = true,  # keep old behavior unless you add a knob here
            mu_incr = Float64(admm_mu_incr),
            mu_decr = Float64(admm_mu_decr),
            tau_incr = Float64(admm_tau_incr),
            tau_decr = Float64(admm_tau_decr),
        )

        out_dist = Lirion.solve(
            file;
            model = Lirion.dc(),
            algorithm = alg,
            optimizer = distributed_optimizer, # can be nothing; solve() will keep it if cfg has it
            outputs = false,
            plots = false,
        )

        _objective_distributed_branch_dc(out_dist)
    catch e
        @warn "Distributed ADMM-Branch-DC failed on $case_name: $e"
        track_failure_branch(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized = $(round(cent_obj, digits=6)) | distributed-branch (total_real_cost) = $(round(dist_obj, digits=6))"

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
    rtol = 0.01,
    atol = 1e-6,
    admm_eps_abs = 1e-4,
    admm_eps_rel = 1e-4,
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
        # Override por caso (ex: pglib_opf_case500_goc)
        # --------------------------------------------
        if case_name == "pglib_opf_case500_goc.m"
            admm_rho_case = 5000.0
            distributed_optimizer_case = Ipopt.Optimizer
            # (se você quiser trocar o central tbm, crie central_optimizer_case e passe abaixo)
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
