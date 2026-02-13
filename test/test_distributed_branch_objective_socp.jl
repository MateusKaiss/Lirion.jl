# ==============================================
# ADMM-(Power coupling)-SOCP vs SOCP (objective) – tests (API solve)
# ==============================================

# -------------------------------------------------
# (Optional) Clarabel robust attributes (as tuple)
# -------------------------------------------------
const CLARABEL_ROBUST_ATTRS = (
    "tol_feas" => 1e-4,
    "tol_gap_abs" => 1e-4,
    "tol_gap_rel" => 1e-4,
    "max_iter" => 1200,
    "dynamic_regularization_enable" => true,
    "dynamic_regularization_eps" => 1e-10,
    "dynamic_regularization_delta" => 1e-4,
    "static_regularization_enable" => true,
    "equilibrate_enable" => true,
    "equilibrate_min_scaling" => 1e-4,
    "equilibrate_max_scaling" => 1e4,
    "direct_solve_method" => :qdldl,
    "iterative_refinement_reltol" => 1e-9,
    "iterative_refinement_abstol" => 1e-9,
    "max_step_fraction" => 0.7,
    "linesearch_backtrack_step" => 0.5,
    "tol_ktratio" => 1e-5,
    "min_terminate_step_length" => 1e-6,
    "min_switch_step_length" => 5e-2,
    "reduced_tol_gap_abs" => 5e-4,
    "reduced_tol_gap_rel" => 5e-4,
    "reduced_tol_feas" => 5e-4,
    "reduced_tol_infeas_abs" => 5e-11,
    "reduced_tol_infeas_rel" => 5e-4,
    "reduced_tol_ktratio" => 1e-4,
)

# -------- helper local (if runtests didn't define) --------
# ENV["LIRION_SOCP_PREFER"] can be "auto" | "gurobi" | "cosmo" | "scs"
if !@isdefined(pick_socp_optimizer)
    const _SOCP_CHAIN_MAP_LOCAL = Dict(
        :auto => (Gurobi.Optimizer, COSMO.Optimizer, SCS.Optimizer),
        :gurobi => (Gurobi.Optimizer, COSMO.Optimizer, SCS.Optimizer),
        :cosmo => (COSMO.Optimizer, Gurobi.Optimizer, SCS.Optimizer),
        :scs => (SCS.Optimizer, Gurobi.Optimizer, COSMO.Optimizer),
    )
    function pick_socp_optimizer(;
        prefer::Symbol = Symbol(lowercase(get(ENV, "LIRION_SOCP_PREFER", "auto"))),
    )
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

# Defaults: uses runtests globals if exist; otherwise local picker
const _CENTRAL_DEFAULT =
    @isdefined(CENTRAL_SOCP_SOLVER) ? CENTRAL_SOCP_SOLVER :
    (
        @isdefined(DEFAULT_SOCP_SOLVER) ? DEFAULT_SOCP_SOLVER :
        pick_socp_optimizer(; prefer = :gurobi)
    )

const _DISTRIB_DEFAULT =
    @isdefined(DISTRIBUTED_SOCP_SOLVER) ? DISTRIBUTED_SOCP_SOLVER :
    (
        @isdefined(DEFAULT_SOCP_SOLVER) ? DEFAULT_SOCP_SOLVER :
        pick_socp_optimizer(; prefer = :auto)
    )

# --- case list ---
const distributed_branch_socp_zone_cases = [
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

function track_failure_branch_socp(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] =
            Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

# ------------------------------------------------------------
# Objective extraction
# ------------------------------------------------------------
@inline function _objective_centralized(out::Lirion.LirionRunResult)
    d = out.result
    d isa AbstractDict || error("Centralized out.result must be Dict, got $(typeof(d))")
    haskey(d, "objective") ||
        error("Centralized result has no key 'objective'. Keys=$(collect(keys(d)))")
    return d["objective"]
end

@inline function _objective_distributed_power_socp(out::Lirion.LirionRunResult)
    r = out.result
    r isa NamedTuple || error("Distributed out.result must be NamedTuple, got $(typeof(r))")
    hasproperty(r, :total_expr_cost) || error(
        "Distributed result has no field :total_expr_cost. Fields=$(propertynames(r))",
    )
    return getproperty(r, :total_expr_cost)
end

function compare_distributed_branch_socp_vs_centralized(
    file::String;
    label::String,
    central_optimizer = _CENTRAL_DEFAULT,
    distributed_optimizer = _DISTRIB_DEFAULT,
    rtol = 0.05,
    atol = 1e-6,

    # ADMM defaults (SOCP + coupling=:power requires ONLY rho_p and rho_q)
    rho_p::Float64 = 1.0e3,
    rho_q::Float64 = 1.0e3,
    admm_max_iter::Int = 350,
    admm_eps_abs::Float64 = 1e-4,
    admm_eps_rel::Float64 = 1e-4,
    admm_adaptive::Union{Nothing, Bool} = nothing,
    admm_mu_incr::Float64 = 3.0,
    admm_mu_decr::Float64 = 3.0,
    admm_tau_incr::Float64 = 1.05,
    admm_tau_decr::Float64 = 0.98,
)
    case_name = basename(file)

    # ================================================================
    # Special case: 2000-bus / 5 zones
    # ================================================================
    if case_name == "pglib_opf_case2000_5zones_goc.m"
        @info "⚙️ [SOCP/:power] Special case $case_name with Clarabel robust (central) and Clarabel rho=1e8 (distributed)."

        cent_obj = try
            opt_cent =
                optimizer_with_attributes(Clarabel.Optimizer, CLARABEL_ROBUST_ATTRS...)
            out_cent = Lirion.solve(
                file;
                model = Lirion.socp(),
                algorithm = Lirion.centralized(),
                optimizer = opt_cent,
                outputs = false,
                plots = false,
            )
            _objective_centralized(out_cent)
        catch e
            @warn "Centralized SOCP failed on $case_name: $e"
            track_failure_branch_socp(label, "Centralized", case_name)
            return
        end

        dist_obj = try
            alg = Lirion.admm(;
                coupling = :power,
                sectioning = :tie_line,
                max_iter = 700,
                rho_p = 1e8,
                rho_q = 1e8,
                eps_abs = 1e-4,
                eps_rel = 1e-4,
                adaptive_rho = true,
                mu_incr = 5.0,
                mu_decr = 5.0,
                tau_incr = 1.1,
                tau_decr = 0.9,
            )
            out_dist = Lirion.solve(
                file;
                model = Lirion.socp(),
                algorithm = alg,
                optimizer = Clarabel.Optimizer,
                outputs = false,
                plots = false,
            )
            _objective_distributed_power_socp(out_dist)
        catch e
            @warn "Distributed ADMM-(power)-SOCP failed on $case_name: $e"
            track_failure_branch_socp(label, "Distributed", case_name)
            return
        end

        @info "File $case_name | centralized SOCP = $(round(cent_obj, digits=6)) | distributed-(power)-SOCP (total_expr_cost) = $(round(dist_obj, digits=6))"
        @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
        return
    end

    # ================================================================
    # General cases
    # ================================================================
    cent_obj = try
        out_cent = Lirion.solve(
            file;
            model = Lirion.socp(),
            algorithm = Lirion.centralized(),
            optimizer = central_optimizer,
            outputs = false,
            plots = false,
        )
        _objective_centralized(out_cent)
    catch e
        @warn "Centralized SOCP failed on $case_name: $e"
        track_failure_branch_socp(label, "Centralized", case_name)
        return
    end

    dist_obj = try
        alg = Lirion.admm(;
            coupling = :power,
            sectioning = :tie_line,
            max_iter = admm_max_iter,
            rho_p = Float64(rho_p),
            rho_q = Float64(rho_q),
            eps_abs = Float64(admm_eps_abs),
            eps_rel = Float64(admm_eps_rel),
            adaptive_rho = isnothing(admm_adaptive) ? true : Bool(admm_adaptive),
            mu_incr = Float64(admm_mu_incr),
            mu_decr = Float64(admm_mu_decr),
            tau_incr = Float64(admm_tau_incr),
            tau_decr = Float64(admm_tau_decr),
        )
        out_dist = Lirion.solve(
            file;
            model = Lirion.socp(),
            algorithm = alg,
            optimizer = distributed_optimizer,
            outputs = false,
            plots = false,
        )
        _objective_distributed_power_socp(out_dist)
    catch e
        @warn "Distributed ADMM-(power)-SOCP failed on $case_name: $e"
        track_failure_branch_socp(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized SOCP = $(round(cent_obj, digits=6)) | distributed-(power)-SOCP (total_expr_cost) = $(round(dist_obj, digits=6)) | " *
          "max_iter=$(admm_max_iter), rhoP=$(rho_p), rhoQ=$(rho_q)"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_branch_socp_objective_tests(;
    case_list = distributed_branch_socp_zone_cases,
    label = "ADMM-(power)-SOCP vs SOCP",
    central_optimizer = _CENTRAL_DEFAULT,
    distributed_optimizer = _DISTRIB_DEFAULT,
    rtol = 0.05,
    atol = 1e-6,
    admm_eps_abs = 1e-4,
    admm_eps_rel = 1e-4,
)
    if isempty(case_list)
        @info "Empty distributed (power-SOCP) case list. No tests will be run."
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

        # defaults for general cases
        rho_p_case = 1e3
        rho_q_case = 1e3
        admm_max_iter_case = 5000
        admm_eps_abs_case = admm_eps_abs
        admm_eps_rel_case = admm_eps_rel

        @info "[$label] Rodando caso $i/$n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rhoP=$(rho_p_case), rhoQ=$(rho_q_case), " *
              "eps_abs=$(admm_eps_abs_case), eps_rel=$(admm_eps_rel_case))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_branch_socp_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer,
                distributed_optimizer = distributed_optimizer,
                rtol = rtol,
                atol = atol,
                rho_p = rho_p_case,
                rho_q = rho_q_case,
                admm_max_iter = admm_max_iter_case,
                admm_eps_abs = admm_eps_abs_case,
                admm_eps_rel = admm_eps_rel_case,
            )
        end
    end
    return true
end
