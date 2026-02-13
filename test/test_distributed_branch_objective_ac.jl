# ==============================================
# ADMM-(Power coupling)-AC vs AC (objective) – tests (API solve)
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

# ------------------------------------------------------------
# Objective extraction (your actual current return formats)
# ------------------------------------------------------------
@inline function _objective_centralized(out::Lirion.LirionRunResult)
    d = out.result
    d isa AbstractDict || error("Centralized out.result must be Dict, got $(typeof(d))")
    haskey(d, "objective") ||
        error("Centralized result has no key 'objective'. Keys=$(collect(keys(d)))")
    return d["objective"]
end

@inline function _objective_distributed_power_ac(out::Lirion.LirionRunResult)
    r = out.result
    r isa NamedTuple || error("Distributed out.result must be NamedTuple, got $(typeof(r))")
    hasproperty(r, :total_expr_cost) || error(
        "Distributed result has no field :total_expr_cost. Fields=$(propertynames(r))",
    )
    return getproperty(r, :total_expr_cost)
end

function compare_distributed_branch_ac_vs_centralized(
    file::String;
    label::String,
    central_optimizer = Ipopt.Optimizer,
    distributed_optimizer = Ipopt.Optimizer,
    rtol = 0.01,
    atol = 1e-6,

    # ---- ADMM params ----
    # New API:
    # - coupling = :power  (replaces old admm_branch for AC)
    # - sectioning = :tie_line (current only)
    # Preferred:
    # - use rho (generic) to fill required rhos for AC+:power => rho_p, rho_q
    # - OR pass both rho_p and rho_q explicitly
    admm_rho = 1.0e5,
    admm_rho_p = nothing,
    admm_rho_q = nothing,
    admm_max_iter = 5000,
    admm_eps_abs = 1e-4,
    admm_eps_rel = 1e-4,
    admm_mu_incr = 2.0,
    admm_mu_decr = 2.0,
    admm_tau_incr = 1.05,
    admm_tau_decr = 0.98,
    admm_adaptive::Union{Nothing, Bool} = nothing,

    # legacy placeholder (was CLI flags). Kept for callsites; ignored here.
    additional_flags::Vector{String} = String[],
)
    case_name = basename(file)

    # --- Centralized AC ---
    cent_obj = try
        out_cent = Lirion.solve(
            file;
            model = Lirion.ac(),
            algorithm = Lirion.centralized(),
            optimizer = central_optimizer,
            outputs = false,
            plots = false,
        )
        _objective_centralized(out_cent)
    catch e
        @warn "Centralized AC failed on $case_name: $e"
        track_failure_branch_ac(label, "Centralized", case_name)
        return
    end

    # --- Distributed ADMM-(power)/AC ---
    dist_obj = try
        alg = if admm_rho_p === nothing && admm_rho_q === nothing
            # generic rho (fills REQUIRED rho_p/rho_q for AC+coupling=:power)
            Lirion.admm(;
                coupling = :power,
                sectioning = :tie_line,
                max_iter = Int(admm_max_iter),
                rho = Float64(admm_rho),
                eps_abs = Float64(admm_eps_abs),
                eps_rel = Float64(admm_eps_rel),
                mu_incr = Float64(admm_mu_incr),
                mu_decr = Float64(admm_mu_decr),
                tau_incr = Float64(admm_tau_incr),
                tau_decr = Float64(admm_tau_decr),
                adaptive_rho = isnothing(admm_adaptive) ? true : Bool(admm_adaptive),
            )
        else
            # explicit rho_p / rho_q
            (admm_rho_p === nothing || admm_rho_q === nothing) &&
                error("If you set rho_p or rho_q, you must provide BOTH (non-nothing).")

            Lirion.admm(;
                coupling = :power,
                sectioning = :tie_line,
                max_iter = Int(admm_max_iter),
                rho_p = Float64(admm_rho_p),
                rho_q = Float64(admm_rho_q),
                eps_abs = Float64(admm_eps_abs),
                eps_rel = Float64(admm_eps_rel),
                mu_incr = Float64(admm_mu_incr),
                mu_decr = Float64(admm_mu_decr),
                tau_incr = Float64(admm_tau_incr),
                tau_decr = Float64(admm_tau_decr),
                adaptive_rho = isnothing(admm_adaptive) ? true : Bool(admm_adaptive),
            )
        end

        out_dist = Lirion.solve(
            file;
            model = Lirion.ac(),
            algorithm = alg,
            optimizer = distributed_optimizer,
            outputs = false,
            plots = false,
        )
        _objective_distributed_power_ac(out_dist)
    catch e
        @warn "Distributed ADMM-(power)-AC failed on $case_name: $e"
        track_failure_branch_ac(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized AC = $(round(cent_obj, digits=6)) | distributed-(power)-AC (total_expr_cost) = $(round(dist_obj, digits=6))"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_branch_ac_objective_tests(;
    case_list = distributed_branch_ac_zone_cases,
    label = "ADMM-(power)-AC vs AC",
    central_optimizer = Ipopt.Optimizer,
    distributed_optimizer = Ipopt.Optimizer,
    rtol = 0.01,
    atol = 1e-6,
    admm_eps_abs = 1e-4,
    admm_eps_rel = 1e-4,
    admm_mu_incr = 3.0,
    admm_mu_decr = 3.0,
    admm_tau_incr = 1.05,
    admm_tau_decr = 0.98,
)
    if isempty(case_list)
        @info "Empty distributed (power-AC) case list. No tests will be run."
        return true
    end

    test_dir = joinpath(@__DIR__, "distributed")
    selected_files =
        [joinpath(test_dir, f) for f in case_list if isfile(joinpath(test_dir, f))]

    n = length(selected_files)
    for (i, file) in enumerate(selected_files)
        case_name = basename(file)

        admm_max_iter_case = 5000
        admm_rho_case = 10000.0   # generic rho (fills rho_p/rho_q)
        # If you ever want explicit:
        admm_rho_p_case = nothing
        admm_rho_q_case = nothing

        admm_eps_abs_case = admm_eps_abs
        admm_eps_rel_case = admm_eps_rel
        mu_incr_case = admm_mu_incr
        mu_decr_case = admm_mu_decr
        tau_incr_case = admm_tau_incr
        tau_decr_case = admm_tau_decr
        adaptive_case = nothing
        extra_flags = String[]  # legacy placeholder (ignored)

        @info "[$label] Rodando caso $i de $n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rho=$(admm_rho_case), " *
              "eps_abs=$(admm_eps_abs_case), eps_rel=$(admm_eps_rel_case), " *
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
                admm_rho = admm_rho_case,
                admm_rho_p = admm_rho_p_case,
                admm_rho_q = admm_rho_q_case,
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
