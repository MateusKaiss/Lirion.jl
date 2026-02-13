# ==============================================
# ADMM-AC vs AC (objective) – tests (API solve) (compat mu/tau)
# Corrigido p/ NOVA API: admm(coupling=:state, sectioning=:tie_line, ...)
# ==============================================

# Casos particionados (2/3/4 zonas)
const distributed_ac_zone_cases = [
    # "pglib_opf_case73_ieee_rts__api.m",
    # "pglib_opf_case73_ieee_rts__sad.m",
    ## "pglib_opf_case57_ieee__api.m",
    # "pglib_opf_case57_ieee__sad.m",

    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
    "pglib_opf_case73_ieee_rts.m",
    "pglib_opf_case118_ieee.m",

    # "pglib_opf_case57_ieee.m", problem
    # "pglib_opf_case200_activ.m",
    # "pglib_opf_case500_goc.m",
]

function track_failure_ac(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] =
            Dict("Centralized" => String[], "Distributed" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

# ------------------------------------------------------------
# Objective extraction (current return formats)
# ------------------------------------------------------------
@inline function _objective_centralized_ac(out::Lirion.LirionRunResult)
    d = out.result
    d isa AbstractDict || error("Centralized out.result must be Dict, got $(typeof(d))")
    haskey(d, "objective") ||
        error("Centralized result has no key 'objective'. Keys=$(collect(keys(d)))")
    return Float64(d["objective"])
end

@inline function _objective_distributed_ac(out::Lirion.LirionRunResult)
    r = out.result
    r isa NamedTuple || error("Distributed out.result must be NamedTuple, got $(typeof(r))")
    hasproperty(r, :total_expr_cost) ||
        error(
            "Distributed result has no field :total_expr_cost. Fields=$(propertynames(r))",
        )
    return Float64(getproperty(r, :total_expr_cost))
end

# ------------------------------------------------------------
# Helper: call solve with/without optimizer kwarg
# ------------------------------------------------------------
@inline function _solve_maybe_optimizer(
    file::String;
    model,
    algorithm,
    optimizer,
    outputs = false,
    plots = false,
)
    if optimizer === nothing
        return Lirion.solve(
            file;
            model = model,
            algorithm = algorithm,
            outputs = outputs,
            plots = plots,
        )
    else
        return Lirion.solve(
            file;
            model = model,
            algorithm = algorithm,
            optimizer = optimizer,
            outputs = outputs,
            plots = plots,
        )
    end
end

function compare_distributed_ac_vs_centralized(
    file::String;
    label::String,
    central_optimizer = Ipopt.Optimizer,
    distributed_optimizer = Ipopt.Optimizer,   # pode ser nothing também
    rtol = 0.02,
    atol = 1e-6,
    admm_rho = 100.0,
    admm_max_iter = 5000,
    admm_eps_abs = 1e-4,
    admm_eps_rel = 1e-4,
    admm_mu_incr = 2.0,
    admm_mu_decr = 2.0,
    admm_tau_incr = 1.05,
    admm_tau_decr = 0.98,
    admm_rho_theta::Union{Nothing, Real} = nothing,
    admm_rho_v::Union{Nothing, Real} = nothing,
    admm_adaptive::Union{Nothing, String} = nothing, # "on" | "off" | nothing
)
    case_name = basename(file)

    # --- Centralizado AC ---
    cent_obj = try
        out_cent = Lirion.solve(
            file;
            model = Lirion.ac(),
            algorithm = Lirion.centralized(),
            optimizer = central_optimizer,
            outputs = false,
            plots = false,
        )
        _objective_centralized_ac(out_cent)
    catch e
        @warn "Centralized AC failed on $case_name: $e"
        track_failure_ac(label, "Centralized", case_name)
        return
    end

    # --- Distribuído ADMM-AC (STATE coupling) ---
    dist_obj = try
        adaptive_bool =
            (admm_adaptive === nothing) ? true :
            (lowercase(admm_adaptive) == "on") ? true :
            false

        # Se passar rho_theta/rho_v explicitamente, use-os; senão use rho (base)
        alg = if (admm_rho_theta !== nothing) || (admm_rho_v !== nothing)
            (admm_rho_theta === nothing || admm_rho_v === nothing) &&
                error("If you set rho_theta or rho_v, you must provide BOTH (non-nothing).")

            Lirion.admm(;
                coupling = :state,
                sectioning = :tie_line,
                max_iter = Int(admm_max_iter),
                rho_theta = Float64(admm_rho_theta),
                rho_v = Float64(admm_rho_v),
                eps_abs = Float64(admm_eps_abs),
                eps_rel = Float64(admm_eps_rel),
                adaptive_rho = adaptive_bool,
                mu_incr = Float64(admm_mu_incr),
                mu_decr = Float64(admm_mu_decr),
                tau_incr = Float64(admm_tau_incr),
                tau_decr = Float64(admm_tau_decr),
            )
        else
            Lirion.admm(;
                coupling = :state,
                sectioning = :tie_line,
                max_iter = Int(admm_max_iter),
                rho = Float64(admm_rho),
                eps_abs = Float64(admm_eps_abs),
                eps_rel = Float64(admm_eps_rel),
                adaptive_rho = adaptive_bool,
                mu_incr = Float64(admm_mu_incr),
                mu_decr = Float64(admm_mu_decr),
                tau_incr = Float64(admm_tau_incr),
                tau_decr = Float64(admm_tau_decr),
            )
        end

        out_dist = _solve_maybe_optimizer(
            file;
            model = Lirion.ac(),
            algorithm = alg,
            optimizer = distributed_optimizer,
            outputs = false,
            plots = false,
        )

        _objective_distributed_ac(out_dist)
    catch e
        @warn "Distributed ADMM-AC failed on $case_name: $e"
        track_failure_ac(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized AC = $(round(cent_obj, digits=6)) | distributed ADMM-AC (total_expr_cost) = $(round(dist_obj, digits=6))"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_ac_objective_tests(;
    case_list = distributed_ac_zone_cases,
    label = "ADMM-AC vs AC",
    central_optimizer = Ipopt.Optimizer,
    distributed_optimizer = Ipopt.Optimizer,
    rtol = 0.01,
    atol = 1e-6,
    admm_eps_abs = 1e-4,
    admm_eps_rel = 1e-4,
    admm_mu_incr = 2.0,
    admm_mu_decr = 2.0,
    admm_tau_incr = 1.05,
    admm_tau_decr = 0.98,
)
    if isempty(case_list)
        @info "Empty distributed (AC) case list. No tests will be run."
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

        # ----------------------------
        # Defaults (sempre definidos)
        # ----------------------------
        admm_max_iter_case = 1000
        admm_rho_case = 1.0

        admm_eps_abs_case = admm_eps_abs
        admm_eps_rel_case = admm_eps_rel

        admm_rho_theta_case = nothing
        admm_rho_v_case = nothing
        admm_adaptive_case = nothing

        admm_mu_incr_case = admm_mu_incr
        admm_mu_decr_case = admm_mu_decr
        admm_tau_incr_case = admm_tau_incr
        admm_tau_decr_case = admm_tau_decr

        central_optimizer_case = central_optimizer
        distributed_optimizer_case = distributed_optimizer

        # ----------------------------
        # Overrides por caso específico
        # ----------------------------
        if case_name == "pglib_opf_case200_activ.m" ||
           case_name == "pglib_opf_case500_goc.m"
            admm_rho_case = 1.0e7
            admm_max_iter_case = 100000
            admm_eps_abs_case = 1e-4
            admm_eps_rel_case = 1e-4
            admm_adaptive_case = "on"

            admm_mu_incr_case = 10.0
            admm_mu_decr_case = 10.0
            admm_tau_incr_case = 1.2
            admm_tau_decr_case = 0.8

            central_optimizer_case = Ipopt.Optimizer
            distributed_optimizer_case = Ipopt.Optimizer

        elseif case_name == "pglib_opf_case73_ieee_rts.m" ||
               case_name == "pglib_opf_case57_ieee.m" ||
               case_name == "pglib_opf_case118_ieee.m"
            admm_rho_case = 1.0e7
            admm_max_iter_case = 100000
            admm_eps_abs_case = 1e-4
            admm_eps_rel_case = 1e-4
            admm_adaptive_case = "on"

            admm_mu_incr_case = 10.0
            admm_mu_decr_case = 10.0
            admm_tau_incr_case = 1.05
            admm_tau_decr_case = 0.98

            central_optimizer_case = Ipopt.Optimizer
            distributed_optimizer_case = Ipopt.Optimizer
        end

        @info "[$label] Rodando caso $i/$n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rho=$(admm_rho_case), " *
              "rho_theta=$(admm_rho_theta_case), rhoV=$(admm_rho_v_case), " *
              "eps_pr=$(admm_eps_abs_case), eps_dual=$(admm_eps_rel_case), " *
              "mu=(incr=$(admm_mu_incr_case),decr=$(admm_mu_decr_case)), " *
              "tau=(incr=$(admm_tau_incr_case),decr=$(admm_tau_decr_case)), " *
              "adaptive=$(admm_adaptive_case))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_ac_vs_centralized(
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
                admm_mu_incr = admm_mu_incr_case,
                admm_mu_decr = admm_mu_decr_case,
                admm_tau_incr = admm_tau_incr_case,
                admm_tau_decr = admm_tau_decr_case,
                admm_rho_theta = admm_rho_theta_case,
                admm_rho_v = admm_rho_v_case,
                admm_adaptive = admm_adaptive_case,
            )
        end
    end

    return true
end
