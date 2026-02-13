# ==============================================
# ADMM-Hybrid-AC vs AC (objective) – tests (API solve)
# Corrigido p/ NOVA API: admm(coupling=:hybrid, sectioning=:tie_line, ...)
# ==============================================

const distributed_hybrid_ac_zone_cases = [
    # "pglib_opf_case73_ieee_rts__api.m",
    # "pglib_opf_case73_ieee_rts__sad.m",
    ## "pglib_opf_case57_ieee__api.m",
    # "pglib_opf_case57_ieee__sad.m",

    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
    "pglib_opf_case73_ieee_rts.m", "pglib_opf_case57_ieee.m",
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

# ------------------------------------------------------------
# Objective extraction (same style as your branch-AC test)
# ------------------------------------------------------------
@inline function _objective_centralized_ac(out::Lirion.LirionRunResult)
    d = out.result
    d isa AbstractDict || error("Centralized out.result must be Dict, got $(typeof(d))")
    haskey(d, "objective") ||
        error("Centralized result has no key 'objective'. Keys=$(collect(keys(d)))")
    return d["objective"]
end

@inline function _objective_distributed_hybrid_ac(out::Lirion.LirionRunResult)
    r = out.result
    r isa NamedTuple || error("Distributed out.result must be NamedTuple, got $(typeof(r))")
    hasproperty(r, :total_expr_cost) || error(
        "Distributed result has no field :total_expr_cost. Fields=$(propertynames(r))",
    )
    return getproperty(r, :total_expr_cost)
end

function compare_distributed_hybrid_ac_vs_centralized(
    file::String;
    label::String,
    central_optimizer = Ipopt.Optimizer,
    distributed_optimizer = Ipopt.Optimizer,
    rtol = 0.01,
    atol = 1e-6,

    # rhos (AC + hybrid => rho_theta, rho_v, rho_p, rho_q)
    rho_base::Union{Nothing, Float64} = 1.0e3,
    rho_theta::Union{Nothing, Float64} = nothing,
    rho_v::Union{Nothing, Float64} = nothing,
    rho_p::Union{Nothing, Float64} = 1.0e5,
    rho_q::Union{Nothing, Float64} = 1.0e5, admm_max_iter::Int = 500,
    admm_eps_abs::Float64 = 1e-7,
    admm_eps_rel::Float64 = 1e-7,
    admm_mu_incr::Float64 = 2.0,
    admm_mu_decr::Float64 = 2.0,
    admm_tau_incr::Float64 = 1.05,
    admm_tau_decr::Float64 = 0.98,
    admm_adaptive::Union{Nothing, Bool} = nothing, # nothing => default true no builder
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
        _objective_centralized_ac(out_cent)
    catch e
        @warn "Centralized AC failed on $case_name: $e"
        track_failure_hybrid_ac(label, "Centralized", case_name)
        return
    end

    # --- Distributed ADMM-Hybrid-AC (NEW API) ---
    dist_obj = try
        # Regra:
        # - Se você passar rho_* explícitos, NÃO passe rho=...
        # - Se quiser preencher tudo com um único valor, passe só rho=...
        use_explicit_rhos =
            (rho_theta !== nothing) || (rho_v !== nothing) ||
            (rho_p !== nothing) || (rho_q !== nothing)

        alg = if use_explicit_rhos
            Lirion.admm(;
                coupling = :hybrid,
                sectioning = :tie_line,
                max_iter = Int(admm_max_iter),
                rho_theta = (rho_theta === nothing ? nothing : Float64(rho_theta)),
                rho_v = (rho_v === nothing ? nothing : Float64(rho_v)),
                rho_p = (rho_p === nothing ? nothing : Float64(rho_p)),
                rho_q = (rho_q === nothing ? nothing : Float64(rho_q)),
                eps_abs = Float64(admm_eps_abs),
                eps_rel = Float64(admm_eps_rel),
                adaptive_rho = isnothing(admm_adaptive) ? true : Bool(admm_adaptive),
                mu_incr = Float64(admm_mu_incr),
                mu_decr = Float64(admm_mu_decr),
                tau_incr = Float64(admm_tau_incr),
                tau_decr = Float64(admm_tau_decr),
            )
        else
            Lirion.admm(;
                coupling = :hybrid,
                sectioning = :tie_line,
                max_iter = Int(admm_max_iter),
                rho = Float64(rho_base === nothing ? 1.0e3 : rho_base),
                eps_abs = Float64(admm_eps_abs),
                eps_rel = Float64(admm_eps_rel),
                adaptive_rho = isnothing(admm_adaptive) ? true : Bool(admm_adaptive),
                mu_incr = Float64(admm_mu_incr),
                mu_decr = Float64(admm_mu_decr),
                tau_incr = Float64(admm_tau_incr),
                tau_decr = Float64(admm_tau_decr),
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
        _objective_distributed_hybrid_ac(out_dist)
    catch e
        @warn "Distributed ADMM-Hybrid-AC failed on $case_name: $e"
        track_failure_hybrid_ac(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized AC = $(round(cent_obj, digits=6)) | distributed-hybrid-AC (total_expr_cost) = $(round(dist_obj, digits=6))"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

function run_distributed_hybrid_ac_objective_tests(;
    case_list = distributed_hybrid_ac_zone_cases,
    label = "ADMM-Hybrid-AC vs AC",
    central_optimizer = Ipopt.Optimizer,
    distributed_optimizer = Ipopt.Optimizer,
    rtol = 0.01,
    atol = 1e-6,
    admm_eps_abs = 1e-4,
    admm_eps_rel = 1e-4,
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
        admm_adaptive_case::Union{Nothing, Bool} = nothing

        # ---- overrides (mantendo tua lógica) ----
        if case_name == "pglib_opf_case118_2zones_ieee.m" ||
           case_name == "pglib_opf_case118_3zones_ieee.m" ||
           case_name == "pglib_opf_case118_5zones_ieee.m"
            admm_max_iter_case = 500
            admm_eps_abs_case = 1e-4
            admm_eps_rel_case = 1e-4
            rho_bus_base = 1.0e5
            rho_flow_base = 1.0e5
            admm_adaptive_case = true

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
            rho_bus_base = 1e4
            rho_flow_base = 1e4
            admm_adaptive_case = true
            admm_mu_incr_case = 3.0
            admm_mu_decr_case = 3.0
            admm_tau_incr_case = 1.1
            admm_tau_decr_case = 0.9
        end

        if case_name == "pglib_opf_case2000_5zones_goc.m"
            admm_max_iter_case = 500
            admm_eps_abs_case = 1e-4
            admm_eps_rel_case = 1e-4
            rho_bus_base = 1.0e4
            rho_flow_base = 1.0e4
            admm_adaptive_case = true
            admm_mu_incr_case, admm_mu_decr_case = 2.0, 2.0
            admm_tau_incr_case, admm_tau_decr_case = 1.1, 0.9
        end

        if case_name == "pglib_opf_case3970_25zones_goc.m"
            admm_max_iter_case = 500
            admm_eps_abs_case = 3e-4
            admm_eps_rel_case = 3e-4
            rho_bus_base = 1.0e4
            rho_flow_base = 1.0e4
            admm_adaptive_case = true
            admm_mu_incr_case, admm_mu_decr_case = 2.0, 2.0
            admm_tau_incr_case, admm_tau_decr_case = 1.3, 0.7
        end

        if case_name == "pglib_opf_case10000_10zones_ac_goc.m"
            admm_max_iter_case = 800
            admm_eps_abs_case = 3e-4
            admm_eps_rel_case = 3e-4
            rho_bus_base = 1.0e3
            rho_flow_base = 1.0e3
            admm_adaptive_case = true
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
                atol = atol, rho_base = nothing,
                rho_theta = rho_bus_base,
                rho_v = rho_bus_base,
                rho_p = rho_flow_base,
                rho_q = rho_flow_base, admm_max_iter = admm_max_iter_case,
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
