# =============================
# Distributed ADMM-DC (bus) vs DC — objective tests (API solve)
# =============================

using Test
using Lirion
using OSQP
using Ipopt

# =============================
# Casos particionados (2/3 zonas)
# =============================
const distributed_zone_cases = [
    # "pglib_opf_case73_ieee_rts__api.m",
    # "pglib_opf_case73_ieee_rts__sad.m",
    ## "pglib_opf_case57_ieee__api.m",
    # "pglib_opf_case57_ieee__sad.m",

    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",

    # "pglib_opf_case57_ieee.m",
    "pglib_opf_case73_ieee_rts.m",
    # "pglib_opf_case118_ieee.m",

    # "pglib_opf_case500_goc.m",
    # "pglib_opf_case200_activ.m",
]

function track_failure_dist(label::String, origin::String, case_name::String)
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

@inline function _objective_distributed_dc_bus(out::Lirion.LirionRunResult)
    r = out.result
    r isa NamedTuple || error("Distributed out.result must be NamedTuple, got $(typeof(r))")
    hasproperty(r, :total_real_cost) ||
        error("Distributed result has no field :total_real_cost. Fields=$(propertynames(r))")
    return Float64(getproperty(r, :total_real_cost))
    # se quiser comparar o outro: Float64(getproperty(r, :total_expr_cost))
end

# ------------------------------------------------------------
# Helper: call solve with/without optimizer kwarg
# ------------------------------------------------------------
@inline function _solve_maybe_optimizer(
    file::String;
    model,
    algorithm,
    optimizer,
    outputs::Bool = false,
    plots::Bool = false,
)
    if optimizer === nothing
        return Lirion.solve(file; model=model, algorithm=algorithm, outputs=outputs, plots=plots)
    else
        return Lirion.solve(file; model=model, algorithm=algorithm, optimizer=optimizer, outputs=outputs, plots=plots)
    end
end

# ------------------------------------------------------------
# Small helper: parse "on/off/nothing" into Bool (default=true)
# ------------------------------------------------------------
@inline function _parse_adaptive(admm_adaptive::Union{Nothing, String})
    admm_adaptive === nothing && return true
    s = lowercase(admm_adaptive)
    if s == "on"
        return true
    elseif s == "off"
        return false
    else
        error("admm_adaptive must be \"on\", \"off\" or nothing; got: $admm_adaptive")
    end
end

function compare_distributed_vs_centralized(
    file::String;
    label::String,
    central_optimizer = get_available_optimizer(OSQP.Optimizer, Ipopt.Optimizer),
    distributed_optimizer::Union{Nothing, Type, Function} = nothing,  # <- FIX AQUI
    rtol = 0.05,
    atol = 1e-6,
    admm_rho = 1_000.0,
    admm_max_iter = 350,
    admm_eps_abs = 1e-6,
    admm_eps_rel = 1e-6,
    admm_alpha::Float64 = 1.0,
    admm_mu_incr::Float64 = 3.0,
    admm_mu_decr::Float64 = 3.0,
    admm_tau_incr::Float64 = 1.05,
    admm_tau_decr::Float64 = 0.98,
    admm_accel::Bool = false,                         # Nesterov
    admm_adaptive::Union{Nothing, String} = nothing,   # "on" | "off" | nothing
)
    case_name = basename(file)

    # --- Centralizado DC ---
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
        track_failure_dist(label, "Centralized", case_name)
        return
    end

    # --- Distribuído ADMM-DC (bus) ---
    dist_obj = try
        alg = Lirion.admm_bus(
            max_iter = Int(admm_max_iter),
            rho = Float64(admm_rho),
            eps_abs = Float64(admm_eps_abs),
            eps_rel = Float64(admm_eps_rel),
            alpha_overrelax = (admm_alpha == 1.0 ? nothing : Float64(admm_alpha)),
            adaptive_rho = _parse_adaptive(admm_adaptive),
            mu_incr = Float64(admm_mu_incr),
            mu_decr = Float64(admm_mu_decr),
            tau_incr = Float64(admm_tau_incr),
            tau_decr = Float64(admm_tau_decr),
            accel_nesterov = Bool(admm_accel),
        )

        out_dist = _solve_maybe_optimizer(
            file;
            model = Lirion.dc(),
            algorithm = alg,
            optimizer = distributed_optimizer,  # now accepts Ipopt.Optimizer
            outputs = false,
            plots = false,
        )

        _objective_distributed_dc_bus(out_dist)
    catch e
        @warn "Distributed ADMM-DC failed on $case_name: $e"
        track_failure_dist(label, "Distributed", case_name)
        return
    end

    @info "File $case_name | centralized = $(round(cent_obj, digits=6)) | distributed (total_real_cost) = $(round(dist_obj, digits=6))"
    @test isapprox(dist_obj, cent_obj; rtol = rtol, atol = atol)
end

# -----------------------------------------------------------
# run_distributed_objective_tests: overrides + params
# -----------------------------------------------------------
function run_distributed_objective_tests(;
    case_list = distributed_zone_cases,
    label = "ADMM-DC (bus) vs DC",
    central_optimizer = get_available_optimizer(OSQP.Optimizer, Ipopt.Optimizer),
    distributed_optimizer::Union{Nothing, Type, Function} = nothing,  # <- FIX AQUI
    rtol = 0.05,
    atol = 1e-6,
    admm_eps_abs = 1e-4,
    admm_eps_rel = 1e-4,
    admm_alpha::Float64 = 1.0,
    admm_accel::Bool = false,
    admm_adaptive::Union{Nothing, String} = nothing,
)
    if isempty(case_list)
        @info "Empty distributed case list. No tests will be run."
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
        # Defaults por zonas
        # ----------------------------
        admm_max_iter_case = 5000
        admm_rho_case = 1000.0

        admm_eps_abs_case = admm_eps_abs
        admm_eps_rel_case = admm_eps_rel
        admm_alpha_case = admm_alpha

        # parâmetros adaptativos: defaults gerais
        admm_mu_incr_case = 3.0
        admm_mu_decr_case = 3.0
        admm_tau_incr_case = 1.05
        admm_tau_decr_case = 0.98
        admm_accel_case = admm_accel
        admm_adaptive_case = admm_adaptive  # nothing | "on" | "off"

        # solvers por caso (não pisa nos globais)
        central_optimizer_case = central_optimizer
        distributed_optimizer_case = distributed_optimizer

        # ----------------------------
        # Overrides por caso específico
        # ----------------------------
        if case_name == "pglib_opf_case73_ieee_rts.m" ||
           case_name == "pglib_opf_case57_ieee.m" ||
           case_name == "pglib_opf_case118_ieee.m" ||
           case_name == "pglib_opf_case200_activ.m"

            admm_rho_case = 100000.0
            admm_max_iter_case = 100000
            admm_eps_abs_case = 1e-4
            admm_eps_rel_case = 1e-4
            admm_adaptive_case = "off"

            central_optimizer_case = Ipopt.Optimizer
            distributed_optimizer_case = Ipopt.Optimizer
        end

        @info "[$label] Rodando caso $i/$n: $case_name " *
              "(max_iter=$(admm_max_iter_case), rho=$(admm_rho_case), " *
              "eps_pr=$(admm_eps_abs_case), eps_dual=$(admm_eps_rel_case), " *
              "alpha=$(admm_alpha_case), mu_incr=$(admm_mu_incr_case), mu_decr=$(admm_mu_decr_case), " *
              "tau_incr=$(admm_tau_incr_case), tau_decr=$(admm_tau_decr_case), " *
              "accel=$(admm_accel_case), adaptive=$(admm_adaptive_case), " *
              "central=$(string(central_optimizer_case)), dist=$(distributed_optimizer_case === nothing ? "default" : string(distributed_optimizer_case)))"

        @testset "$(label) - $(case_name)" begin
            compare_distributed_vs_centralized(
                file;
                label = label,
                central_optimizer = central_optimizer_case,
                distributed_optimizer = distributed_optimizer_case,  # <- FIX: passa o case, não o global
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
                admm_accel = admm_accel_case,
                admm_adaptive = admm_adaptive_case,
            )
        end
    end

    return true
end
