# ============================================================
# =======================  ALADIN AC TESTS  =====================
# ============================================================

const ALADIN_AC_STATE_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
]

const ALADIN_AC_POWER_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
]

const ALADIN_AC_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
]

# ============================================================
# OBJECTIVE EXTRACTORS
# ============================================================
aladin_central_obj_ac(out) = Float64(out.result["objective"])
aladin_dist_obj_ac(out) = Float64(out.result.total_cost)

# ============================================================
# ALADIN AC - STATE
# ============================================================
function run_aladin_ac_state_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = joinpath(@__DIR__, "distributed")

    make_alg(p) = Lirion.aladin(;
        coupling = :state,
        sectioning = :tie_line,
        max_iter = p.max_iter,
        rho = p.rho,
        mu = p.mu,
        rho_growth = p.rho_growth,
        mu_growth = p.mu_growth,
        adaptive = p.adaptive,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
    )

    defaults(_) = (
        max_iter = 1000,
        rho = 1000.0,
        mu = 1000.0,
        rho_growth = 1.05,
        mu_growth = 1.1,
        adaptive = true,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-AC STATE",
        test_dir = test_dir,
        case_list = ALADIN_AC_STATE_CASES,
        model = Lirion.ac(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_ac,
        objective_dist = aladin_dist_obj_ac,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# ALADIN AC - POWER
# ============================================================
function run_aladin_ac_power_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = joinpath(@__DIR__, "distributed")

    make_alg(p) = Lirion.aladin(;
        coupling = :power,
        sectioning = :tie_line,
        max_iter = p.max_iter,
        rho = p.rho,
        mu = p.mu,
        rho_growth = p.rho_growth,
        mu_growth = p.mu_growth,
        adaptive = p.adaptive,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
    )

    defaults(_) = (
        max_iter = 1000,
        rho = 1000.0,
        mu = 1000.0,
        rho_growth = 1.05,
        mu_growth = 1.1,
        adaptive = true,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-AC POWER",
        test_dir = test_dir,
        case_list = ALADIN_AC_POWER_CASES,
        model = Lirion.ac(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_ac,
        objective_dist = aladin_dist_obj_ac,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# ALADIN AC - HYBRID
# ============================================================
function run_aladin_ac_hybrid_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = joinpath(@__DIR__, "distributed")

    make_alg(p) = Lirion.aladin(;
        coupling = :hybrid,
        sectioning = :tie_line,
        max_iter = p.max_iter,
        rho = p.rho,
        mu = p.mu,
        rho_growth = p.rho_growth,
        mu_growth = p.mu_growth,
        adaptive = p.adaptive,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
    )

    defaults(_) = (
        max_iter = 1000,
        rho = 1000.0,
        mu = 1000.0,
        rho_growth = 1.05,
        mu_growth = 1.1,
        adaptive = true,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-AC HYBRID",
        test_dir = test_dir,
        case_list = ALADIN_AC_HYBRID_CASES,
        model = Lirion.ac(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_ac,
        objective_dist = aladin_dist_obj_ac,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end
