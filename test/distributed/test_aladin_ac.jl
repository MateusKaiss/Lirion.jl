# ============================================================
# =======================  ALADIN AC TESTS  =====================
# ============================================================

const ALADIN_AC_STATE_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

const ALADIN_AC_POWER_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

const ALADIN_AC_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

const ALADIN_AC_SPLITBUS_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

const ALADIN_AC_DUMMYBUS_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

# Power-only DummyBus excludes __api (tight bounds can trip infeasibility
# before θ/V discipline kicks in via v² / bus consensus).
const ALADIN_AC_DUMMYBUS_POWER_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

# ============================================================
# OBJECTIVE EXTRACTORS
# ============================================================
aladin_central_obj_ac(out) = Float64(out.objective)
aladin_dist_obj_ac(out) = Float64(out.objective)

# ============================================================
# ALADIN AC - STATE
# ============================================================
function run_aladin_ac_state_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = @__DIR__

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
    test_dir = @__DIR__

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
    test_dir = @__DIR__

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

# ============================================================
# ALADIN AC - SPLIT-BUS HYBRID (Method 3 — Harris et al. 2022, Fig. 1c)
# ============================================================
function run_aladin_ac_splitbus_hybrid_tests(;
    rtol = 0.05,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = @__DIR__

    make_alg(p) = Lirion.aladin(;
        coupling = :hybrid,
        sectioning = Lirion.SplitBus(),
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
        label = "ALADIN-AC SPLITBUS-HYBRID",
        test_dir = test_dir,
        case_list = ALADIN_AC_SPLITBUS_HYBRID_CASES,
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
# ALADIN AC - DUMMY-BUS HYBRID (Method 2 — Harris et al. 2022)
# ============================================================
function run_aladin_ac_dummybus_hybrid_tests(;
    rtol = 0.05,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = @__DIR__

    make_alg(p) = Lirion.aladin(;
        coupling = :hybrid,
        sectioning = Lirion.DummyBus(),
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
        max_iter = 2000,
        rho = 10000.0,
        mu = 10000.0,
        rho_growth = 1.05,
        mu_growth = 1.1,
        adaptive = true,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
    )

    overrides(case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-AC DUMMYBUS-HYBRID",
        test_dir = test_dir,
        case_list = ALADIN_AC_DUMMYBUS_HYBRID_CASES,
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
# ALADIN AC - DUMMY-BUS POWER (Method 2, power-only)
# ============================================================
function run_aladin_ac_dummybus_power_tests(;
    rtol = 0.05,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = @__DIR__

    make_alg(p) = Lirion.aladin(;
        coupling = :power,
        sectioning = Lirion.DummyBus(),
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
        max_iter = 2000,
        rho = 10000.0,
        mu = 10000.0,
        rho_growth = 1.05,
        mu_growth = 1.1,
        adaptive = true,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-AC DUMMYBUS-POWER",
        test_dir = test_dir,
        case_list = ALADIN_AC_DUMMYBUS_POWER_CASES,
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
