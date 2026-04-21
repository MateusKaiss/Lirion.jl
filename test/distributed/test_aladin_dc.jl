# ============================================================
# =======================  ALADIN DC TESTS  =====================
# ============================================================

# ------------------------
# ALADIN DC - STATE
# ------------------------
const ALADIN_DC_STATE_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
    "pglib_opf_case30_ieee.m",
]

const ALADIN_DC_POWER_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
    "pglib_opf_case30_ieee.m",
]

const ALADIN_DC_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
    "pglib_opf_case30_ieee.m",
]

const ALADIN_DC_SPLITBUS_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
    "pglib_opf_case30_ieee.m",
]

const ALADIN_DC_DUMMYBUS_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
    "pglib_opf_case30_ieee.m",
]

# Power-only DummyBus: case3_lmbd__api excluded (same issue as ADMM — without θ
# consensus at the dummy bus, tight bounds can drive early iterations to
# infeasibility).
const ALADIN_DC_DUMMYBUS_POWER_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
    "pglib_opf_case30_ieee.m",
]

# ============================================================
# OBJECTIVE EXTRACTORS
# ============================================================
aladin_central_obj_dc(out) = Float64(out.objective)
aladin_dist_obj_dc(out) = Float64(out.objective)

# ============================================================
# ALADIN DC - STATE
# ============================================================
function run_aladin_dc_state_tests(;
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
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
    )

    defaults(_) = (
        max_iter = 5000,
        rho = 1000.0,
        mu = 1000.0,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-DC STATE",
        test_dir = test_dir,
        case_list = ALADIN_DC_STATE_CASES,
        model = Lirion.dc(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_dc,
        objective_dist = aladin_dist_obj_dc,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# ALADIN DC - POWER
# ============================================================
function run_aladin_dc_power_tests(;
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
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
    )

    defaults(_) = (
        max_iter = 5000,
        rho = 100000.0,
        mu = 1000.0,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-DC POWER",
        test_dir = test_dir,
        case_list = ALADIN_DC_POWER_CASES,
        model = Lirion.dc(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_dc,
        objective_dist = aladin_dist_obj_dc,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# ALADIN DC - HYBRID
# ============================================================
function run_aladin_dc_hybrid_tests(;
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
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
    )

    defaults(_) = (
        max_iter = 5000,
        rho = 10000.0,
        mu = 1000.0,
        eps_abs = 2e-4,
        eps_rel = 2e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-DC HYBRID",
        test_dir = test_dir,
        case_list = ALADIN_DC_HYBRID_CASES,
        model = Lirion.dc(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_dc,
        objective_dist = aladin_dist_obj_dc,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# ALADIN DC - SPLIT-BUS HYBRID (Method 3 — Harris et al. 2022, Fig. 1c)
# ============================================================
function run_aladin_dc_splitbus_hybrid_tests(;
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
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
    )

    defaults(_) = (
        max_iter = 5000,
        rho = 10000.0,
        mu = 1000.0,
        eps_abs = 2e-4,
        eps_rel = 2e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-DC SPLITBUS-HYBRID",
        test_dir = test_dir,
        case_list = ALADIN_DC_SPLITBUS_HYBRID_CASES,
        model = Lirion.dc(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_dc,
        objective_dist = aladin_dist_obj_dc,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# ALADIN DC - DUMMY-BUS HYBRID (Method 2 — Harris et al. 2022)
# ============================================================
function run_aladin_dc_dummybus_hybrid_tests(;
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
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
    )

    defaults(_) = (
        max_iter = 5000,
        rho = 10000.0,
        mu = 1000.0,
        eps_abs = 2e-4,
        eps_rel = 2e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-DC DUMMYBUS-HYBRID",
        test_dir = test_dir,
        case_list = ALADIN_DC_DUMMYBUS_HYBRID_CASES,
        model = Lirion.dc(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_dc,
        objective_dist = aladin_dist_obj_dc,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# ALADIN DC - DUMMY-BUS POWER (Method 2 with power-only coupling)
# ============================================================
function run_aladin_dc_dummybus_power_tests(;
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
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
    )

    defaults(_) = (
        max_iter = 5000,
        rho = 100000.0,
        mu = 1000.0,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-DC DUMMYBUS-POWER",
        test_dir = test_dir,
        case_list = ALADIN_DC_DUMMYBUS_POWER_CASES,
        model = Lirion.dc(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_dc,
        objective_dist = aladin_dist_obj_dc,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end
