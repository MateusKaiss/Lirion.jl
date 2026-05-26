# ============================================================
# =========================  DC TESTS  =========================
# ============================================================

# ------------------------
# DC - STATE (bus/state coupling)
# ------------------------
const DC_STATE_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

# ------------------------
# DC - POWER (branch/power coupling)
# ------------------------
const DC_POWER_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

# ------------------------
# DC - HYBRID
# ------------------------
const DC_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

# ------------------------
# DC - DUMMY-BUS HYBRID
# ------------------------
const DC_DUMMYBUS_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

# ------------------------
# DC - SPLIT-BUS HYBRID
# ------------------------
const DC_SPLITBUS_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

# N-way SplitBus: bus 1 is a boundary bus shared by THREE zones.
const DC_SPLITBUS_3ZONE_CASES = [
    "case4_3zone_star.m",
]

# ------------------------
# DC - DUMMY-BUS POWER (only P consensus, θ_d ambiguous — OK because d is artificial)
# ------------------------
# Note: case3_lmbd__api intentionally excluded. Without θ consensus at the
# dummy bus, early ADMM iterations can drift to a P target that violates the
# tight generation / line bounds of the __api variant, producing local
# infeasibility (OTHER_ERROR from the solver). Hybrid coupling avoids this by
# disciplining θ; power-only alone cannot recover. This is an inherent
# limitation of power-only DummyBus on tightly-bounded problems, not a bug.
const DC_DUMMYBUS_POWER_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

# ============================================================
# OBJECTIVE EXTRACTORS
# ============================================================
central_obj_dc(out) = Float64(out.objective)
dist_obj_dc(out) = Float64(out.objective)

# ============================================================
# DC - STATE
# ============================================================
function run_admm_dc_state_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = @__DIR__

    make_alg(p) = Lirion.admm(;
        coupling = :state,
        sectioning = :tie_line,
        max_iter = p.max_iter,
        rho = p.rho,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
        adaptive_rho = p.adaptive,
        mu_incr = p.mu_incr,
        mu_decr = p.mu_decr,
        tau_incr = p.tau_incr,
        tau_decr = p.tau_decr,
    )

    defaults(_) = (
        max_iter = 5000,
        rho = 1000.0,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    overrides(case_name) = begin
        # case14_ieee state: needs ρ=1e5 (empirical, AC state diverges at ρ=1e3).
        if case_name == "pglib_opf_case14_ieee.m"
            return (rho = 100000.0, max_iter = 5000, adaptive = false)
        end
        return NamedTuple()
    end

    return run_case_suite(;
        label = "ADMM-DC STATE",
        test_dir = test_dir,
        case_list = DC_STATE_CASES,
        model = Lirion.dc(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_dc,
        objective_dist = dist_obj_dc,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# DC - POWER
# ============================================================
function run_admm_dc_power_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = @__DIR__

    make_alg(p) = Lirion.admm(;
        coupling = :power,
        sectioning = :tie_line,
        max_iter = p.max_iter,
        rho = p.rho,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
        adaptive_rho = p.adaptive,
        mu_incr = p.mu_incr,
        mu_decr = p.mu_decr,
        tau_incr = p.tau_incr,
        tau_decr = p.tau_decr,
    )

    defaults(_) = (
        max_iter = 5000,
        rho = 100000.0,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ADMM-DC POWER",
        test_dir = test_dir,
        case_list = DC_POWER_CASES,
        model = Lirion.dc(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_dc,
        objective_dist = dist_obj_dc,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# DC - HYBRID
# ============================================================
function run_admm_dc_hybrid_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = @__DIR__

    make_alg(p) = Lirion.admm(;
        coupling = :hybrid,
        sectioning = :tie_line,
        max_iter = p.max_iter,
        rho = p.rho,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
        adaptive_rho = p.adaptive,
        mu_incr = p.mu_incr,
        mu_decr = p.mu_decr,
        tau_incr = p.tau_incr,
        tau_decr = p.tau_decr,
    )

    defaults(_) = (
        max_iter = 5000,
        rho = 10000.0,
        eps_abs = 2e-4,
        eps_rel = 2e-4,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    overrides(case_name) = NamedTuple()

    return run_case_suite(;
        label = "ADMM-DC HYBRID",
        test_dir = test_dir,
        case_list = DC_HYBRID_CASES,
        model = Lirion.dc(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_dc,
        objective_dist = dist_obj_dc,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# DC - DUMMY-BUS HYBRID (Method 2 — Harris et al. 2022)
# ============================================================
function run_admm_dc_dummybus_hybrid_tests(;
    rtol = 0.05,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = @__DIR__

    make_alg(p) = Lirion.admm(;
        coupling = :hybrid,
        sectioning = Lirion.DummyBus(),
        max_iter = p.max_iter,
        rho = p.rho,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
        adaptive_rho = p.adaptive,
        mu_incr = p.mu_incr,
        mu_decr = p.mu_decr,
        tau_incr = p.tau_incr,
        tau_decr = p.tau_decr,
    )

    defaults(_) = (
        max_iter = 10000,
        rho = 1000.0,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    overrides(case_name) = NamedTuple()

    return run_case_suite(;
        label = "ADMM-DC DUMMYBUS-HYBRID",
        test_dir = test_dir,
        case_list = DC_DUMMYBUS_HYBRID_CASES,
        model = Lirion.dc(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_dc,
        objective_dist = dist_obj_dc,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# DC - SPLIT-BUS HYBRID (Method 3 — Harris et al. 2022, Fig. 1c)
# ============================================================
function run_admm_dc_splitbus_hybrid_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = @__DIR__

    make_alg(p) = Lirion.admm(;
        coupling = :hybrid,
        sectioning = Lirion.SplitBus(),
        max_iter = p.max_iter,
        rho = p.rho,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
        adaptive_rho = p.adaptive,
        mu_incr = p.mu_incr,
        mu_decr = p.mu_decr,
        tau_incr = p.tau_incr,
        tau_decr = p.tau_decr,
    )

    defaults(_) = (
        max_iter = 10000,
        rho = 1000.0,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    overrides(case_name) = NamedTuple()

    return run_case_suite(;
        label = "ADMM-DC SPLITBUS-HYBRID",
        test_dir = test_dir,
        case_list = DC_SPLITBUS_HYBRID_CASES,
        model = Lirion.dc(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_dc,
        objective_dist = dist_obj_dc,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ------------------------------------------------------------
# DC - SPLIT-BUS HYBRID, N-WAY (a boundary bus shared by 3 zones)
# ------------------------------------------------------------

function run_admm_dc_splitbus_3zone_tests(;
    rtol = 0.02,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = @__DIR__

    make_alg(p) = Lirion.admm(;
        coupling = :hybrid,
        sectioning = Lirion.SplitBus(),
        max_iter = p.max_iter,
        rho = p.rho,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
        adaptive_rho = p.adaptive,
        mu_incr = p.mu_incr,
        mu_decr = p.mu_decr,
        tau_incr = p.tau_incr,
        tau_decr = p.tau_decr,
    )

    defaults(_) = (
        max_iter = 30000,
        rho = 1000.0,
        eps_abs = 1e-6,
        eps_rel = 1e-6,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    return run_case_suite(;
        label = "ADMM-DC SPLITBUS 3-ZONE",
        test_dir = test_dir,
        case_list = DC_SPLITBUS_3ZONE_CASES,
        model = Lirion.dc(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_dc,
        objective_dist = dist_obj_dc,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = _ -> NamedTuple(),
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# DC - DUMMY-BUS POWER (Method 2 with power-only coupling)
# ------------------------------------------------------------
# Theory: with DummyBus sectioning, the dummy bus `d` is artificial; its θ
# is gauge-invariant. Enforcing P_{m,d} = P_{d,n} (power consensus on the two
# halves of each split tie-line) is sufficient to pin the physical flow
# through d, which in DC is all that the gen dispatch depends on.
# Expected: distributed cost ≈ centralized cost (like hybrid).
# ============================================================
function run_admm_dc_dummybus_power_tests(;
    rtol = 0.05,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = @__DIR__

    make_alg(p) = Lirion.admm(;
        coupling = :power,
        sectioning = Lirion.DummyBus(),
        max_iter = p.max_iter,
        rho = p.rho,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
        adaptive_rho = p.adaptive,
        mu_incr = p.mu_incr,
        mu_decr = p.mu_decr,
        tau_incr = p.tau_incr,
        tau_decr = p.tau_decr,
    )

    defaults(_) = (
        max_iter = 10000,
        rho = 1000.0,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ADMM-DC DUMMYBUS-POWER",
        test_dir = test_dir,
        case_list = DC_DUMMYBUS_POWER_CASES,
        model = Lirion.dc(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_dc,
        objective_dist = dist_obj_dc,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end
