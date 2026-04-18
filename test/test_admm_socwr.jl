# ============================================================
# =========================  SOCWR TESTS  =======================
# ============================================================

# ------------------------
# SOCWR - STATE
# ------------------------
const SOCWR_STATE_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
    "pglib_opf_case73_ieee_rts.m",
    # "pglib_opf_case200_activ.m",

    # "pglib_opf_case118_ieee.m",
    # "pglib_opf_case57_ieee.m",
    # "pglib_opf_case500_goc.m",
]

# ------------------------
# SOCWR - POWER
# ------------------------
const SOCWR_POWER_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
    "pglib_opf_case73_ieee_rts.m",
    # "pglib_opf_case57_ieee.m",
    # "pglib_opf_case118_ieee.m",
    # "pglib_opf_case200_activ.m",
    # "pglib_opf_case500_goc.m",
]

# ------------------------
# SOCWR - HYBRID
# ------------------------
const SOCWR_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
    "pglib_opf_case73_ieee_rts.m",
    # "pglib_opf_case57_ieee.m",
    # "pglib_opf_case118_ieee.m",
    # "pglib_opf_case200_activ.m",
    # "pglib_opf_case500_goc.m",
]

# ------------------------
# SOCWR - SPLIT-BUS HYBRID
# ------------------------
const SOCWR_SPLITBUS_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__sad.m",
]

# ------------------------
# SOCWR - DUMMY-BUS HYBRID
# ------------------------
const SOCWR_DUMMYBUS_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__sad.m",
]

# ============================================================
# OBJECTIVE EXTRACTORS
# ============================================================
central_obj_socwr(out) = Float64(out.result["objective"])
dist_obj_socwr(out) = Float64(out.result.total_cost)

# ============================================================
# SOCWR - STATE
# ============================================================
function run_admm_socwr_state_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer,
    distributed_optimizer,
)
    test_dir = joinpath(@__DIR__, "distributed")

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
        max_iter = 10_000,
        rho = 1e6,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
        adaptive = true,
        mu_incr = 10.0,
        mu_decr = 10.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    return run_case_suite(;
        label = "ADMM-SOCWR STATE",
        test_dir = test_dir,
        case_list = SOCWR_STATE_CASES,
        model = Lirion.socwr(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_socwr,
        objective_dist = dist_obj_socwr,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = _ -> NamedTuple(),
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# SOCWR - POWER
# ============================================================
function run_admm_socwr_power_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer,
    distributed_optimizer,
)
    test_dir = joinpath(@__DIR__, "distributed")

    make_alg(p) = Lirion.admm(;
        coupling = :power,
        sectioning = :tie_line,
        max_iter = p.max_iter,
        rho_p = p.rho_p,
        rho_q = p.rho_q,
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
        rho_p = 1e3,
        rho_q = 1e3,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    return run_case_suite(;
        label = "ADMM-SOCWR POWER",
        test_dir = test_dir,
        case_list = SOCWR_POWER_CASES,
        model = Lirion.socwr(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_socwr,
        objective_dist = dist_obj_socwr,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = _ -> NamedTuple(),
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# SOCWR - HYBRID
# ============================================================
function run_admm_socwr_hybrid_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer,
    distributed_optimizer,
)
    test_dir = joinpath(@__DIR__, "distributed")

    make_alg(p) = Lirion.admm(;
        coupling = :hybrid,
        sectioning = :tie_line,
        max_iter = p.max_iter,
        rho_v2 = p.rho_v2,
        rho_p = p.rho_p,
        rho_q = p.rho_q,
        rho_r = p.rho_r,
        rho_t = p.rho_t,
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
        rho_v2 = 1e4,
        rho_p = 3e2,
        rho_q = 3e2,
        rho_r = 1e4,
        rho_t = 1e4,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    return run_case_suite(;
        label = "ADMM-SOCWR HYBRID",
        test_dir = test_dir,
        case_list = SOCWR_HYBRID_CASES,
        model = Lirion.socwr(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_socwr,
        objective_dist = dist_obj_socwr,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = _ -> NamedTuple(),
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# SOCWR - SPLIT-BUS HYBRID (Method 3 — Harris et al. 2022, Fig. 1c)
# ============================================================
function run_admm_socwr_splitbus_hybrid_tests(;
    rtol = 0.05,
    atol = 1e-6,
    central_optimizer,
    distributed_optimizer,
)
    test_dir = joinpath(@__DIR__, "distributed")

    make_alg(p) = Lirion.admm(;
        coupling = :hybrid,
        sectioning = Lirion.SplitBus(),
        max_iter = p.max_iter,
        rho_v2 = p.rho_v2,
        rho_p = p.rho_p,
        rho_q = p.rho_q,
        rho_r = p.rho_r,
        rho_t = p.rho_t,
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
        rho_v2 = 1e4,
        rho_p = 3e2,
        rho_q = 3e2,
        rho_r = 1e4,
        rho_t = 1e4,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    return run_case_suite(;
        label = "ADMM-SOCWR SPLITBUS-HYBRID",
        test_dir = test_dir,
        case_list = SOCWR_SPLITBUS_HYBRID_CASES,
        model = Lirion.socwr(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_socwr,
        objective_dist = dist_obj_socwr,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = _ -> NamedTuple(),
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# SOCWR - DUMMY-BUS HYBRID (Method 2 — Harris et al. 2022, Fig. 1b)
# ============================================================
function run_admm_socwr_dummybus_hybrid_tests(;
    rtol = 0.05,
    atol = 1e-6,
    central_optimizer,
    distributed_optimizer,
)
    test_dir = joinpath(@__DIR__, "distributed")

    make_alg(p) = Lirion.admm(;
        coupling = :hybrid,
        sectioning = Lirion.DummyBus(),
        max_iter = p.max_iter,
        rho_v2 = p.rho_v2,
        rho_p = p.rho_p,
        rho_q = p.rho_q,
        rho_r = p.rho_r,
        rho_t = p.rho_t,
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
        rho_v2 = 1e4,
        rho_p = 3e2,
        rho_q = 3e2,
        rho_r = 1e4,
        rho_t = 1e4,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    return run_case_suite(;
        label = "ADMM-SOCWR DUMMYBUS-HYBRID",
        test_dir = test_dir,
        case_list = SOCWR_DUMMYBUS_HYBRID_CASES,
        model = Lirion.socwr(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_socwr,
        objective_dist = dist_obj_socwr,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = _ -> NamedTuple(),
        rtol = rtol,
        atol = atol,
    )
end
