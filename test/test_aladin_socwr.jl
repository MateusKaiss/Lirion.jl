# ============================================================
# =====================  ALADIN SOCWR TESTS  ====================
# ============================================================

const ALADIN_SOCWR_STATE_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
]

const ALADIN_SOCWR_POWER_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
]

const ALADIN_SOCWR_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
]

const ALADIN_SOCWR_SPLITBUS_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
]

# ============================================================
# OBJECTIVE EXTRACTORS
# ============================================================
aladin_central_obj_socwr(out) = Float64(out.result["objective"])
aladin_dist_obj_socwr(out) = Float64(out.result.total_cost)

# ============================================================
# ALADIN SOCWR - STATE
# ============================================================
function run_aladin_socwr_state_tests(;
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
        alpha1 = p.alpha1,
        alpha2 = p.alpha2,
        alpha3 = p.alpha3,
        rho_growth = p.rho_growth,
        mu_growth = p.mu_growth,
        adaptive = p.adaptive,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
    )

    defaults(_) = (
        max_iter = 500,
        rho = 1e4,
        mu = 1e4,
        # Full Newton step as initial value — the canonical Armijo
        # line search in phase2_socwr.jl backtracks α₂ per iteration
        # against the HFD augmented-Lagrangian merit function and
        # accepts α₂ = 1 whenever the active cone faces are projected
        # (the common case after the cone-face patch in
        # sensitivities.jl).
        alpha1 = 1.0,
        alpha2 = 1.0,
        alpha3 = 1.0,
        rho_growth = 1.0,
        mu_growth = 1.0,
        adaptive = false,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-SOCWR STATE",
        test_dir = test_dir,
        case_list = ALADIN_SOCWR_STATE_CASES,
        model = Lirion.socwr(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_socwr,
        objective_dist = aladin_dist_obj_socwr,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# ALADIN SOCWR - POWER
# ============================================================
function run_aladin_socwr_power_tests(;
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
        alpha1 = p.alpha1,
        alpha2 = p.alpha2,
        alpha3 = p.alpha3,
        rho_growth = p.rho_growth,
        mu_growth = p.mu_growth,
        adaptive = p.adaptive,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
    )

    defaults(_) = (
        max_iter = 500,
        rho = 1e4,
        mu = 1e4,
        alpha1 = 1.0,
        alpha2 = 1.0,
        alpha3 = 1.0,
        rho_growth = 1.0,
        mu_growth = 1.0,
        adaptive = false,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-SOCWR POWER",
        test_dir = test_dir,
        case_list = ALADIN_SOCWR_POWER_CASES,
        model = Lirion.socwr(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_socwr,
        objective_dist = aladin_dist_obj_socwr,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# ALADIN SOCWR - HYBRID
# ============================================================
function run_aladin_socwr_hybrid_tests(;
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
        alpha1 = p.alpha1,
        alpha2 = p.alpha2,
        alpha3 = p.alpha3,
        rho_growth = p.rho_growth,
        mu_growth = p.mu_growth,
        adaptive = p.adaptive,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
    )

    defaults(_) = (
        max_iter = 500,
        rho = 1e4,
        mu = 1e4,
        alpha1 = 1.0,
        alpha2 = 1.0,
        alpha3 = 1.0,
        rho_growth = 1.0,
        mu_growth = 1.0,
        adaptive = false,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-SOCWR HYBRID",
        test_dir = test_dir,
        case_list = ALADIN_SOCWR_HYBRID_CASES,
        model = Lirion.socwr(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_socwr,
        objective_dist = aladin_dist_obj_socwr,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# ALADIN SOCWR - SPLIT-BUS HYBRID (Method 3 — Harris et al. 2022, Fig. 1c)
# ============================================================
function run_aladin_socwr_splitbus_hybrid_tests(;
    rtol = 0.05,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = joinpath(@__DIR__, "distributed")

    make_alg(p) = Lirion.aladin(;
        coupling = :hybrid,
        sectioning = Lirion.SplitBus(),
        max_iter = p.max_iter,
        rho = p.rho,
        mu = p.mu,
        alpha1 = p.alpha1,
        alpha2 = p.alpha2,
        alpha3 = p.alpha3,
        rho_growth = p.rho_growth,
        mu_growth = p.mu_growth,
        adaptive = p.adaptive,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
    )

    defaults(_) = (
        max_iter = 500,
        rho = 1e4,
        mu = 1e4,
        alpha1 = 1.0,
        alpha2 = 1.0,
        alpha3 = 1.0,
        rho_growth = 1.0,
        mu_growth = 1.0,
        adaptive = false,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-SOCWR SPLITBUS-HYBRID",
        test_dir = test_dir,
        case_list = ALADIN_SOCWR_SPLITBUS_HYBRID_CASES,
        model = Lirion.socwr(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_socwr,
        objective_dist = aladin_dist_obj_socwr,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = overrides,
        rtol = rtol,
        atol = atol,
    )
end
