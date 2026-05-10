# ============================================================
# =====================  ALADIN SOCWR TESTS  ====================
# ============================================================

const ALADIN_SOCWR_STATE_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

const ALADIN_SOCWR_POWER_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

const ALADIN_SOCWR_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

const ALADIN_SOCWR_SPLITBUS_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

# ============================================================
# OBJECTIVE EXTRACTORS
# ============================================================
aladin_central_obj_socwr(out) = Float64(out.objective)
aladin_dist_obj_socwr(out) = Float64(out.objective)

# ============================================================
# ALADIN SOCWR - STATE
# ============================================================
function run_aladin_socwr_state_tests(;
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

    # case14_ieee: default rho/mu diverges. Paper tune uses rho=100,
    # mu=100, adaptive=false — gap=-0.228%, iters=612.
    overrides(case_name) = begin
        if case_name == "pglib_opf_case14_ieee.m"
            return (rho = 100.0, mu = 100.0, max_iter = 1000)
        end
        return NamedTuple()
    end

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
    test_dir = @__DIR__

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

    # case14_ieee: default rho/mu diverges. Paper tune uses rho=100, mu=10
    # — gap=-0.001%, iters=293.
    overrides(case_name) = begin
        if case_name == "pglib_opf_case14_ieee.m"
            return (rho = 100.0, mu = 10.0, max_iter = 1000)
        end
        return NamedTuple()
    end

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
    test_dir = @__DIR__

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

    # case14_ieee: default rho/mu diverges. Paper tune uses rho=100, mu=10
    # — gap=+0.008%, iters=347.
    overrides(case_name) = begin
        if case_name == "pglib_opf_case14_ieee.m"
            return (rho = 100.0, mu = 10.0, max_iter = 1000)
        end
        return NamedTuple()
    end

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
# ALADIN SOCWR - HYBRID + proximal_mode=:pflow
# ------------------------------------------------------------
# Same consensus as the default Hybrid (shared v² + tie-line p/q_flow
# + r/t); the prox additionally regularises the flow expressions on
# the tie-lines. σ_pflow / σ_qflow stay small because in Hybrid the
# pflow prox stacks on top of σ_v2 / σ_r / σ_t on edge endpoints.
# ============================================================
function run_aladin_socwr_hybrid_pflow_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = @__DIR__

    make_alg(p) = Lirion.aladin(;
        coupling = :hybrid,
        sectioning = :tie_line,
        proximal_mode = :pflow,
        max_iter = p.max_iter,
        rho = p.rho,
        mu = p.mu,
        alpha1 = p.alpha1,
        alpha2 = p.alpha2,
        alpha3 = p.alpha3,
        rho_growth = p.rho_growth,
        mu_growth = p.mu_growth,
        adaptive = p.adaptive,
        sigma_pflow = p.sigma_pflow,
        sigma_qflow = p.sigma_qflow,
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
        sigma_pflow = 1.0,
        sigma_qflow = 1.0,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
    )

    # case14_ieee inherits the same paper tune as the primal Hybrid path.
    overrides(case_name) = begin
        if case_name == "pglib_opf_case14_ieee.m"
            return (rho = 100.0, mu = 10.0, max_iter = 1000)
        end
        return NamedTuple()
    end

    return run_case_suite(;
        label = "ALADIN-SOCWR HYBRID PFLOW",
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
    test_dir = @__DIR__

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
