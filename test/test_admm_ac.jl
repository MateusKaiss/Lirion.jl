# ------------------------
# DC - STATE (bus coupling)
# ------------------------
const AC_STATE_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
    # "pglib_opf_case73_ieee_rts.m",
    # "pglib_opf_case57_ieee.m",
    # "pglib_opf_case118_ieee.m",

    # "pglib_opf_case200_activ.m",
    # "pglib_opf_case500_goc.m",
]

# ------------------------
# AC - POWER (branch coupling)
# ------------------------
const AC_POWER_CASES = [
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
# DC - HYBRID
# ------------------------
const AC_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
    "pglib_opf_case73_ieee_rts.m",

    # "pglib_opf_case118_ieee.m",
    # "pglib_opf_case57_ieee.m",
    # "pglib_opf_case200_activ.m",
    # "pglib_opf_case500_goc.m",
]

# ------------------------
# AC - SPLIT-BUS HYBRID
# ------------------------
const AC_SPLITBUS_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__sad.m",
]

# ------------------------
# AC - DUMMY-BUS HYBRID
# ------------------------
const AC_DUMMYBUS_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__sad.m",
]

# ------------------------
# AC - DUMMY-BUS POWER (only P,Q consensus; θ_d, V_d ambiguous between zones)
# ------------------------
const AC_DUMMYBUS_POWER_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__sad.m",
]

# ============================================================
# OBJECTIVE EXTRACTORS
# ============================================================
central_obj_ac(out) = Float64(out.result["objective"])
dist_obj_ac(out) = Float64(out.result.total_cost)

# ============================================================
# AC - STATE
# ============================================================
function run_admm_ac_state_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = joinpath(@__DIR__, "distributed")

    # ------------------------
    # Algorithm factory
    # ------------------------
    make_alg(p) = Lirion.admm(;
        coupling = :state,
        sectioning = :tie_line,
        max_iter = Int(p.max_iter),
        rho = Float64(p.rho),
        eps_abs = Float64(p.eps_abs),
        eps_rel = Float64(p.eps_rel),
        adaptive_rho = Bool(p.adaptive),
        mu_incr = Float64(p.mu_incr),
        mu_decr = Float64(p.mu_decr),
        tau_incr = Float64(p.tau_incr),
        tau_decr = Float64(p.tau_decr),
    )

    defaults(_) = (
        max_iter = 1000,
        rho = 1000000.0,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
        adaptive = true,
        mu_incr = 2.0,
        mu_decr = 2.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    # ------------------------
    # Overrides por caso
    # ------------------------
    # overrides(case_name) = begin
    #     if case_name == "pglib_opf_case200_activ.m" || case_name == "pglib_opf_case500_goc.m"
    #         return (
    #             rho = 1.0e7,
    #             max_iter = 100000,
    #             eps_abs = 1e-4,
    #             eps_rel = 1e-4,
    #             adaptive = true,
    #             mu_incr = 10.0,
    #             mu_decr = 10.0,
    #             tau_incr = 1.2,
    #             tau_decr = 0.8,
    #             central_optimizer = Ipopt.Optimizer,
    #             distributed_optimizer = Ipopt.Optimizer,
    #         )
    #     elseif case_name == "pglib_opf_case73_ieee_rts.m" ||
    #            case_name == "pglib_opf_case57_ieee.m" ||
    #            case_name == "pglib_opf_case118_ieee.m"
    #         return (
    #             rho = 1.0e7,
    #             max_iter = 100000,
    #             eps_abs = 1e-4,
    #             eps_rel = 1e-4,
    #             # adaptive = true,
    #             # mu_incr = 10.0,
    #             # mu_decr = 10.0,
    #             # tau_incr = 1.05,
    #             # tau_decr = 0.98,
    #             central_optimizer = Ipopt.Optimizer,
    #             distributed_optimizer = Ipopt.Optimizer,
    #         )
    #     end
    #     return NamedTuple()
    # end

    return run_case_suite(;
        label = "ADMM-AC STATE",
        test_dir = test_dir,
        case_list = AC_STATE_CASES,
        model = Lirion.ac(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_ac,
        objective_dist = dist_obj_ac,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = _ -> NamedTuple(),
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# AC - POWER
# ============================================================
function run_admm_ac_power_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = joinpath(@__DIR__, "distributed")

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
        rho = 1e5,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    return run_case_suite(;
        label = "ADMM-AC POWER",
        test_dir = test_dir,
        case_list = AC_POWER_CASES,
        model = Lirion.ac(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_ac,
        objective_dist = dist_obj_ac,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = _ -> NamedTuple(),
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# AC - HYBRID
# ============================================================
function run_admm_ac_hybrid_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = joinpath(@__DIR__, "distributed")

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
        rho = 1e3,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    return run_case_suite(;
        label = "ADMM-AC HYBRID",
        test_dir = test_dir,
        case_list = AC_HYBRID_CASES,
        model = Lirion.ac(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_ac,
        objective_dist = dist_obj_ac,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = _ -> NamedTuple(),
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# AC - SPLIT-BUS HYBRID (Method 3 — Harris et al. 2022, Fig. 1c)
# ============================================================
function run_admm_ac_splitbus_hybrid_tests(;
    rtol = 0.05,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = joinpath(@__DIR__, "distributed")

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

    return run_case_suite(;
        label = "ADMM-AC SPLITBUS-HYBRID",
        test_dir = test_dir,
        case_list = AC_SPLITBUS_HYBRID_CASES,
        model = Lirion.ac(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_ac,
        objective_dist = dist_obj_ac,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = _ -> NamedTuple(),
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# AC - DUMMY-BUS HYBRID (Method 2 — Harris et al. 2022, Fig. 1b)
# ============================================================
function run_admm_ac_dummybus_hybrid_tests(;
    rtol = 0.05,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = joinpath(@__DIR__, "distributed")

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

    return run_case_suite(;
        label = "ADMM-AC DUMMYBUS-HYBRID",
        test_dir = test_dir,
        case_list = AC_DUMMYBUS_HYBRID_CASES,
        model = Lirion.ac(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_ac,
        objective_dist = dist_obj_ac,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = _ -> NamedTuple(),
        rtol = rtol,
        atol = atol,
    )
end

# ============================================================
# AC - DUMMY-BUS POWER (Method 2 with power-only coupling)
# ------------------------------------------------------------
# Theory: dummy bus d is artificial; θ_d, V_d gauge-invariant for cost.
# Enforcing P_{m,d} = P_{d,n} and Q_{m,d} = Q_{d,n} pins the physical
# transfer; artificial gens close the local KCLs. In AC (NLP) the method
# is empirically less stable than SOCWR — tightly-bounded variants can
# produce NaN or infeasibility in some iterations.
# ============================================================
function run_admm_ac_dummybus_power_tests(;
    rtol = 0.05,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
)
    test_dir = joinpath(@__DIR__, "distributed")

    make_alg(p) = Lirion.admm(;
        coupling = :power,
        sectioning = Lirion.DummyBus(),
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
        label = "ADMM-AC DUMMYBUS-POWER",
        test_dir = test_dir,
        case_list = AC_DUMMYBUS_POWER_CASES,
        model = Lirion.ac(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = central_obj_ac,
        objective_dist = dist_obj_ac,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = _ -> NamedTuple(),
        rtol = rtol,
        atol = atol,
    )
end
