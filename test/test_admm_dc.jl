# ============================================================
# =========================  DC TESTS  =========================
# ============================================================

# ------------------------
# DC - STATE (bus/state coupling)
# ------------------------
const DC_STATE_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
    "pglib_opf_case73_ieee_rts.m",
    "pglib_opf_case118_ieee.m",
    "pglib_opf_case57_ieee.m",

    # "pglib_opf_case200_activ.m",
    # "pglib_opf_case500_goc.m",
]

# ------------------------
# DC - POWER (branch/power coupling)
# ------------------------
const DC_POWER_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
    "pglib_opf_case73_ieee_rts.m",
    "pglib_opf_case57_ieee.m",
    "pglib_opf_case118_ieee.m",

    # "pglib_opf_case200_activ.m",
    # "pglib_opf_case500_goc.m",
]

# ------------------------
# DC - HYBRID
# ------------------------
const DC_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case3_lmbd__api.m",
    "pglib_opf_case3_lmbd__sad.m",
    "pglib_opf_case73_ieee_rts.m",
    "pglib_opf_case57_ieee.m",
    "pglib_opf_case200_activ.m",

    # "pglib_opf_case118_ieee.m",
    # "pglib_opf_case500_goc.m",
]

# ============================================================
# OBJECTIVE EXTRACTORS
# ============================================================
central_obj_dc(out) = Float64(out.result["objective"])
dist_obj_dc(out) = Float64(out.result.total_real_cost)

# ============================================================
# DC - STATE
# ============================================================
function run_admm_dc_state_tests(;
    rtol = 0.01,
    atol = 1e-6,
    central_optimizer = nothing,
    distributed_optimizer = nothing,
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
        if case_name == "pglib_opf_case73_ieee_rts.m" ||
           case_name == "pglib_opf_case57_ieee.m" ||
           case_name == "pglib_opf_case118_ieee.m" ||
           case_name == "pglib_opf_case200_activ.m"
            return (
                rho = 100000.0,
                max_iter = 100000,
                eps_abs = 1e-4,
                eps_rel = 1e-4,
                adaptive = false,
                central_optimizer = Ipopt.Optimizer,
                distributed_optimizer = Ipopt.Optimizer,
            )
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
        rho = 100000.0,
        eps_abs = 1e-4,
        eps_rel = 1e-4,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    overrides(case_name) = begin
        if case_name == "pglib_opf_case500_goc.m"
            return (
                rho = 5000.0,
                distributed_optimizer = Ipopt.Optimizer,
            )
        end
        return NamedTuple()
    end

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
        rho = 10000.0,
        eps_abs = 2e-4,
        eps_rel = 2e-4,
        adaptive = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    )

    overrides(case_name) = begin
        if case_name == "pglib_opf_case500_goc.m"
            return (
                rho = 5000.0,
                central_optimizer = Ipopt.Optimizer,
                distributed_optimizer = Ipopt.Optimizer,
            )
        end
        return NamedTuple()
    end

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
