# ------------------------
# DC - STATE (bus coupling)
# ------------------------
const AC_STATE_CASES = [
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
# AC - POWER (branch coupling)
# ------------------------
const AC_POWER_CASES = [
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

# ============================================================
# OBJECTIVE EXTRACTORS
# ============================================================
central_obj_ac(out) = Float64(out.result["objective"])
dist_obj_ac(out) = Float64(out.result.total_expr_cost)

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
        rho = 1.0,
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
    overrides(case_name) = begin
        if case_name == "pglib_opf_case200_activ.m" || case_name == "pglib_opf_case500_goc.m"
            return (
                rho = 1.0e7,
                max_iter = 100000,
                eps_abs = 1e-4,
                eps_rel = 1e-4,
                adaptive = true,
                mu_incr = 10.0,
                mu_decr = 10.0,
                tau_incr = 1.2,
                tau_decr = 0.8,
                central_optimizer = Ipopt.Optimizer,
                distributed_optimizer = Ipopt.Optimizer,
            )
        elseif case_name == "pglib_opf_case73_ieee_rts.m" ||
               case_name == "pglib_opf_case57_ieee.m" ||
               case_name == "pglib_opf_case118_ieee.m"
            return (
                rho = 1.0e7,
                max_iter = 100000,
                eps_abs = 1e-4,
                eps_rel = 1e-4,
                adaptive = true,
                mu_incr = 10.0,
                mu_decr = 10.0,
                tau_incr = 1.05,
                tau_decr = 0.98,
                central_optimizer = Ipopt.Optimizer,
                distributed_optimizer = Ipopt.Optimizer,
            )
        end
        return NamedTuple()
    end

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
        overrides = overrides,
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
        rho = 1e4,
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
