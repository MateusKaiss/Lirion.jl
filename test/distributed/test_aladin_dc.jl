# ============================================================
# =======================  ALADIN DC TESTS  =====================
# ============================================================

# ------------------------
# ALADIN DC - STATE
# ------------------------
const ALADIN_DC_STATE_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

const ALADIN_DC_POWER_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

const ALADIN_DC_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

const ALADIN_DC_SPLITBUS_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

# N-way SplitBus: bus 1 of this 4-bus star borders zones 2 and 3, so it is a
# boundary bus shared by three zones (its own + 2 + 3).
const ALADIN_DC_SPLITBUS_3ZONE_CASES = [
    "case4_3zone_star.m",
]

const ALADIN_DC_DUMMYBUS_HYBRID_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
]

# Power-only DummyBus: case3_lmbd__api excluded (same issue as ADMM — without θ
# consensus at the dummy bus, tight bounds can drive early iterations to
# infeasibility).
const ALADIN_DC_DUMMYBUS_POWER_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case14_ieee.m",
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
# ALADIN DC - HYBRID + proximal_mode=:pflow
# ------------------------------------------------------------
# Hybrid consensus + flow proximal on the tie-lines. σ_pflow stays small
# because it adds to σ_θ_shared on edge-endpoint θ.
# ============================================================
function run_aladin_dc_hybrid_pflow_tests(;
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
        sigma_pflow = p.sigma_pflow,
        eps_abs = p.eps_abs,
        eps_rel = p.eps_rel,
    )

    defaults(_) = (
        max_iter = 5000,
        rho = 10000.0,
        mu = 1000.0,
        sigma_pflow = 1.0,
        eps_abs = 2e-4,
        eps_rel = 2e-4,
    )

    overrides(_case_name) = NamedTuple()

    return run_case_suite(;
        label = "ALADIN-DC HYBRID PFLOW",
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
# ALADIN DC - SPLIT-BUS HYBRID, N-WAY (boundary bus shared by 3 zones)
# ============================================================
function run_aladin_dc_splitbus_3zone_tests(;
    rtol = 0.02,
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

    return run_case_suite(;
        label = "ALADIN-DC SPLITBUS 3-ZONE",
        test_dir = test_dir,
        case_list = ALADIN_DC_SPLITBUS_3ZONE_CASES,
        model = Lirion.dc(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_dc,
        objective_dist = aladin_dist_obj_dc,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = _ -> NamedTuple(),
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

    # case14_ieee: rho=1e5 drives the coordinator QP to DUAL_INFEASIBLE at
    # the first step. Paper tune uses rho=1e3, mu=1e3, adaptive (iters=3).
    overrides(case_name) = begin
        if case_name == "pglib_opf_case14_ieee.m"
            return (rho = 1000.0, mu = 1000.0)
        end
        return NamedTuple()
    end

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

# ============================================================
# ALADIN DC - :alg3_original (HFD Algorithm 3 as published)
# ------------------------------------------------------------
# 1. constructs for DC; 2. the Danskin dual slope hook `_dual_slope`
#    (V_ρ'(α) = (ν−λ)ᵀ(Ay−b)) matches a central finite difference of V_ρ —
#    a wrong sign would silently turn every branch (c) into α₃ = 0;
# 3. the Wolfe branch (c) returns α₃ ∈ [0,1], restores λ and attains at
#    least the golden-section V_ρ; 4. full runs converge, produce only the
#    HFD step patterns (a) (1,1,1) / (b) (1,0,0) / (c) (0,0,α₃), keep λ̄ fixed,
#    keep consistent branch bookkeeping and reach the same optimum as :alg3.
#    (Trajectory identity with :alg3 is NOT asserted: on DC both modes reach
#    branch (c) and legitimately differ there — Wolfe with α₃ = 1 reachable
#    vs golden section; fixed vs ratcheted λ̄.)
# ============================================================
function _aladin_dc_alg3_original_inputs(case_file; coupling, sectioning, glob,
                                         rho, mu, max_iter)
    alg = Lirion.aladin(; coupling, sectioning, max_iter, rho, mu,
        adaptive = false, rho_growth = 1.0, mu_growth = 1.0,
        eps_abs = 1e-4, eps_rel = 1e-4, globalization = glob)
    return Lirion._build_inputs(case_file, alg; model = Lirion.dc(),
        optimizer = nothing, outputs = false)
end

function _aladin_dc_workspace_at_iter1(inputs)
    cfg = inputs.model_type
    fac = Lirion.resolve_optimizer(cfg, inputs.optimizer)
    tag = Lirion.ALADIN{Lirion._variant_type(Val(cfg.coupling)), cfg.base_model}()
    params = Lirion.build_params(cfg)
    ws = Lirion.build_workspace(inputs.system, cfg, fac, tag)
    ws.iter = 1
    Lirion.local_solve!(ws)
    Lirion.solve_coordinator_qp!(ws)
    Lirion.compute_residuals!(ws)
    Lirion.compute_tolerances!(ws, params.eps_abs, params.eps_rel)
    return ws
end

function _aladin_dc_vrho_at!(ws, saved, dir, alpha)
    for (s, s0, d) in zip(values(ws.vars), values(saved), values(dir))
        @. s.lambda = s0 + alpha * d
    end
    return Lirion._vrho_resolve_sum!(ws)
end

function run_aladin_dc_alg3_original_tests()
    test_dir = @__DIR__
    case3 = joinpath(test_dir, "pglib_opf_case3_lmbd.m")
    case14 = joinpath(test_dir, "pglib_opf_case14_ieee.m")

    @testset "ALADIN-DC :alg3_original" begin
        @testset "constructs for DC" begin
            for c in (:state, :power, :hybrid)
                cfg = Lirion.ALADINConfig(Lirion.DC, c, Lirion.TieLine(), 10;
                    globalization = :alg3_original)
                @test cfg.globalization === :alg3_original
            end
        end

        @testset "_dual_slope == d/dα V_ρ(x, λ+α(ν−λ)) (finite difference)" begin
            for (coupling, sectioning) in ((:state, :tie_line), (:power, :tie_line),
                                           (:hybrid, :tie_line), (:hybrid, :dummy_bus),
                                           (:hybrid, :split_bus))
                inputs = _aladin_dc_alg3_original_inputs(case14; coupling, sectioning,
                    glob = :alg3_original, rho = 1000.0, mu = 1000.0, max_iter = 5)
                ws = redirect_stdout(devnull) do
                    _aladin_dc_workspace_at_iter1(inputs)
                end
                saved = map(s -> copy(s.lambda), ws.vars)
                dir = map(s -> s.nu .- s.lambda, ws.vars)
                # re-solve at α = 0: at iteration 1 the models were solved with
                # the pre-copy proximal centre (x_prev := y happens after the
                # solve), so their objective values are NOT V_ρ(x, λ).
                V0 = _aladin_dc_vrho_at!(ws, saved, dir, 0.0)
                g0 = Lirion._dual_slope(ws, dir)
                h = 1e-4
                Vp = _aladin_dc_vrho_at!(ws, saved, dir, h)
                Vm = _aladin_dc_vrho_at!(ws, saved, dir, -h)
                V1 = _aladin_dc_vrho_at!(ws, saved, dir, 1.0)
                fd = (Vp - Vm) / (2h)
                # tolerances at OSQP (eps 1e-6) noise / active-set-kink level
                @test isapprox(g0, fd; rtol = 1e-2, atol = 1e-5 * max(1.0, abs(V0)))
                # concavity of V_ρ in λ: chord below tangent
                @test V1 - V0 <= g0 + 1e-2 * abs(g0) + 1e-5 * max(1.0, abs(V0))

                # Wolfe branch (c) vs golden branch (c) on the same state
                _aladin_dc_vrho_at!(ws, saved, dir, 0.0)
                Lirion._branch_c_dual_linesearch_wolfe!(ws)
                a_w = ws.alpha3
                @test ws.alpha1 == 0.0 && ws.alpha2 == 0.0 && 0.0 <= a_w <= 1.0
                @test all(s.lambda == s0 for (s, s0) in zip(values(ws.vars), values(saved)))
                V_w = _aladin_dc_vrho_at!(ws, saved, dir, a_w)
                _aladin_dc_vrho_at!(ws, saved, dir, 0.0)
                Lirion._branch_c_dual_linesearch!(ws)
                V_g = _aladin_dc_vrho_at!(ws, saved, dir, ws.alpha3)
                if g0 > 0
                    @test V_w >= V_g - 1e-6 * max(1.0, abs(V_g))
                else
                    @test a_w == 0.0      # not an ascent direction: honest stall
                end
            end
        end

        @testset "HFD α patterns, fixed λ̄, bookkeeping, optimum vs :alg3" begin
            for (file, coupling, rho, mu) in ((case3, :state, 1000.0, 1000.0),
                                              (case14, :state, 1000.0, 1000.0),
                                              (case14, :hybrid, 10000.0, 1000.0))
                res = Dict{Symbol, Any}()
                for glob in (:alg3, :alg3_original)
                    inputs = _aladin_dc_alg3_original_inputs(file; coupling,
                        sectioning = :tie_line, glob, rho, mu, max_iter = 5000)
                    _, ws, iters = redirect_stdout(devnull) do
                        Lirion.solve_distributed(inputs)
                    end
                    res[glob] = (ws = ws, iters = iters,
                                 cost = Lirion._total_generation_cost(ws.partitions))
                end
                o = res[:alg3_original]; a = res[:alg3]
                ws = o.ws
                conv = ws.r_primal < ws.eps_primal && ws.r_step_primal < ws.eps_step_primal
                @test conv
                al = (ws.alpha1, ws.alpha2, ws.alpha3)
                @test al == (1.0, 1.0, 1.0) || al == (1.0, 0.0, 0.0) ||
                      (al[1] == 0.0 && al[2] == 0.0 && 0.0 <= al[3] <= 1.0)
                # λ̄ fixed for the whole run (set once, recorded once)
                @test ws.lam_bar > 0.0 && ws.hist.meta[:alg3_lam_bar] == ws.lam_bar
                # branch bookkeeping: one of (a)/(b)/(c) per line search; the
                # line search runs at iterations 1..iters−1 when converged
                m = ws.hist.meta[:alg3_original]
                n_ls = conv ? o.iters - 1 : o.iters
                @test m[:a] + m[:b] + m[:c] == n_ls
                # every branch (c) lands in exactly one outcome key, the
                # Danskin-gradient fallback outcomes (:c_*_grad) included
                @test m[:c] == m[:c_full] + m[:c_wolfe] + m[:c_budget] + m[:c_no_ascent] +
                               m[:c_full_grad] + m[:c_wolfe_grad] + m[:c_budget_grad]
                # same optimum as :alg3 (both converge to tolerance 1e-4)
                @test o.cost ≈ a.cost rtol = 1e-3
            end
        end

        # ---- fixed λ̄, κ̄ majoriser, certification guard, terminal status ----
        @testset "cfg.alg3_lam_bar honoured verbatim on DC" begin
            alg = Lirion.aladin(; coupling = :state, sectioning = :tie_line,
                max_iter = 200, rho = 1000.0, mu = 1000.0, adaptive = false,
                rho_growth = 1.0, mu_growth = 1.0, globalization = :alg3_original,
                alg3_lam_bar = 77.0)
            inputs = Lirion._build_inputs(case3, alg; model = Lirion.dc(),
                optimizer = nothing, outputs = false)
            _, ws, iters = redirect_stdout(devnull) do
                Lirion.solve_distributed(inputs)
            end
            @test ws.lam_bar_fixed
            @test ws.lam_bar == 77.0                      # never re-scaled
            @test ws.hist.meta[:alg3_lam_bar] == 77.0    # recorded = applied
            @test ws.hist.meta[:alg3_lam_bar_fixed] === true
            @test haskey(ws.hist.meta, :alg3_original)   # counters materialised
            # automatic policy still records its own λ̄ (auto ≠ fixed)
            alg0 = Lirion.aladin(; coupling = :state, sectioning = :tie_line,
                max_iter = 200, rho = 1000.0, mu = 1000.0, adaptive = false,
                rho_growth = 1.0, mu_growth = 1.0, globalization = :alg3_original)
            inputs0 = Lirion._build_inputs(case3, alg0; model = Lirion.dc(),
                optimizer = nothing, outputs = false)
            _, ws0, _ = redirect_stdout(devnull) do
                Lirion.solve_distributed(inputs0)
            end
            @test !ws0.lam_bar_fixed
            @test ws0.hist.meta[:alg3_lam_bar_fixed] === false
            @test ws0.lam_bar > 0.0 && ws0.hist.meta[:alg3_lam_bar] == ws0.lam_bar
            # bit-identity guard: other modes never see cfg.alg3_lam_bar
            for g in (:alg3, :full_step, :armijo_l2)
                cfg = Lirion.ALADINConfig(Lirion.DC, :state, Lirion.TieLine(), 3;
                    globalization = g, alg3_lam_bar = 7.0)
                tag = Lirion.ALADIN{Lirion._variant_type(Val(cfg.coupling)), cfg.base_model}()
                wsg = Lirion.build_workspace(inputs.system, cfg,
                    Lirion.resolve_optimizer(cfg, nothing), tag)
                @test wsg.lam_bar == 0.0 && !wsg.lam_bar_fixed
            end
        end

        @testset "κ̄ majoriser includes variable-bound duals (DC)" begin
            inputs = _aladin_dc_alg3_original_inputs(case14; coupling = :state,
                sectioning = :tie_line, glob = :alg3_original,
                rho = 1000.0, mu = 1000.0, max_iter = 5)
            ws = redirect_stdout(devnull) do
                _aladin_dc_workspace_at_iter1(inputs)
            end
            m0 = Lirion._dc_max_local_dual(ws)
            m1 = Lirion._dc_max_local_dual(ws; include_bounds = true)
            @test m1 >= m0                       # a majoriser can only grow
            @test m0 == Lirion._dc_max_local_dual(ws; include_bounds = false)
            # at least one box bound is active on case14 (gen limits): the
            # bound-dual walk must actually see a nonzero multiplier
            maxb = 0.0
            bdual(ref) = try; abs(Lirion.JuMP.dual(ref)); catch; 0.0; end
            for zone in ws.zones, v in zone.all_vars
                Lirion.JuMP.has_lower_bound(v) &&
                    (maxb = max(maxb, bdual(Lirion.JuMP.LowerBoundRef(v))))
                Lirion.JuMP.has_upper_bound(v) &&
                    (maxb = max(maxb, bdual(Lirion.JuMP.UpperBoundRef(v))))
            end
            @test m1 == max(m0, maxb)
        end

        @testset "no certification at iteration 1 (vacuous r_step)" begin
            inputs = _aladin_dc_alg3_original_inputs(case3; coupling = :state,
                sectioning = :tie_line, glob = :armijo_l2,
                rho = 1000.0, mu = 1000.0, max_iter = 5)
            ws = redirect_stdout(devnull) do
                _aladin_dc_workspace_at_iter1(inputs)
            end
            params = Lirion.build_params(inputs.model_type)
            @test ws.r_step_primal == 0.0        # x_prev ← y makes it vacuous
            ws.r_primal = 0.0                    # force the primal test to pass
            ws.iter = 1
            @test !Lirion.stop_if_converged(ws, params)   # guard
            ws.iter = 2
            @test Lirion.stop_if_converged(ws, params)    # same residuals, iter ≥ 2
            ws.r_primal = Inf
            @test !Lirion.stop_if_converged(ws, params)
        end

        @testset "terminal status, ALMOST counters, manifest" begin
            @test Lirion._aladin_terminal_status(true, true) === :converged
            @test Lirion._aladin_terminal_status(true, false) === :converged
            @test Lirion._aladin_terminal_status(false, true) === :stalled
            @test Lirion._aladin_terminal_status(false, false) === :iteration_limit
            for (max_iter, expect) in ((2, :iteration_limit), (5000, :converged))
                inputs = _aladin_dc_alg3_original_inputs(case3; coupling = :state,
                    sectioning = :tie_line, glob = :alg3_original,
                    rho = 1000.0, mu = 1000.0, max_iter = max_iter)
                parts, ws, iters = redirect_stdout(devnull) do
                    Lirion.solve_distributed(inputs)
                end
                meta = ws.hist.meta
                @test meta[:terminal_status] === expect
                @test meta[:zero_step_iters] >= 0
                @test length(meta[:almost_zone_solves]) == iters
                @test length(meta[:almost_qp_solves]) == iters
                @test all(>=(0), meta[:almost_zone_solves])
                @test all(>=(0), meta[:almost_qp_solves])
                r = redirect_stdout(devnull) do
                    Lirion.build_lirion_result_distributed(
                        parts, ws, inputs.model_type, iters, time_ns(), inputs.system)
                end
                @test r.status === expect
                if expect === :iteration_limit
                    # :stalled is a relabel of :iteration_limit only
                    meta[:terminal_status] = :stalled
                    r2 = redirect_stdout(devnull) do
                        Lirion.build_lirion_result_distributed(
                            parts, ws, inputs.model_type, iters, time_ns(), inputs.system)
                    end
                    @test r2.status === :stalled
                    meta[:terminal_status] = :iteration_limit
                end
                rm = Lirion.aladin_run_manifest(ws, inputs.model_type)
                @test rm["GLOBALIZATION"] == "alg3_original"
                @test rm["TERMINAL_STATUS"] == String(expect)
                @test rm["LAM_BAR_APPLIED"] == ws.lam_bar
                @test rm["ALMOST_ZONE_SOLVES"] == sum(meta[:almost_zone_solves])
                @test rm["ALMOST_QP_SOLVES"] == sum(meta[:almost_qp_solves])
                # λ̄ violation counter (auto λ̄ → 0), effective Σ (DC wires
                # cfg.sigma_* → defaults 10/1/1), local-model state — and the
                # per-zone solutions the postprocessing reads are zone.y
                @test rm["LAM_BAR_VIOLATIONS"] == 0
                @test rm["SIGMA_EFFECTIVE"]["sigma_theta_shared"] == ws.sigma_theta_shared == 10.0
                @test rm["SIGMA_EFFECTIVE"]["sigma_pg"] == 1.0
                @test rm["PROXIMAL_MODE"] == "primal"
                @test rm["LOCAL_MODELS"] == String(meta[:local_models])
                @test meta[:local_models] in (:accepted, :restored)
                expect === :converged && @test meta[:local_models] === :accepted
                for zone in ws.zones
                    @test maximum(abs.(Lirion.JuMP.value.(zone.all_vars) .- zone.y)) <= 1e-5
                end
            end
            m = Lirion.aladin_runtime_manifest()
            for key in ("AC_QP_NO_PROX", "SOCWR_QP_PROX_SCALE", "SOCWR_ACTIVE_DUAL_TOL",
                        "SOCWR_ACTIVE_CONE_TOL", "MODIFIED_GRADIENT",
                        "ALG3_PHI_LOCAL_VIOL", "ALG3_ORIG_SAFETY", "ALG3_ORIG_FLOOR",
                        "BOUND_DETECT_MODE", "QP_HESS_REG", "GAUGE_REG_MODE",
                        "ACTIVE_DUAL_TOL", "LIRION_GAUGEFIX")
                @test haskey(m, key)
            end
            @test m["AC_QP_NO_PROX"] == Lirion._AC_QP_NO_PROX[]
            @test m["SOCWR_QP_PROX_SCALE"] == 1.0
            @test m["MODIFIED_GRADIENT"] == Lirion._MODIFIED_GRADIENT[]
        end

        # A branch-(c) V_ρ probe leaves the JuMP local models at a trial λ;
        # after the loop they must be put back on the accepted iterate
        # (`_restore_local_models!`). Iteration 2 is driven by hand (iteration
        # 1 has the x_prev ← y copy).
        @testset "local models back on the accepted iterate after a (c) probe" begin
            alg = Lirion.aladin(; coupling = :state, sectioning = :tie_line,
                max_iter = 10, rho = 1000.0, mu = 1000.0, adaptive = true,
                rho_growth = 2.0, mu_growth = 1.0, rho_upper = 1e9, mu_upper = 1e9,
                eps_abs = 1e-4, eps_rel = 1e-4, globalization = :alg3_original)
            inputs = Lirion._build_inputs(case3, alg; model = Lirion.dc(),
                optimizer = nothing, outputs = false)
            cfg = inputs.model_type
            params = Lirion.build_params(cfg)
            ws = redirect_stdout(devnull) do
                _aladin_dc_workspace_at_iter1(inputs)
            end
            rho0 = ws.rho[]
            vals(w) = [Lirion.JuMP.value.(z.all_vars) for z in w.zones]
            dist(a, b) = maximum(maximum(abs.(a[i] .- b[i])) for i in eachindex(a))
            # finish iteration 1 exactly as run_algorithm! does (ρ doubles)
            ws.alpha1 = 1.0; ws.alpha3 = cfg.alpha3
            redirect_stdout(devnull) do
                Lirion._line_search_step!(ws)
            end
            Lirion.primal_update!(ws); Lirion.dual_update!(ws)
            Lirion.adapt_rho_mu!(ws, params)
            @test ws.rho[] == 2 * rho0
            # iteration 2 up to the residual record: models = accepted solve
            ws.iter = 2
            redirect_stdout(devnull) do
                Lirion.local_solve!(ws)
                ws.hist.meta[:local_models_probed] = false
                Lirion.solve_coordinator_qp!(ws)
                Lirion.compute_residuals!(ws)
                Lirion.compute_tolerances!(ws, params.eps_abs, params.eps_rel)
                Lirion._record_iter!(ws)
            end
            @test ws.hist.rhos[:rho][end] == 2 * rho0
            y_acc = [copy(z.y) for z in ws.zones]
            x_acc = [copy(z.x_prev) for z in ws.zones]
            @test dist(vals(ws), y_acc) == 0.0          # zone.y came from these values
            # a branch-(c) probe: re-solve at a trial λ — models move, zone.y not
            lam_acc = [copy(s.lambda) for s in values(ws.vars)]
            for s in values(ws.vars)
                s.lambda .+= 50.0
            end
            redirect_stdout(devnull) do
                Lirion._vrho_resolve_sum!(ws)
            end
            for (s, l) in zip(values(ws.vars), lam_acc)
                copyto!(s.lambda, l)
            end
            @test ws.hist.meta[:local_models_probed] === true
            @test dist(vals(ws), y_acc) > 1e-4           # probe point ≠ accepted
            @test all(z.y == y0 for (z, y0) in zip(ws.zones, y_acc))
            # end of a branch-(c) iteration: α = (0, 0, α₃), then ρ grows again
            ws.alpha1 = 0.0; ws.alpha2 = 0.0; ws.alpha3 = 0.5
            Lirion.primal_update!(ws); Lirion.dual_update!(ws)
            Lirion.adapt_rho_mu!(ws, params)
            @test ws.rho[] == 4 * rho0
            lam_end = [copy(s.lambda) for s in values(ws.vars)]
            @test all(s.lambda_prev == l for (s, l) in zip(values(ws.vars), lam_acc))
            outcome = redirect_stdout(devnull) do
                Lirion._finalize_local_models!(ws)
            end
            @test outcome === :restored
            @test ws.hist.meta[:local_models] === :restored
            # per-zone solutions are back on the accepted iterate (same QP,
            # re-solved at ρ = 2ρ₀ — OSQP tolerance 1e-6)
            @test dist(vals(ws), y_acc) <= 1e-5
            # ...and nothing the loop owns moved: x, λ, λ_prev, ρ, zone.y
            @test all(z.x_prev == x0 for (z, x0) in zip(ws.zones, x_acc))
            @test all(s.lambda == l for (s, l) in zip(values(ws.vars), lam_end))
            @test all(s.lambda_prev == l for (s, l) in zip(values(ws.vars), lam_acc))
            @test ws.rho[] == 4 * rho0
            @test all(z.y == y0 for (z, y0) in zip(ws.zones, y_acc))
            @test Lirion.aladin_run_manifest(ws, cfg)["LOCAL_MODELS"] == "restored"
            # guard: x moved (α₁ ≠ 0) → no restore point → :probe (flagged)
            ws.hist.meta[:local_models_probed] = true
            ws.alpha1 = 1.0
            @test Lirion._finalize_local_models!(ws) === :probe
            @test Lirion.aladin_run_manifest(ws, cfg)["LOCAL_MODELS"] == "probe"
            ws.alpha1 = 0.0
            # no probe in the final iteration → :accepted, no solve at all
            ws.hist.meta[:local_models_probed] = false
            @test Lirion._finalize_local_models!(ws) === :accepted
        end
    end
end
