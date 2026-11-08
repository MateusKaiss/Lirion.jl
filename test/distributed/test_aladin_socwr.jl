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

# N-way SplitBus: hub bus 1 of this 6-bus case borders zones 2 and 3, so it is a
# boundary bus shared by three zones. Impedances/reactive headroom are tuned to
# keep the SOCWR relaxation well-conditioned.
const ALADIN_SOCWR_SPLITBUS_3ZONE_CASES = [
    "case6_3zone_hub.m",
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
        # the socwr line search backtracks α₂ per iteration
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
# Hybrid consensus + flow proximal on the tie-lines. σ_pflow / σ_qflow
# stay small because they stack on σ_v2 / σ_r / σ_t on edge endpoints.
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

# ============================================================
# ALADIN SOCWR - :alg3_original (HFD Algorithm 3 verbatim, exact cone penalty)
# ------------------------------------------------------------
# (i)   the mode constructs for SOCWR (all couplings);
# (ii)  the conic distance used in Φ is the Euclidean distance to the cone:
#       closed forms vs hand-computed projections, degree-1 homogeneity, and
#       the RSOC ↔ SOC 45° isometry;
# (iii) on a live SOCWR iterate the exact penalty is degree-1 in the cone
#       vector while the `:alg3` squared residual is degree-2 (the point of
#       the fix), and y itself is cone-feasible (zero penalty);
# (iv)  the Danskin slope used by the Wolfe branch (c) matches a central
#       finite difference of V_ρ along λ_QP − λ;
# (v)   case3 (state/power/hybrid) converges under :alg3_original and every
#       step is one of HFD's patterns (1,1,1) / (1,0,0) / (0,0,α₃∈[0,1]).
# ============================================================
function _aladin_socwr_alg3_original_ws(path; coupling, rho, mu, max_iter)
    alg = Lirion.aladin(; coupling, sectioning = :tie_line, max_iter, rho, mu,
        adaptive = false, rho_growth = 1.0, mu_growth = 1.0,
        globalization = :alg3_original)
    inputs = Lirion._build_inputs(path, alg; model = Lirion.socwr(),
        optimizer = nothing, outputs = false)
    cfg = inputs.model_type
    of = Lirion.resolve_optimizer(cfg, inputs.optimizer)
    tag = Lirion.ALADIN{Lirion._variant_type(Val(cfg.coupling)), cfg.base_model}()
    return Lirion.build_workspace(inputs.system, cfg, of, tag), Lirion.build_params(cfg)
end

function run_aladin_socwr_alg3_original_tests()
    test_dir = @__DIR__
    case3 = joinpath(test_dir, "pglib_opf_case3_lmbd.m")

    @testset "ALADIN-SOCWR :alg3_original (HFD Alg. 3, exact cone penalty)" begin
        # (i) constructs
        for cpl in (:state, :power, :hybrid)
            cfg = Lirion.ALADINConfig(Lirion.SOCWR, cpl, Lirion.TieLine(), 10;
                globalization = :alg3_original)
            @test cfg.globalization === :alg3_original
        end

        # (ii) cone distance closed forms
        @testset "dist(·, K) closed forms" begin
            s2 = sqrt(2.0)
            @test Lirion._soc_dist([2.0, 1.0, 1.0]) == 0.0            # inside
            @test Lirion._soc_dist([0.0, 3.0, 4.0]) ≈ 5.0 / s2        # boundary-projected
            @test Lirion._soc_dist([-6.0, 3.0, 4.0]) ≈ sqrt(36 + 25)  # polar → apex
            # RSOC {2tu ≥ ‖x‖², t,u ≥ 0}: (1,1,x) with ‖x‖ = √2 is on the boundary
            @test Lirion._rsoc_dist([1.0, 1.0, 1.0, 1.0]) ≈ 0.0 atol = 1e-14
            @test Lirion._rsoc_dist([1.0, 1.0, 2.0, 0.0]) ≈
                  Lirion._soc_dist([2.0 / s2, 0.0, 2.0, 0.0])          # 45° isometry
            v = [0.3, 0.7, 1.9, -0.4]
            d = Lirion._rsoc_dist(v)
            @test d > 0
            for tau in (0.5, 2.0, 7.0)
                @test Lirion._rsoc_dist(tau .* v) ≈ tau * d rtol = 1e-12  # degree 1
                @test Lirion._soc_dist(tau .* [0.0, 3.0, 4.0]) ≈ tau * 5.0 / s2 rtol = 1e-12
            end
        end

        # (iii)+(iv) live iterate
        ws, params = _aladin_socwr_alg3_original_ws(case3;
            coupling = :state, rho = 1e4, mu = 1e4, max_iter = 500)
        ws.iter = 1
        redirect_stdout(devnull) do
            Lirion.local_solve!(ws)
            Lirion.solve_coordinator_qp!(ws)
            Lirion.compute_residuals!(ws)
            Lirion.compute_tolerances!(ws, params.eps_abs, params.eps_rel)
        end
        accy = (z, c) -> z.y[c]
        @test Lirion._socwr_cone_dist_l1(ws, accy) < 1e-6         # y cone-feasible
        rset = [Set(z.r_cols) for z in ws.zones]
        zidx = Dict(z => i for (i, z) in enumerate(ws.zones))
        acc = (z, c) -> (c in rset[zidx[z]] ? z.y[c] + 0.05 : z.y[c])
        # cone vectors at the violating point, scaled: exact = degree 1, :alg3 = degree 2
        zone = ws.zones[1]
        v2c = zone.var_to_col
        valfn = v -> acc(zone, v2c[v])
        found = false
        for (Fc, S) in Lirion.JuMP.list_of_constraint_types(zone.local_model)
            (S <: Lirion.MOI.RotatedSecondOrderCone) || continue
            for ci in Lirion.JuMP.all_constraints(zone.local_model, Fc, S)
                f = Float64[Lirion.JuMP.value(valfn, g) for g in Lirion.JuMP.constraint_object(ci).func]
                d = Lirion._rsoc_dist(f)
                d > 1e-9 || continue
                q = f[3]^2 + f[4]^2 - 2 * f[1] * f[2]
                g3 = 3.0 .* f
                @test Lirion._rsoc_dist(g3) ≈ 3.0 * d rtol = 1e-12
                @test (g3[3]^2 + g3[4]^2 - 2 * g3[1] * g3[2]) ≈ 9.0 * q rtol = 1e-12
                found = true
                break
            end
            found && break
        end
        @test found                                                # a violated cone existed
        @test Lirion._socwr_cone_dist_l1(ws, acc) > 0
        # (iv) Danskin slope vs central finite difference
        states = collect(values(ws.vars))
        saved = [copy(s.lambda) for s in states]
        dirs = map(s -> s.nu .- s.lambda, ws.vars)
        dvec = collect(values(dirs))
        setl!(a) = (for (s, s0, dd) in zip(states, saved, dvec); @. s.lambda = s0 + a * dd; end)
        V(a) = (setl!(a); redirect_stdout(devnull) do; Lirion._vrho_resolve_sum!(ws); end)
        V(0.5)
        g = Lirion._dual_slope(ws, dirs)
        h = 1e-4
        fd = (V(0.5 + h) - V(0.5 - h)) / (2h)
        for (s, s0) in zip(states, saved); copyto!(s.lambda, s0); end
        @test isapprox(g, fd; rtol = 1e-3, atol = 1e-6)

        # (iv-b) the Wolfe zoom itself: overshoot the dual direction ×4 so the
        # concave φ(α) = V_ρ(x, λ + α·4d) peaks at α ≈ 1/4 — α₃ = 1 must be
        # rejected and an interior point returned (strong Wolfe, or the best
        # Armijo/ascent point if the eval budget runs out); λ restored. The
        # shared `_branch_c_dual_linesearch_wolfe!` returns a status Symbol and
        # writes α₃ into ws.alpha3 (ws.alpha1 = ws.alpha2 = 0).
        for s in states
            @. s.nu = s.lambda + 4.0 * (s.nu - s.lambda)
        end
        saved4 = [copy(s.lambda) for s in states]
        dvec4 = [4.0 .* dd for dd in dvec]
        setl4!(a) = (for (s, s0, dd) in zip(states, saved4, dvec4); @. s.lambda = s0 + a * dd; end)
        V4(a) = (setl4!(a); redirect_stdout(devnull) do; Lirion._vrho_resolve_sum!(ws); end)
        st = redirect_stdout(devnull) do
            Lirion._branch_c_dual_linesearch_wolfe!(ws)
        end
        @test st in (:wolfe, :armijo_only, :ascent_only)
        a3 = ws.alpha3
        @test 0.0 < a3 < 1.0
        @test (ws.alpha1, ws.alpha2) == (0.0, 0.0)
        @test all(s.lambda == s0 for (s, s0) in zip(states, saved4))
        @test V4(a3) > V4(0.0)                        # V_ρ strictly increased
        for (s, s0) in zip(states, saved); copyto!(s.lambda, s0); end

        # (v) case14 (paper tune): state + hybrid converge under :alg3_original,
        # power reaches the centralized objective within 1 %; every step is an
        # HFD pattern (1,1,1) / (1,0,0) / (0,0,α₃ ∈ [0,1]).
        case14 = joinpath(test_dir, "pglib_opf_case14_ieee.m")
        central14 = 2175.683887902792     # Clarabel SOCWR centralized (runtests)
        for (cpl, rho, mu, need_conv) in ((:state, 100.0, 100.0, true),
                                          (:power, 100.0, 10.0, false),
                                          (:hybrid, 100.0, 10.0, true))
            ws, params = _aladin_socwr_alg3_original_ws(case14;
                coupling = cpl, rho, mu, max_iter = 1000)
            conv = false
            ok_pattern = true
            n_c = 0
            redirect_stdout(devnull) do
                a3cfg = ws.alpha3
                for iter in 1:params.max_iter
                    ws.iter = iter
                    Lirion.local_solve!(ws)
                    Lirion.solve_coordinator_qp!(ws)
                    Lirion.compute_residuals!(ws)
                    Lirion.compute_tolerances!(ws, params.eps_abs, params.eps_rel)
                    Lirion._record_iter!(ws)
                    if Lirion.stop_if_converged(ws, params)
                        conv = true
                        break
                    end
                    ws.alpha1 = 1.0
                    ws.alpha3 = a3cfg
                    Lirion._line_search_step!(ws)
                    a = (ws.alpha1, ws.alpha2, ws.alpha3)
                    ok_pattern &= a == (1.0, 1.0, 1.0) || a == (1.0, 0.0, 0.0) ||
                        (a[1] == 0.0 && a[2] == 0.0 && 0.0 <= a[3] <= 1.0)
                    ws.alg3_branch === :c && (n_c += 1)
                    Lirion.primal_update!(ws)
                    Lirion.dual_update!(ws)
                    Lirion.adapt_rho_mu!(ws, params)
                end
            end
            cost = sum(p.generation_cost for p in ws.partitions.parts)
            need_conv && @test conv
            @test ok_pattern
            @test n_c > 0                                  # branch (c) exercised
            @test abs(cost - central14) / central14 < 0.01
            @test ws.lam_bar > 0 && ws.kap_bar > 0
        end

        # (vi) paper-literal fixed λ̄ = κ̄ via alg3_lam_bar
        @test_throws ArgumentError Lirion.ALADINConfig(Lirion.SOCWR, :state,
            Lirion.TieLine(), 10; globalization = :alg3_original, alg3_lam_bar = -1.0)
        alg = Lirion.aladin(; coupling = :state, sectioning = :tie_line, max_iter = 1000,
            rho = 100.0, mu = 100.0, adaptive = false, rho_growth = 1.0, mu_growth = 1.0,
            globalization = :alg3_original, alg3_lam_bar = 2e4)
        inputs = Lirion._build_inputs(case14, alg; model = Lirion.socwr(),
            optimizer = nothing, outputs = false)
        _, ws, iters = redirect_stdout(devnull) do
            Lirion.solve_distributed(inputs)
        end
        @test ws.lam_bar == 2e4 && ws.kap_bar == 2e4
        cost = sum(p.generation_cost for p in ws.partitions.parts)
        @test abs(cost - central14) / central14 < 0.01

        # (vii) SOCWR knobs `_SOCWR_QP_PROX_SCALE` (QP-side ρσ scale; 0.0 =
        # prox out of B) and `_SOCWR_ACTIVE_{DUAL,CONE}_TOL`. Defaults
        # (1.0 / 1e-6 / 1e-6) reproduce the legacy diagonal and active set
        # bit-for-bit; the knobs are recorded in the manifest.
        @testset "SOCWR prox-in-B scale + active-set tolerance knobs" begin
            @test Lirion._SOCWR_QP_PROX_SCALE[] == 1.0
            @test Lirion._SOCWR_ACTIVE_DUAL_TOL[] == 1e-6
            @test Lirion._SOCWR_ACTIVE_CONE_TOL[] == 1e-6
            ws, params = _aladin_socwr_alg3_original_ws(case3;
                coupling = :state, rho = 1e4, mu = 1e4, max_iter = 5)
            ws.iter = 1
            redirect_stdout(devnull) do
                Lirion.local_solve!(ws)
            end
            legacy = [copy(z.hess_diag) for z in ws.zones]
            C_legacy = [Lirion.build_active_jacobian(z) for z in ws.zones]
            try
                Lirion._SOCWR_QP_PROX_SCALE[] = 0.0
                for z in ws.zones
                    Lirion.refresh_lagrangian_hessian!(z, ws)
                    @test z.hess_diag == max.(z.base_hess_diag, 1e-6)  # no ρσ
                end
                @test Lirion.aladin_runtime_manifest()["SOCWR_QP_PROX_SCALE"] == 0.0
                # legacy scale restores the legacy diagonal bit-for-bit
                Lirion._SOCWR_QP_PROX_SCALE[] = 1.0
                for (z, h) in zip(ws.zones, legacy)
                    Lirion.refresh_lagrangian_hessian!(z, ws)
                    @test z.hess_diag == h
                end
                @test any(z -> any(z.hess_diag .> max.(z.base_hess_diag, 1e-6)), ws.zones)
                # activity thresholds: defaults ≡ the old hard-coded 1e-6; a
                # looser cone threshold can only drop cone rows
                for (z, C) in zip(ws.zones, C_legacy)
                    @test Lirion.build_active_jacobian(z; dual_tol = 1e-6, cone_tol = 1e-6) == C
                end
                Lirion._SOCWR_ACTIVE_CONE_TOL[] = Inf
                for (z, C) in zip(ws.zones, C_legacy)
                    @test size(Lirion.build_active_jacobian(z), 1) <= size(C, 1)
                end
                @test Lirion.aladin_runtime_manifest()["SOCWR_ACTIVE_CONE_TOL"] == Inf
                Lirion._SOCWR_ACTIVE_CONE_TOL[] = 1e-6
                # canonical capability: a full SOCWR :alg3_original run with the
                # prox OUT of the QP Hessian executes end to end
                Lirion._SOCWR_QP_PROX_SCALE[] = 0.0
                alg = Lirion.aladin(; coupling = :state, sectioning = :tie_line,
                    max_iter = 30, rho = 1e4, mu = 1e4, adaptive = false,
                    rho_growth = 1.0, mu_growth = 1.0, globalization = :alg3_original)
                inputs = Lirion._build_inputs(case3, alg; model = Lirion.socwr(),
                    optimizer = nothing, outputs = false)
                _, wsc, iters = redirect_stdout(devnull) do
                    Lirion.solve_distributed(inputs)
                end
                @test 1 <= iters <= 30
                @test wsc.hist.meta[:terminal_status] in (:converged, :iteration_limit, :stalled)
                @test Lirion.aladin_run_manifest(wsc, inputs.model_type)["GLOBALIZATION"] ==
                      "alg3_original"
            finally
                Lirion._SOCWR_QP_PROX_SCALE[] = 1.0
                Lirion._SOCWR_ACTIVE_DUAL_TOL[] = 1e-6
                Lirion._SOCWR_ACTIVE_CONE_TOL[] = 1e-6
            end
        end
    end
end

# ============================================================
# ALADIN SOCWR - SPLIT-BUS HYBRID, N-WAY (boundary bus shared by 3 zones)
# ============================================================
function run_aladin_socwr_splitbus_3zone_tests(;
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
        max_iter = 8000,
        rho = 10000.0,
        mu = 1000.0,
        eps_abs = 2e-4,
        eps_rel = 2e-4,
    )

    return run_case_suite(;
        label = "ALADIN-SOCWR SPLITBUS 3-ZONE",
        test_dir = test_dir,
        case_list = ALADIN_SOCWR_SPLITBUS_3ZONE_CASES,
        model = Lirion.socwr(),
        central_optimizer = central_optimizer,
        distributed_optimizer = distributed_optimizer,
        objective_central = aladin_central_obj_socwr,
        objective_dist = aladin_dist_obj_socwr,
        make_algorithm = make_alg,
        defaults = defaults,
        overrides = _ -> NamedTuple(),
        rtol = rtol,
        atol = atol,
    )
end
