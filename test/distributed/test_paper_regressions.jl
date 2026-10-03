# ============================================================
# ==============  PAPER REGRESSION TESTS  ====================
# ============================================================
# Measurement-layer behaviours the paper results depend on:
#   1. derived stationarity ruler: eps_step_primal == n_prox · eps_abs
#   2. structural mirror: _extra_prox_support! covers exactly the columns
#      _add_extra_prox_grad! touches (AC :pflow flow-proximal)
#   3. :alg3 globalization runs and only ever produces the HFD step
#      patterns (a) α=(1,1,1), (b) α=(1,0,0), (c) α=(0,0,α₃∈[0,1])
#   4. globalization availability across base models
#   5. :alg3_original (HFD Alg. 3 as printed, AC): registration, λ̄ policy
#      (monotone / fixed), HFD step patterns, branch-(c) Wolfe status log
# ============================================================

function _paper_reg_solve(alg; model)
    case3 = joinpath(@__DIR__, "pglib_opf_case3_lmbd.m")
    inputs = Lirion._build_inputs(case3, alg; model = model,
        optimizer = nothing, outputs = false)
    _, ws, iters = redirect_stdout(devnull) do
        Lirion.solve_distributed(inputs)
    end
    return ws, iters
end

function run_paper_regression_tests()
    @testset "Paper regressions (derived ruler / prox mirror / alg3 / guards)" begin
        # ---- 1. derived stationarity ruler --------------------------
        @testset "eps_step_primal == n_prox * eps_abs" begin
            alg = Lirion.aladin(; coupling = :state, sectioning = :tie_line,
                max_iter = 5, rho = 100.0, mu = 100.0, adaptive = false,
                rho_growth = 1.0, mu_growth = 1.0,
                eps_abs = 1e-4, eps_rel = 1e-4)
            ws, _ = _paper_reg_solve(alg; model = Lirion.ac())
            np = Lirion._structural_nprox(ws)
            @test np > 0
            @test ws.eps_step_primal ≈ np * 1e-4
        end

        # ---- 2. structural mirror of the AC :pflow prox gradient ----
        @testset "_extra_prox_support! mirrors _add_extra_prox_grad!" begin
            alg = Lirion.aladin(; coupling = :hybrid, sectioning = :tie_line,
                max_iter = 3, rho = 100.0, mu = 100.0, adaptive = false,
                rho_growth = 1.0, mu_growth = 1.0, proximal_mode = :pflow)
            ws, _ = _paper_reg_solve(alg; model = Lirion.ac())
            n_support_total = 0
            for (pos, zone) in enumerate(ws.zones)
                sup = Set{Int}()
                Lirion._extra_prox_support!(sup, ws, pos, zone)
                n_support_total += length(sup)
                # deterministic generic perturbation (no Random dep in the test
                # target): flow(y) − flow_prev ≠ 0 and every linearisation
                # derivative nonzero (measure-zero exceptions)
                zone.y .+= 0.01 .* (1.2 .+ sin.((1.0:zone.n_vars) .+ pos))
                grad = zeros(zone.n_vars)
                Lirion._add_extra_prox_grad!(grad, ws, pos, zone)
                touched = Set(findall(!iszero, grad))
                @test touched == sup
            end
            @test n_support_total > 0   # pflow hook must be active here
        end

        # ---- 3. :alg3 runs (AC) and yields only HFD step patterns ---
        # α is sampled at EVERY iteration 1..8 by capping max_iter successively
        # (the workspace exposes only the last step's α; the caps walk the same
        # deterministic trajectory, so each run observes iteration k's α).
        @testset ":alg3 produces α ∈ {(1,1,1), (1,0,0), (0,0,α₃)} per iteration" begin
            for cap in 1:8
                alg = Lirion.aladin(; coupling = :state, sectioning = :tie_line,
                    max_iter = cap, rho = 1000.0, mu = 1000.0, adaptive = false,
                    rho_growth = 1.0, mu_growth = 1.0, globalization = :alg3)
                ws, iters = _paper_reg_solve(alg; model = Lirion.ac())
                iters < cap && break   # converged before the cap — trajectory done
                a = (ws.alpha1, ws.alpha2, ws.alpha3)
                @test a == (1.0, 1.0, 1.0) || a == (1.0, 0.0, 0.0) ||
                      (a[1] == 0.0 && a[2] == 0.0 && 0.0 <= a[3] <= 1.0)
            end
        end

        # ---- 4. globalization availability across models ------------
        # :alg3 and :full_step are now implemented for AC, DC and SOCWR
        # (aladin/{ac,dc,socwr}/line_search.jl); only unsupported base
        # models are barred.
        @testset "globalization :alg3/:full_step available for AC/DC/SOCWR" begin
            for bm in (Lirion.AC, Lirion.DC, Lirion.SOCWR),
                g in (:alg3, :full_step, :armijo_l2, :alg3_original)
                cfg = Lirion.ALADINConfig(
                    bm, :state, Lirion.TieLine(), 10; globalization = g)
                @test cfg isa Lirion.ALADINConfig
                @test cfg.globalization === g
            end
        end

        # ---- 5. :alg3_original (HFD Alg. 3 as printed) — AC ----------
        # Registers for AC (fixed or auto λ̄), is barred where no
        # `_alg3_original_step!` / `_dual_slope` exists, runs on case3 and
        # only ever produces the HFD step patterns; branch (c) outcomes are
        # logged in hist.meta[:alg3_original] with a status from the Wolfe
        # search and α₃ = 1 whenever the full dual step was accepted.
        @testset ":alg3_original registers (AC) with λ̄ policy" begin
            cfg = Lirion.ALADINConfig(
                Lirion.AC, :state, Lirion.TieLine(), 10;
                globalization = :alg3_original)
            @test cfg.globalization === :alg3_original
            @test cfg.alg3_lam_bar == 0.0
            cfg = Lirion.ALADINConfig(
                Lirion.AC, :hybrid, Lirion.SplitBus(), 10;
                globalization = :alg3_original, alg3_lam_bar = 250.0)
            @test cfg.alg3_lam_bar == 250.0
            @test_throws ArgumentError Lirion.ALADINConfig(
                Lirion.AC, :state, Lirion.TieLine(), 10; alg3_lam_bar = -1.0)
            # :alg3_original now available for AC/DC/SOCWR (construction tested
            # above); other modes never see cfg.alg3_lam_bar (bit-identity guard)
            for g in (:alg3, :full_step, :armijo_l2)
                cfg = Lirion.ALADINConfig(
                    Lirion.AC, :state, Lirion.TieLine(), 10;
                    globalization = g, alg3_lam_bar = 7.0)
                @test Lirion._alg3_original_lam_bar(cfg) == 0.0
            end
        end

        @testset ":alg3_original runs (AC) — HFD step patterns + (c) log" begin
            c_statuses = (:full_step, :wolfe, :armijo_only, :ascent_only, :no_ascent)
            for cap in 1:8
                alg = Lirion.aladin(; coupling = :state, sectioning = :tie_line,
                    max_iter = cap, rho = 1000.0, mu = 1000.0, adaptive = false,
                    rho_growth = 1.0, mu_growth = 1.0,
                    globalization = :alg3_original)
                ws, iters = _paper_reg_solve(alg; model = Lirion.ac())
                iters < cap && break   # converged before the cap
                a = (ws.alpha1, ws.alpha2, ws.alpha3)
                @test a == (1.0, 1.0, 1.0) || a == (1.0, 0.0, 0.0) ||
                      (a[1] == 0.0 && a[2] == 0.0 && 0.0 <= a[3] <= 1.0)
                log = ws.hist.meta[:alg3_original]
                @test length(log) == cap
                branch, a3, lam_bar, status_c = log[end]
                @test branch in (:a, :b, :c)
                @test lam_bar >= 1.0
                @test a3 == ws.alpha3
                if branch === :c
                    @test status_c in c_statuses
                    @test ws.alg3_c_status === status_c
                    status_c === :full_step && @test a3 == 1.0
                    status_c === :no_ascent && @test a3 == 0.0
                else
                    @test status_c === :none
                end
            end
            # λ̄ is monotone non-decreasing along the run (Theorem-2 policy)
            alg = Lirion.aladin(; coupling = :state, sectioning = :tie_line,
                max_iter = 8, rho = 1000.0, mu = 1000.0, adaptive = false,
                rho_growth = 1.0, mu_growth = 1.0, globalization = :alg3_original)
            ws, _ = _paper_reg_solve(alg; model = Lirion.ac())
            lams = [e[3] for e in ws.hist.meta[:alg3_original]]
            @test all(diff(lams) .>= 0.0)
            # fixed λ̄ is honoured verbatim
            alg = Lirion.aladin(; coupling = :state, sectioning = :tie_line,
                max_iter = 3, rho = 1000.0, mu = 1000.0, adaptive = false,
                rho_growth = 1.0, mu_growth = 1.0, globalization = :alg3_original,
                alg3_lam_bar = 123.0)
            ws, _ = _paper_reg_solve(alg; model = Lirion.ac())
            @test ws.lam_bar_fixed
            @test all(e[3] == 123.0 for e in ws.hist.meta[:alg3_original])
            # the κ̄ majoriser of :alg3_original walks the variable-bound
            # multipliers too (Φ penalises box violations); :alg3 keeps the
            # constraint-rows-only value
            m0 = Lirion._ac_max_local_dual(ws)
            m1 = Lirion._ac_max_local_dual(ws; include_bounds = true)
            @test m1 >= m0
            @test m0 == Lirion._ac_max_local_dual(ws; include_bounds = false)
            # terminal status + run manifest present for AC as well
            @test ws.hist.meta[:terminal_status] in (:converged, :iteration_limit, :stalled)
            rm = Lirion.aladin_run_manifest(ws)          # cfg optional: read from ws
            @test rm["GLOBALIZATION"] == "alg3_original"
            @test rm["LAM_BAR_FIXED"] === true && rm["LAM_BAR_APPLIED"] == 123.0
            @test rm["ALG3_LAM_BAR_CFG"] == 123.0
        end

        # ---- 6. fixed-λ̄ violation counter + manifest fields (AC) -------
        # The AC fixed-λ̄ path counts + warns when |λ|∞ exceeds the user-fixed
        # λ̄; λ̄ itself is never re-scaled. The configured λ̄ is recorded next
        # to the applied one, and the effective Σ is in the run manifest.
        @testset "fixed λ̄ violation counter (AC)" begin
            case3 = joinpath(@__DIR__, "pglib_opf_case3_lmbd.m")
            lam_fixed = 1.0e12                       # never exceeded by a real λ
            alg = Lirion.aladin(; coupling = :state, sectioning = :tie_line,
                max_iter = 3, rho = 1000.0, mu = 1000.0, adaptive = false,
                rho_growth = 1.0, mu_growth = 1.0, globalization = :alg3_original,
                alg3_lam_bar = lam_fixed)
            inputs = Lirion._build_inputs(case3, alg; model = Lirion.ac(),
                optimizer = nothing, outputs = false)
            @test inputs.model_type isa Lirion.ALADINConfig
            _, ws, _ = redirect_stdout(devnull) do
                Lirion.solve_distributed(inputs)
            end
            @test ws.lam_bar_fixed && ws.lam_bar == lam_fixed
            @test get(ws.hist.meta, :lam_bar_violations, 0) == 0
            rm = Lirion.aladin_run_manifest(ws, inputs.model_type)
            @test rm["LAM_BAR_VIOLATIONS"] == 0
            @test rm["ALG3_LAM_BAR_CFG"] == lam_fixed == rm["LAM_BAR_APPLIED"]
            # force |λ|∞ > λ̄ and re-run the λ̄ policy: counted, warned, λ̄ intact
            saved = [copy(s.lambda) for s in values(ws.vars)]
            ws.vars.theta.lambda .= 2.0 * lam_fixed
            @test_logs (:warn, r"exceeds the fixed L1 penalty weight") match_mode = :any begin
                @test Lirion._alg3_original_lam_bar!(ws) == lam_fixed
            end
            @test ws.hist.meta[:lam_bar_violations] == 1
            @test ws.lam_bar == lam_fixed                            # never re-scaled
            @test_logs (:warn, r"exceeds the fixed L1 penalty weight") match_mode = :any begin
                @test Lirion._alg3_original_lam_bar!(ws) == lam_fixed    # second hit
            end
            @test ws.hist.meta[:lam_bar_violations] == 2
            @test Lirion.aladin_run_manifest(ws, inputs.model_type)["LAM_BAR_VIOLATIONS"] == 2
            for (s, s0) in zip(values(ws.vars), saved)
                copyto!(s.lambda, s0)
            end
            @test Lirion._alg3_original_lam_bar!(ws) == lam_fixed
            @test ws.hist.meta[:lam_bar_violations] == 2             # no false hit
            # the automatic (non-fixed) policy never touches the counter
            alg0 = Lirion.aladin(; coupling = :state, sectioning = :tie_line,
                max_iter = 3, rho = 1000.0, mu = 1000.0, adaptive = false,
                rho_growth = 1.0, mu_growth = 1.0, globalization = :alg3_original)
            ws0, _ = _paper_reg_solve(alg0; model = Lirion.ac())
            @test !ws0.lam_bar_fixed
            @test !haskey(ws0.hist.meta, :lam_bar_violations)
            @test Lirion.aladin_run_manifest(ws0)["LAM_BAR_VIOLATIONS"] == 0
        end

        @testset "effective Σ + local-model state in the run manifest (AC)" begin
            for (coupling, mode) in ((:state, :primal), (:power, :pflow), (:hybrid, :primal))
                alg = Lirion.aladin(; coupling, sectioning = :tie_line,
                    max_iter = 2, rho = 1000.0, mu = 1000.0, adaptive = false,
                    rho_growth = 1.0, mu_growth = 1.0, proximal_mode = mode)
                ws, _ = _paper_reg_solve(alg; model = Lirion.ac())
                rm = Lirion.aladin_run_manifest(ws)
                sig = rm["SIGMA_EFFECTIVE"]
                @test sig isa Dict{String, Float64}
                expected = [String(f) for f in fieldnames(typeof(ws))
                            if startswith(String(f), "sigma_")]
                @test sort(collect(keys(sig))) == sort(expected)
                for (k, v) in sig
                    @test v == getproperty(ws, Symbol(k))
                end
                # AC never wires cfg.sigma_*: the recorded values are the
                # workspace defaults actually applied (that is the point)
                @test sig["sigma_theta"] == 10.0 && sig["sigma_vm"] == 5.0
                @test sig["sigma_pg"] == 1.0 && sig["sigma_qg"] == 1.0
                @test rm["PROXIMAL_MODE"] == String(mode)
                @test rm["LOCAL_MODELS"] in ("accepted", "restored")
                @test ws.hist.meta[:local_models] in (:accepted, :restored)
            end
        end

        # ---- 7. config guards -----------------------------------------
        # Sectioning × coupling combinations whose exchange the ALADIN
        # workspaces do not couple must be rejected (DummyBus + :state,
        # SplitBus + :state / :power). Unwired parameters must not be
        # silently accepted either.
        @testset "sectioning × coupling guard" begin
            for bm in (Lirion.AC, Lirion.DC, Lirion.SOCWR)
                for c in (:state, :power)
                    @test_throws ArgumentError Lirion.ALADINConfig(bm, c, Lirion.SplitBus(), 10)
                end
                @test Lirion.ALADINConfig(bm, :hybrid, Lirion.SplitBus(), 10) isa Lirion.ALADINConfig
            end
            for bm in (Lirion.AC, Lirion.DC)
                @test_throws ArgumentError Lirion.ALADINConfig(bm, :state, Lirion.DummyBus(), 10)
                for c in (:power, :hybrid)
                    @test Lirion.ALADINConfig(bm, c, Lirion.DummyBus(), 10) isa Lirion.ALADINConfig
                end
            end
        end

        @testset "unwired parameters fail loudly" begin
            @test_throws ArgumentError Lirion.ALADINConfig(
                Lirion.DC, :state, Lirion.TieLine(), 10; active_set_tol = 1e-5)
            for bm in (Lirion.AC, Lirion.SOCWR), kw in ((; sigma_pg = 2.0),
                (; sigma_other = 2.0), (; sigma_theta_shared = 5.0))
                @test_throws ArgumentError Lirion.ALADINConfig(
                    bm, :state, Lirion.TieLine(), 10; kw...)
            end
            # DC keeps them live
            cfg = Lirion.ALADINConfig(Lirion.DC, :state, Lirion.TieLine(), 10;
                sigma_pg = 2.0, sigma_other = 3.0, sigma_theta_shared = 5.0)
            @test (cfg.sigma_pg, cfg.sigma_other, cfg.sigma_theta_shared) == (2.0, 3.0, 5.0)
        end

        # `_ac_prox` in AC :pflow must carry the whole proximal term the zones
        # minimise (quad_static part + nonlinear flow proximal). Reference =
        # the local JuMP objectives themselves:
        #   Σᵢ [Pᵢ(y) − Pᵢ(x)],  Pᵢ := objectiveᵢ − fᵢ − (λ-terms), Pᵢ(x) = min.
        # With λ = 0 (first iteration) the λ-terms vanish, so
        #   Σᵢ [(objᵢ − fᵢ)(y) − (objᵢ − fᵢ)(x)] == _ac_prox(ws).
        @testset "_ac_prox carries the :pflow flow proximal" begin
            alg = Lirion.aladin(; coupling = :power, sectioning = :tie_line,
                max_iter = 5, rho = 1000.0, mu = 1000.0, proximal_mode = :pflow,
                sigma_pflow = 1.0, sigma_qflow = 1.0)
            case3 = joinpath(@__DIR__, "pglib_opf_case3_lmbd.m")
            inputs = Lirion._build_inputs(case3, alg; model = Lirion.ac(),
                optimizer = nothing, outputs = false)
            cfg = inputs.model_type
            ws = Lirion.build_workspace(inputs.system, cfg,
                Lirion.resolve_optimizer(cfg, inputs.optimizer),
                Lirion.ALADIN{Lirion.Power, Lirion.AC}())
            redirect_stdout(devnull) do
                ws.iter = 2                      # keep the seeded x_prev (≠ y)
                Lirion.local_solve!(ws)
            end
            @test all(all(iszero, st.lambda) for st in values(ws.vars))
            ref = 0.0
            for (i, part) in enumerate(ws.partitions.parts)
                z = ws.zones[i]
                pen(pt) = Lirion.JuMP.value(v -> pt[z.var_to_col[v]],
                              Lirion.JuMP.objective_function(part.local_model)) -
                          Lirion.JuMP.value(v -> pt[z.var_to_col[v]], part.obj_base)
                ref += pen(z.y) - pen(z.x_prev)
            end
            @test Lirion._ac_flow_prox_value(ws) > 0.0
            @test Lirion._ac_prox(ws) ≈ ref rtol = 1e-8
        end
    end
end
