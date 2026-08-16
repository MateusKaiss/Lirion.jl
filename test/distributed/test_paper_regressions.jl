# ============================================================
# ==============  PAPER REGRESSION TESTS (BL-2)  =============
# ============================================================
# Locks the four measurement-layer behaviours the paper campaign
# depends on (audit round 2026-08-04):
#   1. derived stationarity ruler: eps_step_primal == n_prox · eps_abs
#   2. structural mirror: _extra_prox_support! covers exactly the columns
#      _add_extra_prox_grad! touches (AC :pflow flow-proximal)
#   3. :alg3 globalization runs and only ever produces the HFD step
#      patterns (a) α=(1,1,1), (b) α=(1,0,0), (c) α=(0,0,α₃∈[0,1])
#   4. B1 guard: globalization=:alg3 with DC/SOCWR throws (AC-only impl)
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
    @testset "Paper regressions (derived ruler / prox mirror / alg3 / B1)" begin
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

        # ---- 4. B1 guard: :alg3 is AC-only --------------------------
        @testset "B1: globalization=:alg3 throws for DC/SOCWR" begin
            @test_throws ArgumentError Lirion.ALADINConfig(
                Lirion.DC, :state, Lirion.TieLine(), 10; globalization = :alg3)
            @test_throws ArgumentError Lirion.ALADINConfig(
                Lirion.SOCWR, :state, Lirion.TieLine(), 10; globalization = :alg3)
            cfg = Lirion.ALADINConfig(
                Lirion.AC, :state, Lirion.TieLine(), 10; globalization = :alg3)
            @test cfg isa Lirion.ALADINConfig
        end
    end
end
