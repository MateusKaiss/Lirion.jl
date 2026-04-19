# ============================================================
# Renewables — invariants for curtailable generators.
#
# A generator is tagged renewable when column 26 of `mpc.gen` is `'R'`.
# The model adds a `curtailment[g, t]` variable and the envelope
# equality  p_gen[g, t] + curtailment[g, t] == pmax_ts[g, t].
# When `mpc.gen_timeseries` is present, `pmax_ts[g, t]` overrides the
# static `pmax[g]` nameplate per period.
#
# Tests here cover:
#   1. Static availability (no gen_timeseries, curt_cost=0) — curtailment
#      variable exists, envelope holds, objective matches the un-tagged
#      baseline (since curtailment is free).
#   2. Time-varying availability — pmax_ts matches the input profile;
#      envelope holds at every period.
#   3. Zero availability — at a period where pmax_ts=0, both p_gen and
#      curtailment are forced to 0.
#   4. Curtailment cost term — when curt_cost > 0 and supply > demand,
#      the curtailment cost appears in the objective and equals
#      sbase · curt_cost · Σ curtailment.
#   5. Distributed parity — ADMM DC on a tagged single-period renewable
#      case matches centralized DC.
# ============================================================

import Lirion
import JuMP

const _RENEW_DIR = joinpath(@__DIR__, "renewables")

# Helpers ---------------------------------------------------------------

@inline _case(name::AbstractString) = joinpath(_RENEW_DIR, name)

function _solve(file, model; optimizer = Ipopt.Optimizer)
    return Lirion.solve(
        file;
        model = model, algorithm = Lirion.centralized(),
        optimizer = optimizer, outputs = false, plots = false,
    )
end

# Invariant 1 — Static availability (curtailment var + envelope) -------

function _check_envelope_holds(result::Lirion.LirionResult; atol = 1e-5)
    m = result.model
    T = get(m.ext, :time_periods, 1)
    @test haskey(m, :curtailment)

    gen = result.system.gen
    curt = m[:curtailment]
    pgen = m[:p_generation]
    for g in gen.curtailable, t in 1:T
        lhs = JuMP.value(pgen[g, t]) + JuMP.value(curt[g, t])
        rhs = gen.pmax_ts[g, t]
        @test isapprox(lhs, rhs; atol = atol)
    end
    return nothing
end

# Invariant 4 — Curtailment cost enters the objective -------------------

function _check_curtailment_cost_term(result::Lirion.LirionResult; atol = 1e-4)
    m = result.model
    gen = result.system.gen
    @test haskey(m, :curtailment_cost_expr)
    T = get(m.ext, :time_periods, 1)
    curt = m[:curtailment]
    sbase = result.system.baseMVA
    expected = sum(
        gen.curtailment_cost[g] * sbase * JuMP.value(curt[g, t])
        for g in gen.curtailable, t in 1:T
    )
    @test isapprox(JuMP.value(m[:curtailment_cost_expr]), expected; atol = atol)
    return nothing
end

# Public runner ---------------------------------------------------------

function run_renewables_tests(; optimizer = Ipopt.Optimizer)
    @testset "Renewables" begin
        @testset "1. Static availability (curt_cost=0)" begin
            # The untagged baseline: the same topology/data without any
            # plant_type="R" tag — objective is the pure economic dispatch
            # optimum. With curt_cost=0 and availability==nameplate, the
            # renewable-tagged run must reach the same objective.
            base = _solve(
                joinpath(@__DIR__, "matpower_cases", "pglib_opf_case3_lmbd.m"),
                Lirion.dc(); optimizer = optimizer,
            )
            tagged = _solve(_case("case3_renew_static.m"), Lirion.dc();
                optimizer = optimizer)

            @test tagged.status in (:optimal, :local_optimal)
            _check_envelope_holds(tagged)
            @test isapprox(tagged.objective, base.objective; rtol = 1e-5)
        end

        @testset "2. Time-varying availability" begin
            r = _solve(_case("case3_renew_ts.m"), Lirion.dc();
                optimizer = optimizer)
            @test r.status in (:optimal, :local_optimal)
            # Availability profile in pu (baseMVA=100): 50/150/200 MW → 0.5/1.5/2.0.
            expected_pmax_ts = [0.5, 1.5, 2.0]
            for t in 1:3
                @test isapprox(r.system.gen.pmax_ts[2, t], expected_pmax_ts[t];
                    atol = 1e-10)
            end
            _check_envelope_holds(r)
        end

        @testset "3. Zero availability forces zero dispatch" begin
            r = _solve(_case("case3_renew_zero.m"), Lirion.dc();
                optimizer = optimizer)
            @test r.status in (:optimal, :local_optimal)
            m = r.model
            # At t=1, pmax_ts[2,1] == 0 → both p_gen and curtailment must
            # be 0 (the envelope equality is 0 = 0).
            @test isapprox(JuMP.value(m[:p_generation][2, 1]), 0.0; atol = 1e-6)
            @test isapprox(JuMP.value(m[:curtailment][2, 1]), 0.0; atol = 1e-6)
            _check_envelope_holds(r)
        end

        @testset "4. Curtailment cost term in objective" begin
            r = _solve(_case("case3_renew_curtcost.m"), Lirion.dc();
                optimizer = optimizer)
            @test r.status in (:optimal, :local_optimal)
            _check_envelope_holds(r)
            _check_curtailment_cost_term(r)
            # Since renewable capacity (1000 MW) dwarfs demand (315 MW),
            # some curtailment is unavoidable; with curt_cost=50 /MW the
            # objective must EXCEED the untagged baseline (≈ 5696) by at
            # least sbase · curt_cost · min_forced_curtailment.
            min_curt_mw = 1000.0 - 315.0  # lower bound on forced curtailment
            @test r.objective >= 5696.0 + 50.0 * min_curt_mw - 1.0
        end

        @testset "5. Extra gen_timeseries columns (pmin/qmax/qmin)" begin
            # case3_gen_ts_bounds: gen 1 must-run 300 MW at t=1, gen 2
            # qmax tightened to 50 MVAr at t=2.
            r = _solve(_case("case3_gen_ts_bounds.m"), Lirion.ac();
                optimizer = optimizer)
            @test r.status in (:optimal, :local_optimal)
            g = r.system.gen
            m = r.model

            # gen 1: static pmin = 0, pmin_ts[1, 1] = 3.0 (pu), pmin_ts[1, 2] = 0.
            @test isapprox(g.pmin_ts[1, 1], 3.0; atol = 1e-10)
            @test isapprox(g.pmin_ts[1, 2], 0.0; atol = 1e-10)
            @test JuMP.value(m[:p_generation][1, 1]) >= 3.0 - 1e-6

            # gen 2: static qmax = 10.0 (pu), qmax_ts[2, 2] = 0.5.
            @test isapprox(g.qmax_ts[2, 1], 10.0; atol = 1e-10)
            @test isapprox(g.qmax_ts[2, 2], 0.5; atol = 1e-10)
            @test JuMP.value(m[:q_generation][2, 2]) <= 0.5 + 1e-6

            # NaN sentinel: other entries must fall back to static.
            @test isapprox(g.pmax_ts[1, 1], g.pmax[1]; atol = 1e-10)
            @test isapprox(g.qmin_ts[2, 1], g.qmin[2]; atol = 1e-10)
        end

        @testset "6. Distributed ADMM DC matches centralized" begin
            # Distributed ADMM is single-period only today, so use the
            # static renewable case. Both runs must land on the same
            # objective (within the usual ADMM tolerance).
            file = _case("case3_renew_static.m")
            cent = _solve(file, Lirion.dc(); optimizer = optimizer)

            admm = Lirion.solve(
                file;
                model = Lirion.dc(),
                algorithm = Lirion.admm(;
                    coupling = :hybrid,
                    sectioning = :tie_line,
                    max_iter = 500,
                    rho = 1e3,
                ),
                optimizer = optimizer,
                outputs = false, plots = false,
            )
            @test isapprox(admm.objective, cent.objective; rtol = 5e-3)
        end
    end
    return nothing
end
