# ============================================================
# ==========  ALADIN coordinator-QP dual-sign guard  =========
# ============================================================

using JuMP

function run_aladin_qp_dual_sign_tests(; atol = 1e-5)
    @testset "ALADIN coordinator-QP dual sign convention" begin
        h = 2.0      # diagonal Hessian (zone cost curvature)
        mu = 1.0     # ALADIN consensus-slack penalty
        r0 = 0.3     # current consensus residual y1 - y2
        lambda = 0.5 # current coupling multiplier

        # Solve with the SAME optimizer/attributes the coordinator QP uses.
        qp = Model(Lirion._DEFAULT_QP_FACTORY)
        set_silent(qp)

        @variable(qp, dy1)
        @variable(qp, dy2)
        @variable(qp, s)

        @objective(qp, Min,
            0.5 * h * dy1^2 + 0.5 * h * dy2^2 + lambda * s + 0.5 * mu * s^2,
        )

        # _bus_coupling_expr! form: (y1 - y2) + (dy1 - dy2) == s
        cons = @constraint(qp, r0 + dy1 - dy2 == s)

        optimize!(qp)
        @test termination_status(qp) in (
            Lirion.MOI.OPTIMAL,
            Lirion.MOI.LOCALLY_SOLVED,
            Lirion.MOI.ALMOST_OPTIMAL,
        )

        # This is the exact recovery used in every solve_coordinator_qp!.
        nu = -dual(cons)

        nu_true = (mu * r0 + lambda) / (1 + 2 * mu / h)  # = 0.4

        # Value AND sign must match: a flipped MOI/OSQP convention (or dropping
        # the negation) yields -nu_true and fails here.
        @test isapprox(nu, nu_true; atol = atol)
        @test nu > 0                       # sign guard for this instance

        # Primal sanity (ties the multiplier to the QP it came from).
        @test isapprox(value(dy1), -nu_true / h; atol = atol)
        @test isapprox(value(dy2), nu_true / h; atol = atol)
        @test isapprox(value(s), (nu_true - lambda) / mu; atol = atol)
    end
    return nothing
end
