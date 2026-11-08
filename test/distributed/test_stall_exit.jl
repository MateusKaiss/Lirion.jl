# ============================================================
# ===========  Frozen early exit status  =====================
# ============================================================
# A run frozen for `_STALL_EXIT_STEPS[]` consecutive iterations (α = (0,0,0),
# ρ, μ unchanged) exits before max_iter with terminal status :stalled, and
# the recorded terminal status is authoritative in LirionResult.status.
#   (1) a frozen :alg3_original run exits early with :stalled (never :converged);
#   (2) recorded :stalled / :iteration_limit win over the iteration-count
#       heuristic, which is only the fallback when nothing was recorded.
# The freeze in (1) is forced by a temporary `_line_search_step!` method on
# the DC :state workspace type, removed with Base.delete_method in `finally`.
# ============================================================

function _stall_exit_inputs(; max_iter::Int, eps::Float64)
    case3 = joinpath(@__DIR__, "pglib_opf_case3_lmbd.m")
    alg = Lirion.aladin(; coupling = :state, sectioning = :tie_line,
        max_iter = max_iter, rho = 1000.0, mu = 1000.0, adaptive = false,
        eps_abs = eps, eps_rel = eps, globalization = :alg3_original)
    return Lirion._build_inputs(case3, alg; model = Lirion.dc(),
        optimizer = nothing, outputs = false)
end

function run_stall_exit_tests()
    @testset "frozen early exit reports :stalled (never :converged)" begin
        # ---- (1) genuine frozen run, cut by the early exit ------------------
        @testset "early exit → hist :stalled → LirionResult.status :stalled" begin
            max_iter = 40
            inputs = _stall_exit_inputs(; max_iter, eps = 1e-12)   # cannot certify
            old_steps = Lirion._STALL_EXIT_STEPS[]
            Lirion._STALL_EXIT_STEPS[] = 3
            forced = nothing
            try
                # force α = (0,0,0) on every globalization call (most specific method)
                Core.eval(Lirion, quote
                    function _line_search_step!(ws::DCStateWorkspaceALADIN)
                        ws.alpha1 = 0.0; ws.alpha2 = 0.0; ws.alpha3 = 0.0
                        return nothing
                    end
                end)
                forced = which(Lirion._line_search_step!, Tuple{Lirion.DCStateWorkspaceALADIN})
                parts, ws, it = redirect_stdout(devnull) do
                    Base.invokelatest(Lirion.solve_distributed, inputs)
                end
                @test it < max_iter                                   # cut early
                @test it == Lirion._STALL_EXIT_STEPS[]                # frozen from iteration 1
                @test get(ws.hist.meta, :early_stall_exit, 0) == it
                @test ws.hist.meta[:terminal_status] === :stalled
                @test ws.hist.meta[:zero_step_iters] == it
                @test ws.alpha1 == 0.0 && ws.alpha2 == 0.0 && ws.alpha3 == 0.0
                r = Base.invokelatest(Lirion.build_lirion_result_distributed,
                    parts, ws, inputs.model_type, it, time_ns(), inputs.system)
                @test r.status === :stalled
                @test r.status !== :converged
                runm = Lirion.aladin_run_manifest(ws, inputs.model_type)
                @test runm["TERMINAL_STATUS"] == "stalled"
                @test runm["EARLY_STALL_EXIT"] == it
                @test runm["ZERO_STEP_ITERS"] == it
            finally
                Lirion._STALL_EXIT_STEPS[] = old_steps
                forced === nothing || Base.delete_method(forced)
            end
            # the forced method is gone: the union method is dispatched again
            @test which(Lirion._line_search_step!, Tuple{Lirion.DCStateWorkspaceALADIN}).sig !==
                  Tuple{typeof(Lirion._line_search_step!), Lirion.DCStateWorkspaceALADIN}
        end

        # ---- (2) the status mapping of build_lirion_result_distributed --------
        @testset "recorded terminal status is authoritative" begin
            max_iter = 3
            inputs = _stall_exit_inputs(; max_iter, eps = 1e-12)
            parts, ws, it = redirect_stdout(devnull) do
                Lirion.solve_distributed(inputs)
            end
            @test it == max_iter
            meta = ws.hist.meta
            saved = get(meta, :terminal_status, nothing)
            try
                res(itd) = Lirion.build_lirion_result_distributed(
                    parts, ws, inputs.model_type, itd, time_ns(), inputs.system).status
                # recorded status wins regardless of the iteration count
                meta[:terminal_status] = :stalled
                @test res(1) === :stalled
                @test res(max_iter) === :stalled
                meta[:terminal_status] = :iteration_limit
                @test res(1) === :iteration_limit
                meta[:terminal_status] = :converged
                @test res(max_iter) === :converged
                # no recorded status: the iteration-count heuristic is the fallback
                delete!(meta, :terminal_status)
                @test res(1) === :converged
                @test res(max_iter) === :iteration_limit
            finally
                saved === nothing ? delete!(meta, :terminal_status) :
                                    (meta[:terminal_status] = saved)
            end
        end
    end
end
