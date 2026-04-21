# ------------------------------------------------------------
# CORE LOOP
# ------------------------------------------------------------
function run_algorithm!(ws, params::ADMMParams)
    tag = _print_tag(ws)

    for iter in 1:params.max_iter
        ws.iter = iter

        _print_iter_header!(tag, iter)

        local_solve!(ws)
        save_prev!(ws)
        consensus!(ws)
        dual_update!(ws)
        compute_residuals!(ws)
        compute_tolerances!(ws, params.eps_abs, params.eps_rel)

        _record_iter!(ws)
        _print_admm_residuals!(ws)

        if stop_if_converged(ws)
            _print_converged!()
            break
        end

        # Warm-up period for adaptive ρ. Two pathologies motivate this:
        #
        #  (a) At iter 1 `consensus_prev` is zeros, so the dual residual
        #      s = ρ·‖z − 0‖ is artificially huge and would cascade a
        #      chain of shrink steps (observed ρ: 1e3 → 1 in 500 iters
        #      on case14–case57 before this guard was added).
        #  (b) In iters 2..~5 the subproblems are still moving their
        #      primals sharply; s ≈ ρ·‖Δz‖ is still an order of magnitude
        #      above its eventual value while r has barely changed. The
        #      rule fires "s ≫ r → shrink" repeatedly, again decaying ρ
        #      into a region too small to enforce state consensus
        #      (angles drift at tie-line interfaces).
        #
        # Waiting `adapt_warmup` iterations before the first adaptation
        # lets the dual residual settle to its steady-state scale so the
        # residual-balancing rule measures genuine imbalance rather than
        # startup transients. Default 5 matches empirical convergence
        # profiles on case14 / case30 / case57 DC-state.
        if iter > params.adapt_warmup
            rho_update!(ws, params)
            refresh_rho!(ws)
        end
    end

    return nothing
end

function run_algorithm(
    sys::DistributedSystem,
    config::ADMMConfig;
    optimizer_factory = nothing,
)::LirionResult
    t0 = time_ns()
    tag = ADMM{_variant_type(Val(config.coupling)), config.base_model}()

    params = build_params(config)
    ws = build_workspace(sys, config, optimizer_factory, tag)
    initial_rhos = _get_rhos(ws)
    run_algorithm!(ws, params)

    tag_str = _print_tag(ws)
    final_rhos = _get_rhos(ws)
    rho_lines = [
        "$(k) inicial: $(initial_rhos[k]) | $(k) final: $(final_rhos[k])"
        for k in keys(final_rhos)
    ]
    if _FAST_USE_LOGS
        _print_summary!(
            tag_str, "ADMM",
            _finalize(ws.partitions, t0, ws.iter);
            rho_lines = rho_lines,
        )
    end

    return build_lirion_result_distributed(
        ws.partitions, ws, config, ws.iter, t0, sys,
    )
end
