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

        # Skip adaptive-rho at iter 1: consensus_prev is initialized to
        # zeros, so the dual residual s = ρ·‖z − 0‖ is artificially large
        # compared to the primal residual r = ‖x − z‖. Firing the
        # residual-balancing rule on that spurious s triggers a cascade
        # of shrink steps that decays ρ by orders of magnitude over the
        # remainder of the run (observed ρ: 1e3 → 1 in 500 iters on
        # case14–case57). From iter 2 onwards consensus_prev holds a real
        # previous consensus vector, so the rule measures genuine progress.
        if iter >= 2
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
