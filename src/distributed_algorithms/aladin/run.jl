algorithm_tag(::ALADINAlgorithm) = :aladin

function run_algorithm!(ws::AbstractALADINWorkspace, params::ALADINParams)
    tag = _print_tag(ws)

    for iter in 1:params.max_iter
        ws.iter = iter
        _print_iter_header!(tag, iter)

        local_solve!(ws)
        solve_coordinator_qp!(ws)

        compute_residuals!(ws)
        compute_tolerances!(ws, params.eps_abs, params.eps_rel)

        _record_iter!(ws)
        _print_aladin_residuals!(ws)

        if stop_if_converged(ws)
            _print_converged!()
            break
        end

        _line_search_step!(ws)
        primal_update!(ws)
        dual_update!(ws)

        adapt_rho_mu!(ws, params)
    end
    return nothing
end

# ---------------------------------------------------------------
# ENTRY POINT
# ---------------------------------------------------------------

function run_algorithm(
    sys::DistributedSystem,
    config::ALADINConfig;
    optimizer_factory = nothing,
)::LirionResult
    t0 = time_ns()
    tag = ALADIN{_variant_type(Val(config.coupling)), config.base_model}()

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
    push!(rho_lines, "mu: $(ws.mu[])")

    if _FAST_USE_LOGS
        _print_summary!(
            tag_str, "ALADIN",
            _finalize(ws.partitions, t0, ws.iter);
            rho_lines = rho_lines,
        )
    end

    return build_lirion_result_distributed(
        ws.partitions, ws, config, ws.iter, t0, sys,
    )
end
