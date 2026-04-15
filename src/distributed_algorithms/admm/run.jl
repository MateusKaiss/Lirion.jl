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

        rho_update!(ws, params)
        refresh_rho!(ws)
    end

    return nothing
end

function run_algorithm(
    sys::DistributedSystem,
    config::ADMMConfig;
    optimizer_factory = nothing,
)::LirionRunResult
    t0 = time_ns()
    tag = ADMM{_variant_type(Val(config.coupling)), config.base_model}()

    params = build_params(config)
    ws = build_workspace(sys, config, optimizer_factory, tag)
    initial_rhos = _get_rhos(ws)
    run_algorithm!(ws, params)
    result = _finalize(ws.partitions, t0, ws.iter)

    tag_str = _print_tag(ws)
    final_rhos = _get_rhos(ws)

    rho_lines = [
        "$(k) inicial: $(initial_rhos[k]) | $(k) final: $(final_rhos[k])"
        for k in keys(final_rhos)
    ]
    _FAST_USE_LOGS && _print_summary!(
        tag_str,
        "ADMM",
        result;
        rho_lines = rho_lines,
    )

    return LirionRunResult(;
        result = result,
        model = nothing,
        system = sys,
        diagnostics = ws.hist,
    )
end
