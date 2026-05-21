# ============================================================
# Core loop
# ============================================================

@inline stop_if_converged(ws::AbstractALADINWorkspace) =
    (ws.r_primal < ws.eps_primal) &&
    (ws.r_step_primal < ws.eps_step_primal) &&
    (ws.r_dual < ws.eps_dual)

@inline _line_search_step!(::AbstractALADINWorkspace) = nothing

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

# ============================================================
# Entry point
# ============================================================

function run_algorithm(
    sys::DistributedSystem,
    config::ALADINConfig;
    optimizer_factory = nothing,
)
    t0_solver = time_ns()
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
            _finalize(ws.partitions, t0_solver, ws.iter);
            rho_lines = rho_lines,
        )
    end

    return ws.partitions, ws, ws.iter
end
