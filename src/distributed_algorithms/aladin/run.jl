# =======================================================================
# ALADIN — Top-level runner
# =======================================================================

@inline function aladin_tag(cfg::ALADINConfig)
    V = _variant_type(Val(cfg.coupling))
    M = cfg.base_model
    return ALADIN{V, M}()
end

algorithm_tag(::ALADINAlgorithm) = :aladin

# ---------------------------------------------------------------
# FINALIZE
# ---------------------------------------------------------------

@inline function _finalize_aladin(partitions, t0_ns::Integer, iter_done::Integer)
    total_cost = _total_generation_cost(partitions)
    elapsed_s = _elapsed_s(t0_ns)
    return (;
        total_cost,
        partitions,
        iter_done = Int(iter_done),
        elapsed_s,
    )
end

function run_aladin!(ws::AbstractALADINWorkspace, params::ALADINParams)
    tag = _aladin_print_tag(ws)

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

function run_aladin(
    sys::DistributedSystem,
    config::ALADINConfig,
    tag::ALADIN;
    optimizer_factory = nothing,
)::LirionRunResult
    t0 = time_ns()

    params = build_aladin_params(config)
    ws = build_aladin_workspace(sys, config, optimizer_factory, tag)
    initial_rhos = _get_rhos(ws)
    run_aladin!(ws, params)
    result = _finalize_aladin(ws.partitions, t0, ws.iter)

    tag_str = _aladin_print_tag(ws)
    final_rhos = _get_rhos(ws)

    rho_lines = [
        "$(k) inicial: $(initial_rhos[k]) | $(k) final: $(final_rhos[k])"
        for k in keys(final_rhos)
    ]
    push!(rho_lines, "mu: $(ws.mu[])")

    _FAST_USE_LOGS && _print_summary!(
        tag_str,
        "ALADIN",
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
