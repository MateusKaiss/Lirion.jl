# ------------------------------------------------------------
# ADMM type-level tag used by run_admm dispatch
# ------------------------------------------------------------
@inline function admm_tag(cfg::DistributedOPFModel)
    V = _variant_type(Val(cfg.coupling))
    M = cfg.base_model
    return ADMM{V, M}()
end

# ------------------------------------------------------------
# FINALIZE
# ------------------------------------------------------------
@inline function _finalize_admm(partitions, t0_ns::Integer, iter_done::Integer)
    total_cost = _admm_totals(partitions)
    elapsed_s = _elapsed_s(t0_ns)
    return (;
        total_cost,
        partitions,
        iter_done = Int(iter_done),
        elapsed_s,
    )
end

# ------------------------------------------------------------
# CORE LOOP
# ------------------------------------------------------------
function run_admm!(ws, params::ADMMParams)
    tag = _admm_print_tag(ws)

    for iter in 1:params.max_iter
        ws.iter = iter

        _print_admm_iter_header!(tag, iter)

        local_solve!(ws)
        save_prev!(ws)
        consensus!(ws)
        dual_update!(ws)
        compute_residuals!(ws)
        compute_tolerances!(ws, params.eps_abs, params.eps_rel)

        _record_iter!(ws)
        _print_admm_residuals!(ws)

        if stop_if_converged(ws)
            _print_admm_converged!()
            break
        end

        rho_update!(ws, params)
        refresh_rho!(ws)
    end

    return nothing
end

function run_admm(
    sys::DistributedSystem,
    config::DistributedOPFModel,
    tag::ADMM;
    optimizer_factory = nothing,
)::LirionRunResult
    t0 = time_ns()

    params = build_admm_params(config)
    ws = build_admm_workspace(sys, config, optimizer_factory, tag)
    initial_rhos = _get_rhos(ws)
    run_admm!(ws, params)
    result = _finalize_admm(ws.partitions, t0, ws.iter)

    tag_str = _admm_print_tag(ws)
    final_rhos = _get_rhos(ws)

    rho_lines = [
        "$(k) inicial: $(initial_rhos[k]) | $(k) final: $(final_rhos[k])"
        for k in keys(final_rhos)
    ]
    _FAST_USE_LOGS && _print_admm_summary!(
        tag_str,
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
