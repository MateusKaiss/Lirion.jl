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
    total_cost = _aladin_totals(partitions)
    elapsed_s = _elapsed_s(t0_ns)
    return (;
        total_cost,
        partitions,
        iter_done = Int(iter_done),
        elapsed_s,
    )
end

# ---------------------------------------------------------------
# CORE LOOP
# ---------------------------------------------------------------

function run_aladin!(ws::AbstractALADINWorkspace, params::ALADINParams)
    return error("run_aladin! not yet implemented (Milestone 1c)")
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
    run_aladin!(ws, params)
    result = _finalize_aladin(ws.partitions, t0, ws.iter)

    return LirionRunResult(;
        result = result,
        model = nothing,
        system = sys,
        diagnostics = NoDiagnostics(),
    )
end
