# ============================================================
# Zone warm-start (OPT-IN, OFF by default)
# ============================================================

const _ZONE_WARMSTART_ENV = "LIRION_ZONE_WARMSTART"

@inline function _zone_warmstart_env_enabled()::Bool
    v = lowercase(strip(get(ENV, _ZONE_WARMSTART_ENV, "")))
    return v == "1" || v == "true" || v == "yes" || v == "on"
end

@inline zone_warmstart_enabled(config::ADMMConfig)::Bool =
    config.zone_warmstart || _zone_warmstart_env_enabled()

# Returns true when the warm start was applied, false on fallback to cold start.
function zone_warmstart!(ws::AbstractADMMWorkspace)
    tag = _print_tag(ws)
    parts = ws.partitions.parts
    n = length(parts)
    failures = Vector{Union{Nothing, String}}(nothing, n)
    t0 = time_ns()

    Threads.@threads for i in eachindex(parts)
        part = parts[i]
        model = part.local_model
        try
            # Base objective only: no (ρ/2)‖·‖² term, no linear dual term.
            @objective(model, Min, part.obj_base)
            JuMP.optimize!(model)
            # Same acceptance ruler as the loop (iter = 0 marks the pre-solve).
            _require_optimal!(model; iter = 0, zone_id = part.id)
        catch err
            failures[i] = sprint(showerror, err)
        end
    end

    failed = [parts[i].id for i in 1:n if failures[i] !== nothing]
    if !isempty(failed)
        detail = join(
            ("zone $(parts[i].id): $(failures[i])" for i in 1:n if failures[i] !== nothing),
            "\n",
        )
        @warn "$(tag) zone warm-start FAILED in $(length(failed))/$(n) zones — " *
              "falling back to the cold start (z = 0, y = 0).\n$(detail)"
        _FAST_USE_LOGS && println(
            "$(tag) zone warm-start: fallback to cold start (failed zones $(failed))",
        )
        ws.hist.meta[:zone_warmstart] = (
            applied = false, n_zones = n, failed_zones = failed,
            elapsed_s = _elapsed_s(t0),
        )
        return false
    end

    # z⁰ from the zone solutions; `theta` is gauge-dependent and stays cold.
    theta_cold = haskey(ws.vars, :theta) ? copy(ws.vars.theta.consensus) : nothing
    consensus!(ws)
    theta_cold === nothing || copyto!(ws.vars.theta.consensus, theta_cold)
    save_prev!(ws)   # z⁻¹ := z⁰

    _FAST_USE_LOGS && println(
        "$(tag) zone warm-start: applied ($(n) zones, $(_pretty_time(_elapsed_s(t0))))",
    )
    ws.hist.meta[:zone_warmstart] = (
        applied = true, n_zones = n, failed_zones = Int[],
        elapsed_s = _elapsed_s(t0), theta_kept_cold = theta_cold !== nothing,
    )
    return true
end
