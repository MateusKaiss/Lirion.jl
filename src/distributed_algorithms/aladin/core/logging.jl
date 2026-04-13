# ------------------------------------------------------------
# PRINT TAG
# ------------------------------------------------------------

_aladin_print_tag(::DCStateWorkspaceALADIN) = "[ALADIN-state/DC]"
_aladin_print_tag(::DCPowerWorkspaceALADIN) = "[ALADIN-power/DC]"
_aladin_print_tag(::DCHybridWorkspaceALADIN) = "[ALADIN-hybrid/DC]"
_aladin_print_tag(::AbstractALADINWorkspace) = "[ALADIN]"

# ------------------------------------------------------------
# ITERATION PRINTS
# ------------------------------------------------------------

@inline function _print_aladin_iter_header!(tag::AbstractString, iter::Int)
    _FAST_USE_LOGS && println("\n$tag Iteration $iter")
    return nothing
end

@inline function _print_aladin_residuals!(ws::AbstractALADINWorkspace)
    _FAST_USE_LOGS && println(
        "‖r_primal‖ = $(ws.r_primal) | eps_pri  = $(ws.eps_primal)",
    )
    _FAST_USE_LOGS && println(
        "‖r_progress‖ = $(ws.r_progress) | eps_prog = $(ws.eps_progress)",
    )
    _FAST_USE_LOGS && println(
        "‖Δλ‖ = $(ws.lambda_step)  (diagnostic)",
    )
    return nothing
end

@inline function _print_aladin_converged!()
    _FAST_USE_LOGS && println("Stopping criterion met.")
    return nothing
end

# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------

function _print_aladin_summary!(
    tag::AbstractString,
    result;
    rho_lines = (),
)
    println("\n$tag Resumo:")
    println("→ Total generation cost = $(result.total_cost)")
    for line in rho_lines
        println(line)
    end
    println("Iterações realizadas: $(result.iter_done)")
    println("Tempo ALADIN: $(_pretty_time(result.elapsed_s))")
    return nothing
end
