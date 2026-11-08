# ============================================================
# Shared
# ============================================================

@inline function _print_iter_header!(tag::AbstractString, iter::Int)
    _FAST_USE_LOGS && println("\n$tag Iteration $iter")
    return nothing
end

@inline function _print_converged!()
    _FAST_USE_LOGS && println("Stopping criterion met.")
    return nothing
end

function _print_summary!(
    tag::AbstractString,
    algo_name::AbstractString,
    result;
    rho_lines = (),
)
    _FAST_USE_LOGS || return nothing
    println("\n$tag Summary:")
    println("→ Total generation cost = $(result.total_cost)")
    for line in rho_lines
        println(line)
    end
    println("Iterations: $(result.iter_done)")
    println("$algo_name time: $(_pretty_time(result.elapsed_s))")
    return nothing
end

# ============================================================
# ADMM
# ============================================================

_print_tag(::DCStateWorkspace) = "[ADMM-state/DC]"
_print_tag(::DCPowerWorkspace) = "[ADMM-power/DC]"
_print_tag(::DCHybridWorkspace) = "[ADMM-hybrid/DC]"
_print_tag(::ACStateWorkspace) = "[ADMM-state/AC]"
_print_tag(::ACPowerWorkspace) = "[ADMM-power/AC]"
_print_tag(::ACHybridWorkspace) = "[ADMM-hybrid/AC]"
@inline function _socwr_zone_label(ws)
    ps = ws.partitions.parts
    return isempty(ps) ? "SOCWR" : string(nameof(ps[1].zone_model))
end
_print_tag(ws::SOCWRStateWorkspace) = "[ADMM-state/$(_socwr_zone_label(ws))]"
_print_tag(ws::SOCWRPowerWorkspace) = "[ADMM-power/$(_socwr_zone_label(ws))]"
_print_tag(ws::SOCWRHybridWorkspace) = "[ADMM-hybrid/$(_socwr_zone_label(ws))]"
_print_tag(::AbstractADMMWorkspace) = "[ADMM]"

@inline function _print_admm_residuals!(ws::AbstractADMMWorkspace)
    _FAST_USE_LOGS && println("‖r‖ = $(ws.r_primal) | eps_pri = $(ws.eps_primal)")
    _FAST_USE_LOGS && println("‖s‖ = $(ws.r_dual) | eps_du  = $(ws.eps_dual)")
    return nothing
end

# ============================================================
# ALADIN
# ============================================================

_print_tag(::DCStateWorkspaceALADIN) = "[ALADIN-state/DC]"
_print_tag(::DCPowerWorkspaceALADIN) = "[ALADIN-power/DC]"
_print_tag(::DCHybridWorkspaceALADIN) = "[ALADIN-hybrid/DC]"

_print_tag(::ACStateWorkspaceALADIN) = "[ALADIN-state/AC]"
_print_tag(::ACPowerWorkspaceALADIN) = "[ALADIN-power/AC]"
_print_tag(::ACHybridWorkspaceALADIN) = "[ALADIN-hybrid/AC]"

_print_tag(ws::SOCWRStateWorkspaceALADIN) = "[ALADIN-state/$(_socwr_zone_label(ws))]"
_print_tag(ws::SOCWRPowerWorkspaceALADIN) = "[ALADIN-power/$(_socwr_zone_label(ws))]"
_print_tag(ws::SOCWRHybridWorkspaceALADIN) = "[ALADIN-hybrid/$(_socwr_zone_label(ws))]"

_print_tag(::AbstractALADINWorkspace) = "[ALADIN]"

@inline function _print_aladin_residuals!(ws::AbstractALADINWorkspace)
    _FAST_USE_LOGS && println(
        "‖r_primal‖      = $(ws.r_primal)      | eps_pri  = $(ws.eps_primal)",
    )
    _FAST_USE_LOGS && println(
        "‖r_step_primal‖ = $(ws.r_step_primal) | eps_step = $(ws.eps_step_primal)",
    )
    _FAST_USE_LOGS && println(
        "‖r_dual‖        = $(ws.r_dual)        | eps_dual = $(ws.eps_dual)",
    )
    _FAST_USE_LOGS && println(
        "‖Δλ‖            = $(ws.lambda_step)   (diagnostic)",
    )
    return nothing
end
