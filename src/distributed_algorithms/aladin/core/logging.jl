_print_tag(::DCStateWorkspaceALADIN) = "[ALADIN-state/DC]"
_print_tag(::DCPowerWorkspaceALADIN) = "[ALADIN-power/DC]"
_print_tag(::DCHybridWorkspaceALADIN) = "[ALADIN-hybrid/DC]"

_print_tag(::ACStateWorkspaceALADIN) = "[ALADIN-state/AC]"
_print_tag(::ACPowerWorkspaceALADIN) = "[ALADIN-power/AC]"
_print_tag(::ACHybridWorkspaceALADIN) = "[ALADIN-hybrid/AC]"

_print_tag(::SOCWRStateWorkspaceALADIN) = "[ALADIN-state/SOCWR]"
_print_tag(::SOCWRPowerWorkspaceALADIN) = "[ALADIN-power/SOCWR]"
_print_tag(::SOCWRHybridWorkspaceALADIN) = "[ALADIN-hybrid/SOCWR]"

_print_tag(::AbstractALADINWorkspace) = "[ALADIN]"

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
