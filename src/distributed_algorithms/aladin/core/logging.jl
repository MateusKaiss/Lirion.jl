_aladin_print_tag(::DCStateWorkspaceALADIN) = "[ALADIN-state/DC]"
_aladin_print_tag(::DCPowerWorkspaceALADIN) = "[ALADIN-power/DC]"
_aladin_print_tag(::DCHybridWorkspaceALADIN) = "[ALADIN-hybrid/DC]"

_aladin_print_tag(::ACStateWorkspaceALADIN) = "[ALADIN-state/AC]"
_aladin_print_tag(::ACPowerWorkspaceALADIN) = "[ALADIN-power/AC]"
_aladin_print_tag(::ACHybridWorkspaceALADIN) = "[ALADIN-hybrid/AC]"

_aladin_print_tag(::SOCWRStateWorkspaceALADIN) = "[ALADIN-state/SOCWR]"
_aladin_print_tag(::SOCWRPowerWorkspaceALADIN) = "[ALADIN-power/SOCWR]"
_aladin_print_tag(::SOCWRHybridWorkspaceALADIN) = "[ALADIN-hybrid/SOCWR]"

_aladin_print_tag(::AbstractALADINWorkspace) = "[ALADIN]"

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
