_print_tag(::DCStateWorkspace) = "[ADMM-state/DC]"
_print_tag(::DCPowerWorkspace) = "[ADMM-power/DC]"
_print_tag(::DCHybridWorkspace) = "[ADMM-hybrid/DC]"
_print_tag(::ACStateWorkspace) = "[ADMM-state/AC]"
_print_tag(::ACPowerWorkspace) = "[ADMM-power/AC]"
_print_tag(::ACHybridWorkspace) = "[ADMM-hybrid/AC]"
_print_tag(::SOCWRStateWorkspace) = "[ADMM-state/SOCWR]"
_print_tag(::SOCWRPowerWorkspace) = "[ADMM-power/SOCWR]"
_print_tag(::SOCWRHybridWorkspace) = "[ADMM-hybrid/SOCWR]"
_print_tag(::AbstractADMMWorkspace) = "[ADMM]"

@inline function _print_admm_residuals!(ws::AbstractADMMWorkspace)
    _FAST_USE_LOGS && println("‖r‖ = $(ws.r_primal) | eps_pri = $(ws.eps_primal)")
    _FAST_USE_LOGS && println("‖s‖ = $(ws.r_dual) | eps_du  = $(ws.eps_dual)")
    return nothing
end
