_admm_print_tag(::DCStateWorkspace) = "[ADMM-state/DC]"
_admm_print_tag(::DCPowerWorkspace) = "[ADMM-power/DC]"
_admm_print_tag(::DCHybridWorkspace) = "[ADMM-hybrid/DC]"
_admm_print_tag(::ACStateWorkspace) = "[ADMM-state/AC]"
_admm_print_tag(::ACPowerWorkspace) = "[ADMM-power/AC]"
_admm_print_tag(::ACHybridWorkspace) = "[ADMM-hybrid/AC]"
_admm_print_tag(::SOCWRStateWorkspace) = "[ADMM-state/SOCWR]"
_admm_print_tag(::SOCWRPowerWorkspace) = "[ADMM-power/SOCWR]"
_admm_print_tag(::SOCWRHybridWorkspace) = "[ADMM-hybrid/SOCWR]"
_admm_print_tag(::AbstractADMMWorkspace) = "[ADMM]"

@inline function _print_admm_residuals!(ws::AbstractADMMWorkspace)
    _FAST_USE_LOGS && println("‖r‖ = $(ws.r_primal) | eps_pri = $(ws.eps_primal)")
    _FAST_USE_LOGS && println("‖s‖ = $(ws.r_dual) | eps_du  = $(ws.eps_dual)")
    return nothing
end
