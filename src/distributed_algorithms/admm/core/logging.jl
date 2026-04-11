@inline function _pretty_time(sec::Real)
    if sec < 60
        return "$(round(sec; digits = 3)) s"
    elseif sec < 3600
        return "$(round(sec / 60; digits = 2)) min"
    else
        return "$(round(sec / 3600; digits = 2)) h"
    end
end

# ------------------------------------------------------------
# PRINT TAG
# ------------------------------------------------------------
_admm_print_tag(::DCStateWorkspace) = "[ADMM-state/DC]"
_admm_print_tag(::DCPowerWorkspace) = "[ADMM-power/DC]"
_admm_print_tag(::DCHybridWorkspace) = "[ADMM-hybrid/DC]"
_admm_print_tag(::ACStateWorkspace) = "[ADMM-state/AC]"
_admm_print_tag(::ACPowerWorkspace) = "[ADMM-power/AC]"
_admm_print_tag(::ACHybridWorkspace) = "[ADMM-hybrid/AC]"
_admm_print_tag(::SOCPStateWorkspace) = "[ADMM-state/SOCP]"
_admm_print_tag(::SOCPPowerWorkspace) = "[ADMM-power/SOCP]"
_admm_print_tag(::SOCPHybridWorkspace) = "[ADMM-hybrid/SOCP]"
_admm_print_tag(::AbstractADMMWorkspace) = "[ADMM]"

# ------------------------------------------------------------
# ITERATION PRINTS
# ------------------------------------------------------------
@inline function _print_admm_iter_header!(tag::AbstractString, iter::Int)
    _FAST_USE_LOGS && println("\n$tag Iteration $iter")
    return nothing
end

@inline function _print_admm_residuals!(ws::AbstractADMMWorkspace)
    _FAST_USE_LOGS && println("‖r‖ = $(ws.r_primal) | eps_pri = $(ws.eps_primal)")
    _FAST_USE_LOGS && println("‖s‖ = $(ws.r_dual) | eps_du  = $(ws.eps_dual)")
    return nothing
end

@inline function _print_admm_converged!()
    _FAST_USE_LOGS && println("Stopping criterion met.")
    return nothing
end

# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------
function _print_admm_summary!(
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
    println("Tempo ADMM: $(_pretty_time(result.elapsed_s))")
    return nothing
end
