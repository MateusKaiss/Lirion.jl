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
    println("\n$tag Resumo:")
    println("→ Total generation cost = $(result.total_cost)")
    for line in rho_lines
        println(line)
    end
    println("Iterações realizadas: $(result.iter_done)")
    println("Tempo $algo_name: $(_pretty_time(result.elapsed_s))")
    return nothing
end
