# Running NFA

The network flow approximation, run with `model = nfa()`. For the formulation, see [NFA](@ref).

All examples use `pglib_opf_case24_ieee_rts.m` from the [PGLib-OPF library](https://github.com/power-grid-lib/pglib-opf).

```julia
using Lirion

result = solve("test/matpower_cases/pglib_opf_case24_ieee_rts.m"; model = nfa())

result.objective   # optimal cost
result.status      # solver status
result.bus         # per-bus results
result.gen         # generator dispatch
```

NFA carries active power only, so the reactive columns of `result.gen` and `result.bus` are not populated.

## Choosing a solver

NFA runs on OSQP by default. To swap in another solver, or to set options on the default one, pass it through the `optimizer` keyword, which takes any JuMP optimizer.

```julia
using JuMP, OSQP

result = solve("test/matpower_cases/pglib_opf_case24_ieee_rts.m";
    model = nfa(),
    optimizer = optimizer_with_attributes(OSQP.Optimizer, "max_iter" => 1000),
)
```

For what `result` holds and how to work with it, see [Results](@ref).
