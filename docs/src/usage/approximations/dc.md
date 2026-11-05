# Running DC

The linear DC power flow, run with `model = dc()`. For the formulation, see [DC](@ref).

All examples use `pglib_opf_case24_ieee_rts.m` from the [PGLib-OPF library](https://github.com/power-grid-lib/pglib-opf).

```julia
using Lirion

result = solve("pglib_opf_case24_ieee_rts.m"; model = dc())

result.objective   # optimal cost
result.status      # solver status
result.bus         # bus angles
result.gen         # generator dispatch
```

## Choosing a solver

DC runs on OSQP by default. To swap in another solver, or to set options on the default one, pass it through the `optimizer` keyword, which takes any JuMP optimizer.

```julia
using JuMP, OSQP

result = solve("pglib_opf_case24_ieee_rts.m";
    model = dc(),
    optimizer = optimizer_with_attributes(OSQP.Optimizer, "max_iter" => 1000),
)
```

For what `result` holds and how to work with it, see [Results](@ref).
