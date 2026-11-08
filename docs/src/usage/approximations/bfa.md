# Running BFA

The branch-flow (LinDistFlow) approximation, run with `model = bfa()`. For the formulation, see [BFA (LinDistFlow)](@ref).

All examples use `pglib_opf_case24_ieee_rts.m` from the [PGLib-OPF library](https://github.com/power-grid-lib/pglib-opf).

```julia
using Lirion

result = solve("test/matpower_cases/pglib_opf_case24_ieee_rts.m"; model = bfa())

result.objective   # optimal cost
result.status      # solver status
result.bus         # squared voltage magnitude per bus
result.gen         # generator dispatch
```

## Choosing a solver

BFA runs on Clarabel by default. To swap in another solver, or to set options on the default one, pass it through the `optimizer` keyword, which takes any JuMP optimizer.

```julia
using JuMP, Clarabel

result = solve("test/matpower_cases/pglib_opf_case24_ieee_rts.m";
    model = bfa(),
    optimizer = optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000),
)
```

For what `result` holds and how to work with it, see [Results](@ref).
