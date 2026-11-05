# Running SOCWR

The Jabr second-order-cone relaxation of the AC power flow, written in W-space and run with `model = socwr()`. Like the other relaxations it enlarges the feasible set to make the problem convex, so its optimum is a lower bound on the true AC cost rather than a feasible AC operating point. For the formulation, see [SOCWR (Jabr conic relaxation)](@ref).

All examples use `pglib_opf_case24_ieee_rts.m` from the [PGLib-OPF library](https://github.com/power-grid-lib/pglib-opf).

```julia
using Lirion

result = solve("pglib_opf_case24_ieee_rts.m"; model = socwr())

result.objective   # optimal cost (a lower bound on the AC optimum)
result.status      # solver status
result.bus         # per-bus solution
result.gen         # generator dispatch
```

## Choosing a solver

SOCWR runs on Clarabel by default. To swap in another solver, or to set options on the default one, pass it through the `optimizer` keyword, which takes any JuMP optimizer.

```julia
using JuMP, Clarabel

result = solve("pglib_opf_case24_ieee_rts.m";
    model = socwr(),
    optimizer = optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000),
)
```

For what `result` holds and how to work with it, see [Results](@ref).
