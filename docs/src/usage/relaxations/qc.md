# Running QC

The quadratic-convex relaxation of the AC power flow, run with `model = qc()`. It is a relaxation like SOCWR, so its optimum is a lower bound on the true AC cost; it tightens the voltage-product blocks with McCormick envelopes instead of a single cone. For the formulation, see [QC (Quadratic Convex)](@ref).

All examples use `pglib_opf_case24_ieee_rts.m` from the [PGLib-OPF library](https://github.com/power-grid-lib/pglib-opf).

```julia
using Lirion

result = solve("test/matpower_cases/pglib_opf_case24_ieee_rts.m"; model = qc())

result.objective   # optimal cost (a lower bound on the AC optimum)
result.status      # solver status
result.bus         # per-bus solution
result.gen         # generator dispatch
```

## Choosing a solver

QC runs on Ipopt by default. To swap in another solver, or to set options on the default one, pass it through the `optimizer` keyword, which takes any JuMP optimizer.

```julia
using JuMP, Ipopt

result = solve("test/matpower_cases/pglib_opf_case24_ieee_rts.m";
    model = qc(),
    optimizer = optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000),
)
```

For what `result` holds and how to work with it, see [Results](@ref).
