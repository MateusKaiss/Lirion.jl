# Running IVR

The current–voltage model in rectangular coordinates, run with `model = ivr()`.
For the formulation, see [IVR](@ref).

All examples use `pglib_opf_case57_ieee.m` from the [PGLib-OPF library](https://github.com/power-grid-lib/pglib-opf).

```julia
using Lirion
result = solve("pglib_opf_case57_ieee.m"; model = ivr())

result.objective   # optimal cost
result.status      # solver status
result.bus         # per-bus solution (DataFrame)
result.gen         # generator dispatch (DataFrame)
```

## Choosing a solver

IVR runs on Ipopt by default. To swap in another solver, or to set options on the default one, pass it through the `optimizer` keyword, which takes any JuMP optimizer.

```julia
using JuMP, Ipopt

result = solve("pglib_opf_case57_ieee.m";
    model = ivr(),
    optimizer = optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000),
)
```

For what `result` holds and how to work with it, see [Results](@ref).
