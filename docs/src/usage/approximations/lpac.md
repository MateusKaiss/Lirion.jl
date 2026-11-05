# Running LPAC

The LPAC model, which gives voltage magnitude and reactive power back to a linear model, run with `model = lpac()`. For the formulation, see [LPAC](@ref).

All examples use `pglib_opf_case57_ieee.m` from the [PGLib-OPF library](https://github.com/power-grid-lib/pglib-opf).

```julia
using Lirion

result = solve("pglib_opf_case57_ieee.m"; model = lpac())

result.objective   # optimal cost
result.status      # solver status
result.bus         # voltage and angle per bus
result.gen         # generator dispatch (active and reactive)
```

## Choosing a solver

LPAC runs on Ipopt by default. To swap in another solver, or to set options on the default one, pass it through the `optimizer` keyword, which takes any JuMP optimizer.

```julia
using JuMP, Ipopt

result = solve("pglib_opf_case57_ieee.m";
    model = lpac(),
    optimizer = optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000),
)
```

For what `result` holds and how to work with it, see [Results](@ref).
