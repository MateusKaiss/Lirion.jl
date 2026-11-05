# Running ACR

The AC model in rectangular voltage coordinates, run with `model = acr()`.
For the formulation, see [AC Rectangular (ACR)](@ref).

All examples use `pglib_opf_case24_ieee_rts.m` from the [PGLib-OPF library](https://github.com/power-grid-lib/pglib-opf).

```julia
using Lirion
result = solve("pglib_opf_case24_ieee_rts.m"; model = acr())

result.objective   # optimal cost
result.status      # solver status
result.bus         # per-bus solution (DataFrame)
result.gen         # generator dispatch (DataFrame)
```

## Choosing a solver

ACR runs on Ipopt by default. To swap in another solver, or to set options on the default one, pass it through the `optimizer` keyword, which takes any JuMP optimizer.

```julia
using JuMP, Ipopt

result = solve("pglib_opf_case24_ieee_rts.m";
    model = acr(),
    optimizer = optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000),
)
```

For what `result` holds and how to work with it, see [Results](@ref).
