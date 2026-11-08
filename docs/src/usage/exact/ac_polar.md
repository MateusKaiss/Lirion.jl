# Running AC Polar

The AC model in polar coordinates, run with `model = ac()`.
For the formulation, see [AC Polar](@ref).

All examples use `pglib_opf_case24_ieee_rts.m` from the [PGLib-OPF library](https://github.com/power-grid-lib/pglib-opf).

```julia
using Lirion
result = solve("test/matpower_cases/pglib_opf_case24_ieee_rts.m"; model = ac())

result.objective   # optimal cost
result.status      # solver status
result.bus         # per-bus solution (DataFrame)
result.gen         # generator dispatch (DataFrame)
```

## Choosing a solver

AC Polar runs on Ipopt by default. To swap in another solver, or to set options on the default one, pass it through the `optimizer` keyword, which takes any JuMP optimizer.

```julia
using JuMP, Ipopt

result = solve("test/matpower_cases/pglib_opf_case24_ieee_rts.m";
    model = ac(),
    optimizer = optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000),
)
```

For what `result` holds and how to work with it, see [Results](@ref).
