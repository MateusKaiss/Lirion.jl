# Running SDPWRM

The semidefinite W-matrix relaxation of the AC power flow, run with `model = sdpwrm()`. It is a relaxation like SOCWR and QC, so its optimum is a lower bound on the true AC cost; it collects the voltage-product blocks into one matrix and imposes a single positive-semidefinite condition over the whole network. For the formulation, see [SDPWRM](@ref).

All examples use `pglib_opf_case24_ieee_rts.m` from the [PGLib-OPF library](https://github.com/power-grid-lib/pglib-opf).

```julia
using Lirion

result = solve("pglib_opf_case24_ieee_rts.m"; model = sdpwrm())

result.objective   # optimal cost (a lower bound on the AC optimum)
result.status      # solver status
result.bus         # per-bus solution
result.gen         # generator dispatch
```

## Choosing a solver

SDPWRM runs on Clarabel by default. To swap in another solver, or to set options on the default one, pass it through the `optimizer` keyword, which takes any JuMP optimizer.

```julia
using JuMP, Clarabel

result = solve("pglib_opf_case24_ieee_rts.m";
    model = sdpwrm(),
    optimizer = optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000),
)
```

For what `result` holds and how to work with it, see [Results](@ref).
