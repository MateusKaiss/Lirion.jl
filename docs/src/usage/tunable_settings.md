# Tunable settings

!!! note "This page covers centralized runs"
    These settings apply to a plain centralized solve. For a distributed run, the penalty, coupling, sectioning, and iteration limit live under [ADMM settings](@ref) and [ALADIN settings](@ref).

In a centralized run the main thing you tune is the solver.

## The solver

`solve` picks a default solver from the model. Here is the full list.

| Model     | Default solver |
| --------- | -------------- |
| `DC`      | OSQP           |
| `AC`      | Ipopt          |
| `ACR`     | Ipopt          |
| `IVR`     | Ipopt          |
| `LPAC`    | Ipopt          |
| `SOCWR`   | Clarabel       |
| `SOCWRNL` | Ipopt          |
| `QC`      | Ipopt          |
| `SDPWRM`  | Clarabel       |
| `SOCBF`   | Clarabel       |
| `BFA`     | Clarabel       |
| `NFA`     | OSQP           |

Override the default, or set options on it, with the `optimizer` keyword. It takes any JuMP optimizer, so you can wrap it with `optimizer_with_attributes` to pass solver options.

```julia
using JuMP, Ipopt

result = solve("pglib_opf_case24_ieee_rts.m";
    model = ac(),
    optimizer = optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000, "tol" => 1e-8),
)
```

You can point `optimizer` at any solver that handles the problem class the model produces (LP, QP, QCQP, SOCP, SDP, or NLP). DC is a linear program, SOCWR is a second-order cone program, and AC is a nonlinear program, so a solver that fits one may not fit another.
