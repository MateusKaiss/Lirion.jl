# Installation Guide

This guide explains how to install Julia, Lirion, and the required solvers.

---

## Install Julia

Lirion is written in [Julia](https://julialang.org).

Download and install Julia from:
[https://julialang.org/downloads/](https://julialang.org/downloads/)


If you are new to Julia, we recommend reading the official
[Getting Started guide](https://docs.julialang.org/en/v1/manual/getting-started/).


This tutorial explains how to install Julia, start the REPL, manage packages, and run scripts.

---

## Install Lirion

Lirion is installed using Julia’s built-in package manager.

Start Julia and run:

```julia
julia> import Pkg

julia> Pkg.add("Lirion")
```

## Default Solvers

Lirion uses external optimization solvers to solve OPF models.

We highly recommend that users study and select a solver that best fits their problem characteristics, numerical requirements, and performance expectations. See [Install a Solver and Other Popular Solvers](@ref) below.

However, if no solver is explicitly specified, Lirion automatically selects a default solver depending on the OPF formulation:

| OPF Model | Default Solver | Problem Class |
|-----------|---------------|---------------|
| **DC OPF** | `OSQP` | Convex Quadratic Programming |
| **AC OPF** | `Ipopt` | Nonlinear Programming (NLP) |
| **SOCWR OPF** | `Clarabel` | Conic Optimization (SOCWR) |

These default solvers are installed automatically as dependencies when Lirion is installed.

For distributed solves, ALADIN adds one more default: the coordination QP that links the zones runs on `Clarabel`, unless you override it with `qp_optimizer` (see [ALADIN settings](usage/distributed/aladin.md)).

---

## Install a Solver and Other Popular Solvers

To install any solver package, use Julia’s built-in package manager:

```julia
import Pkg
Pkg.add("HiGHS")
```
In the following sections of this documentation, we will show how to use both the default and custom-installed solvers with Lirion’s centralized and distributed OPF formulations.

Below is a list of well-known solvers compatible with the Julia optimization ecosystem:

| Solver | Julia Package | License | Typical Use |
|--------|--------------|----------|------------|
| HiGHS | `HiGHS.jl` | MIT | LP / QP (fast open-source alternative) |
| Gurobi | `Gurobi.jl` | Commercial | LP / QP / MIP / SOCWR |
| CPLEX | `CPLEX.jl` | Commercial | LP / QP / MIP |
| Mosek | `MosekTools.jl` | Commercial | LP / QP / SOCWR / SDP |
| SCS | `SCS.jl` | MIT | Large-scale conic problems |
| COSMO | `COSMO.jl` | Apache | Conic optimization |
| KNITRO | `KNITRO.jl` | Commercial | Large-scale nonlinear optimization |
| MadNLP | `MadNLP.jl` | MIT | Large-scale nonlinear programming |

---

### Complete Solver List

For a comprehensive and continuously updated list of supported solvers in the Julia ecosystem, see the official JuMP documentation:

[JuMP – Supported Solvers](https://jump.dev/JuMP.jl/stable/installation/#supported-solvers)
