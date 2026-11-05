```@raw html
<img src="assets/logo.png" alt="Lirion.jl logo"
     style="max-width: 520px; width: 100%; margin: 0 auto 1.25rem; display: block;" />

```


```@meta
# This prevents certain rendering issues in Documenter
```

# Introduction

Lirion is a Julia package for Optimal Power Flow (OPF). You can solve a case centrally, as a single optimization problem, or you can distribute it. In the distributed setting the grid is split into zones. Each zone solves its own subproblem, and the zones coordinate only through the variables they share at their boundaries.

Centralized, Lirion covers a wide range of models. The exact AC formulation comes in polar and rectangular coordinates, and in the current-voltage form. The convex relaxations are SOCWR (Jabr), QC, and SDP. The linear and approximate models are DC, LPAC, NFA, and BFA (LinDistFlow). For the distributed setting, Lirion provides two coordination methods, **ADMM** and **ALADIN**, currently for the DC, AC, and second-order-cone models. Either way, you choose how the grid is partitioned and which variables the zones agree on.

## How this documentation is structured

The documentation has two main parts, plus supporting material:

- **[Getting Started](installation.md)** is the practical side: how to run each centralized model, how to pick a solver, the tunable settings, the distributed workflow with ADMM and ALADIN, and how to read the results.
- **[Theory](teoria/fundamentos.md)** is the mathematics behind every formulation, centralized and distributed, including how a distributed problem is cut into zones and coupled at the boundary.
- **[Performance](performance.md)** collects the benchmark results, and **[References](references.md)** is the bibliography.

## Centralized Models

Lirion implements a broad family of centralized OPF formulations, grouped by how they treat the AC physics:

**Exact AC**
- AC OPF in polar coordinates (`AC`)
- AC OPF in rectangular coordinates (`ACR`)
- Current–voltage rectangular (`IVR`)

**Convex relaxations**
- Second-order cone / Jabr (`SOCWR`), and its smooth NLP variant (`SOCWRNL`)
- Quadratic convex (`QC`)
- Semidefinite (`SDPWRM`)
- SOC branch-flow (`SOCBF`)

**Linear and approximate**
- DC OPF (`DC`)
- LPAC (`LPAC`)
- Network-flow approximation (`NFA`)
- LinDistFlow / branch-flow approximation (`BFA`)

---

## Distributed Formulations

Beyond centralized models, Lirion provides distributed solvers based on the **Alternating Direction Method of Multipliers (ADMM)** and the **Augmented Lagrangian Alternating Direction Inexact Newton (ALADIN)** method.

Distributed formulations are structured along three main dimensions:

- The base model — currently DC, AC, and the smooth second-order-cone relaxation (`SOCWRNL`);[^socwr]
- The choice of **coupling variables** (e.g. voltage–angle or power);
- The **system partitioning strategy** (tie-line, split-bus, or dummy-bus sectioning).

[^socwr]: The conic second-order-cone model (`SOCWR`) is also available and works under ADMM. Under ALADIN the second-order-cone solvers currently limit it, so for distributed solves it is superseded by the smooth `SOCWRNL`; those limitations are discussed in [Sectioning and coupling](teoria/distribuido/seccionamento_acoplamento.md#SOCWR,-and-its-limitations).

---
