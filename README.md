<p align="center">
  <img src="docs/src/assets/logo.png" alt="Lirion.jl" 
  width="600"/>
</p>

<h3 align="center">
  Large-scale Integrated Research on Innovative Optimization of Networks

</h3>

<p align="center">
  A Julia package for solving Optimal Power Flow (OPF) problems using both centralized and distributed algorithms.
</p>

<p align="center">
  <a href="https://mateuskaiss.github.io/Lirion.jl/"><img src="https://img.shields.io/badge/docs-blue.svg" alt="Documentation"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
</p>

## Why Lirion?  

Lirion was built to bridge implementation clarity and theoretical rigor in Optimal Power Flow research.

It provides:

- **Twelve centralized OPF formulations**:
  - **DC** — linear active-power approximation (LP)
  - **AC** — full nonlinear AC in polar coordinates (NLP)
  - **ACR** — full nonlinear AC in rectangular coordinates (NLP)
  - **IVR** — current-voltage rectangular formulation (NLP)
  - **LPAC** — LP-based AC approximation (QCQP)
  - **SOCWR** — W-rectangular second-order cone relaxation (SOCP)
  - **SOCWRNL** — smooth nonlinear form of the SOCWR relaxation (NLP)
  - **QC** — quadratic-convex relaxation (QCQP)
  - **SDPWRM** — semidefinite W-matrix relaxation (SDP)
  - **SOCBF** — branch-flow second-order cone relaxation (SOCP)
  - **BFA** — linearized branch-flow approximation / LinDistFlow (LP)
  - **NFA** — network-flow approximation, active-power only (LP)
- **Distributed OPF** — two decomposition algorithms over DC / AC / SOCWR / SOCWRNL, each with its own variants:
  - **ADMM** — alternating direction method of multipliers
    - fixed or **adaptive penalty (ρ)** with residual balancing
  - **ALADIN** — augmented Lagrangian based alternating direction inexact Newton
    - **proximal mode** (`:primal` / `:pflow`) and **stopping criterion** (`:boyd` / `:hfd`)
    - optional **adaptive ρ / μ** growth
  - **Coupling variables** — which boundary quantities neighboring zones exchange:
    - **`:state`** — voltage state only: θ (DC), θ + V (AC), or v², R, T (SOCWR / SOCWRNL)
    - **`:power`** — branch active/reactive flows only
    - **`:hybrid`** — both state and power.
  - **Network sectioning** — how the grid is split into zones:
    - **`:tie_line`** — zones overlap on shared tie-line buses
    - **`:dummy_bus`** — each tie-line is split by an intermediate dummy bus
    - **`:split_bus`** — boundary buses are cloned across zones with injection consensus

> **Supported combinations.**
> - **`:state`** — works only with `:tie_line`
> - **`:power`** — works with `:tie_line`, or with `:dummy_bus` for DC / AC / SOCWRNL (the conic **SOCWR** does not support `:dummy_bus` under ALADIN)
> - **`:hybrid`** — works with all three (`:tie_line`, `:dummy_bus`, `:split_bus`), except the conic **SOCWR**, which does not support `:dummy_bus` under ALADIN (its smooth sibling **SOCWRNL** does)

---

## Installation

### Install Julia

Lirion is written in [Julia](https://julialang.org).

If you do not have Julia installed, download it from:

https://julialang.org/downloads/

---

### Install Lirion

Lirion can be installed using Julia’s built-in package manager.

Start Julia and run:

```julia
import Pkg
Pkg.add("Lirion")
```

Then load the package:

```julia
using Lirion
```

---

## Quick Examples

The centralized and distributed examples below use `pglib_opf_case24_ieee_rts.m` from
the **PGLib-OPF** benchmark library, shipped with Lirion. They assume you cloned the
repository and started Julia from its root (`git clone … && cd Lirion.jl && julia
--project`); otherwise pass an absolute path to any MATPOWER file. The documentation
covers the formulations and the distributed algorithms in depth.

If you don't pass a solver, Lirion uses a per-formulation default:

- **OSQP** — DC, NFA
- **Ipopt** — AC, ACR, IVR, LPAC, QC, SOCWRNL
- **Clarabel** — SOCWR, SOCBF, SDPWRM, BFA

---

### Centralized

A centralized solve takes a case file and a model:

```julia
using Lirion

result = solve("test/matpower_cases/pglib_opf_case24_ieee_rts.m"; model = ac())

result.objective
```

Switch the formulation through the `model` argument; any of the twelve works:

```julia
model = dc()       # DC linear approximation        (LP)
model = ac()       # AC polar (V, θ)                (NLP)
model = acr()      # AC rectangular (vr, vi)        (NLP)
model = ivr()      # Current-voltage rectangular    (NLP)
model = lpac()     # LP-based AC approximation      (QCQP)
model = socwr()    # W-rectangular SOCP relaxation  (SOCP)
model = socwrnl()  # Smooth-NLP SOCWR relaxation    (NLP)
model = qc()       # Quadratic-convex relaxation    (QCQP)
model = sdpwrm()   # Semidefinite W-matrix          (SDP)
model = socbf()    # Branch-flow SOCP relaxation    (SOCP)
model = bfa()      # Linearized branch flow         (LP)
model = nfa()      # Network-flow approximation     (LP)
```

### Distributed

A distributed run uses the same `solve`; swap `centralized()` for `admm(...)` or
`aladin(...)` and pick the **coupling** and the **sectioning**. The case must carry
zones (a MATPOWER zone column) — see **Automatic partitioning** below to make one from
any case; the example uses a shipped three-zone partition of case24.

```julia
using Lirion

result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = dc(),
    algorithm = admm(
        coupling = :hybrid, sectioning = :tie_line,
        rho = 1e3,                       # penalty weight
        eps_abs = 1e-4, eps_rel = 1e-4,  # convergence tolerance
        max_iter = 1500,
    ),
)

result.status        # :converged / :stalled / :iteration_limit
result.objective
result.diagnostics   # iteration history (residuals, ρ, ...)
```

Swap `admm` for `aladin` to run ALADIN, and change `coupling` / `sectioning` to set
what the zones exchange and how the grid is cut:

```julia
# Algorithm
algorithm = admm(...)        # ADMM
algorithm = aladin(...)      # ALADIN

# Coupling — which boundary quantities the zones agree on
coupling = :state            # voltage state (θ / V / v²,R,T)
coupling = :power            # branch active/reactive flows
coupling = :hybrid           # both

# Network sectioning — how the grid is split into zones
sectioning = :tie_line       # zones share boundary buses
sectioning = :dummy_bus      # tie-lines split by dummy buses
sectioning = :split_bus      # boundary buses cloned across zones
```

The distributed solvers run on DC, AC, SOCWR, and SOCWRNL. Not every coupling pairs
with every sectioning, and the penalty (`rho`) and tolerance settings decide whether,
and how closely, a run converges, so they are not innocent defaults. The Distributed
pages in the documentation cover the accepted combinations, the per-model penalties,
and the convergence settings; AC, SOCWR, and SOCWRNL runs in particular need larger
iteration budgets than DC.

#### Automatic partitioning

To cut a single-zone case into zones automatically, use impedance-weighted **Metis**
partitioning (balanced zones, connectivity-repaired). It needs the optional extension
dependencies:

```julia
using Lirion, Metis, Graphs, SimpleWeightedGraphs

# Split case118 into 4 zones and write the zoned case next to the input
zoned = create_partition("test/matpower_cases/pglib_opf_case118_ieee.m"; k = 4)
```

The path it returns is then solved like any distributed case:

```julia
result = solve(zoned;
    model = dc(),
    algorithm = admm(
        coupling = :hybrid, sectioning = :tie_line,
        rho = 1e4, eps_abs = 1e-4, eps_rel = 1e-4, max_iter = 1000,
    ),
)
```

## Results

`solve` returns a `LirionResult`: a short summary plus the solution
as `DataFrame`s following the MATPOWER schema.

```julia
result.objective   # total generation cost
result.status      # :optimal / :local_optimal (centralized)
result.bus         # per-bus voltages / angles
result.gen         # generator dispatch (pg, qg)
result.branch      # branch flows (pf, qf, pt, qt)
```

A distributed run fills in the per-zone view and the convergence history, and its
`status` reports the algorithm's own termination:

```julia
result.status      # :converged / :stalled / :iteration_limit (distributed)
result.zones       # per-zone (bus, gen, branch, cost)
result.consensus   # boundary-variable agreement across zones
result.diagnostics # iteration history (residuals, ρ, ...)
```

Write everything to CSV (needs `using CSV`):

```julia
Lirion.to_csvs(result, "results/")
```

## Performance & Threading

Lirion has two very different workloads, and each wants a different threading setup.
The rule is **never nest parallelism**: go *deep* (centralized: one solver, BLAS on all
cores) or *wide* (distributed: many zones, one BLAS thread each).

**Centralized.** Leave the defaults. Don't set `OPENBLAS_NUM_THREADS=1`, or the dense
factorizations lose their threads.

**Distributed.** Set `JULIA_NUM_THREADS` before you launch Julia (it can't change at
runtime; the `-t N` flag works too), and pin BLAS to one thread:

```bash
export JULIA_NUM_THREADS=$(getconf _NPROCESSORS_ONLN)   # macOS: $(sysctl -n hw.ncpu)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
julia --project
```

On Windows PowerShell, set `$env:JULIA_NUM_THREADS = $env:NUMBER_OF_PROCESSORS` and the
three BLAS variables to `"1"`. These counts are logical CPUs (hyperthreads included),
which is what you want here: one Julia thread per zone.

Lirion scopes BLAS to one thread for the whole distributed solve, so the variables above
are belt-and-suspenders. Call `Lirion.threadinfo()` to print the Julia, BLAS, and CPU
thread counts and check it took. For the rationale and verification tips, see the
**Performance & Threading** page in the documentation.

---

## Documentation

The [Lirion documentation](https://mateuskaiss.github.io/Lirion.jl/) comes in two parts.

The **user guide** is the practical half: installation, the centralized and
distributed workflows, solver configuration, worked examples, and the API reference.

The **theory** half builds each model from the ground up instead of just listing it.
It grows out of the master's research (still ongoing) of the lead developer,
Mateus Kaiss, and works through the formulations from first principles:

- the **twelve OPF formulations** (DC, AC, ACR, IVR, LPAC, SOCWR, SOCWRNL, QC, SDPWRM,
  SOCBF, BFA, NFA), and where each one comes from;
- the **distributed methods** — how the grid is sectioned into zones and what the
  zones exchange at their boundaries;
- the **coupling-variable design** and the reasoning behind the network decomposition.


---

## Scientific Outputs & Deliverables

This software is part of an ongoing research program. The following works have been produced as part of this project's development:

1. **BSc Thesis** — *Desenvolvimento de Ferramenta para Determinação do Fluxo de Potência Ótimo em Redes Elétricas*. Mateus Kaiss, UFPR (2025).
2. **White Paper** — *Lirion.jl: A Unified Framework for Centralized and Certified Distributed Optimal Power Flow*. Mateus Kaiss et al. (in preparation).

## Citing Lirion

If you find Lirion useful in your work, we kindly request that you cite it:

```bibtex
@software{kaiss2026lirion,
  author  = {Kaiss, Mateus and Küster, Kristie Kaminski and Gebbran, Daniel},
  title   = {{Lirion.jl}: A Unified Framework for Centralized and Certified Distributed Optimal Power Flow},
  year    = {2026},
  version = {0.1.0},
  url     = {https://github.com/MateusKaiss/Lirion.jl}
}
```

A white paper describing Lirion is in preparation; this entry will be updated with the publication details once it is available.

## Acknowledgments

Lirion.jl is the result of a collective research effort under the coordination of the Federal University of Parana. The project was conceptualized by Mateus Kaiss and Daniel Gebbran, evolving towards a comprehensive open-source framework. The key collaborators of this work are:

- **Mateus Kaiss**: Researcher and Lead Developer. Responsible for R&D, implementation, maintenance, and technical architecture of the project.

- **Daniel Gebbran, PhD**: Research Lead. Established the scientific workflow & logic, and provides continuous high-level technical direction, including system architecture design, collaborative debugging, and feature implementation oversight.

- **Prof. Dr. Kristie Kaminski Küster**: Academic Advisor. Provides critical guidance in power systems engineering, reviews the project's logic and implementation, and coordinates the research group's scientific strategy.


## Julia community

Lirion was developed with deep respect for the Julia open-source optimization ecosystem. In particular, **PowerModels.jl** has been an important source of inspiration throughout the developer’s academic training, ever since the start of his undergraduate studies.

Lirion does not aim to replicate PowerModels, but rather to explore complementary research directions, especially in distributed optimization and new methods for OPF in electrical power systems — driven by the creator’s strong interest in applied mathematics, optimization theory, and scientific programming.

**PowerModels.jl**:
https://github.com/lanl-ansi/PowerModels.jl

## Funding

The authors gratefully acknowledge Brazil's Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq) for its support through grant 447271/2024-5 — established in 2024 and executed from 2025 through present — led by Daniel Gebbran as the Principal Investigator. This support has been essential in providing the institutional framework, scientific scope, and computational resources essential for the conceptualization and development of Lirion.jl, enabling the continued development, scientific expansion, and open-source release of this work.

## License

Lirion.jl is released under the MIT License.

See the [LICENSE](LICENSE) file for details.
