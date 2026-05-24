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

## Why Lirion?  

Lirion was built to bridge implementation clarity and theoretical rigor in Optimal Power Flow research.

It provides:

- **Eleven centralized OPF formulations**:
  - **DC** — linear active-power approximation (LP)
  - **AC** — full nonlinear AC in polar coordinates (NLP)
  - **ACR** — full nonlinear AC in rectangular coordinates (NLP)
  - **IVR** — current-voltage rectangular formulation (NLP)
  - **LPAC** — linear-programming AC approximation (QCQP)
  - **SOCWR** — W-rectangular second-order cone relaxation (SOCP)
  - **QC** — quadratic-convex relaxation (QCQP)
  - **SDPWRM** — semidefinite W-matrix relaxation (SDP)
  - **SOCBF** — branch-flow second-order cone relaxation (SOCP)
  - **BFA** — linearized branch-flow approximation / LinDistFlow (LP)
  - **NFA** — network-flow approximation, active-power only (LP)
- **Distributed OPF** — two decomposition algorithms over DC / AC / SOCWR, each with its own variants:
  - **ADMM** — alternating direction method of multipliers
    - fixed or **adaptive penalty (ρ)** with residual balancing
  - **ALADIN** — augmented Lagrangian based alternating direction inexact Newton
    - **proximal mode** (`:primal` / `:pflow`) and **stopping criterion** (`:boyd` / `:hfd`)
    - optional **adaptive ρ / μ** growth
  - **Coupling variables** — which boundary quantities neighboring zones exchange:
    - **`:state`** — voltage state only: θ (DC), θ + V (AC), or v², R, T (SOCWR)
    - **`:power`** — branch active/reactive flows only
    - **`:hybrid`** — both state and power.
  - **Network sectioning** — how the grid is split into zones:
    - **`:tie_line`** — zones overlap on shared tie-line buses
    - **`:dummy_bus`** — each tie-line is split by an intermediate dummy bus
    - **`:split_bus`** — boundary buses are cloned across zones with injection consensus

> **Supported combinations.**
> - **`:state`** — works only with `:tie_line`
> - **`:power`** — works with `:tie_line` or `:dummy_bus`
> - **`:hybrid`** — works with all three (`:tie_line`, `:dummy_bus`, `:split_bus`), except under **ALADIN + SOCWR**, which does not support `:dummy_bus`

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

Below are examples showing the **centralized** and **distributed** (ADMM / ALADIN) workflows in Lirion.

The example cases below are based on the **PGLib-OPF** benchmark library and are included in the repository for convenience.

If no solver is specified, Lirion automatically selects its default solver per formulation:

- **OSQP** — DC, NFA
- **Ipopt** — AC, ACR, IVR, LPAC, QC
- **Clarabel** — SOCWR, SOCBF, SDPWRM, BFA

For detailed installation instruction and solver selection guidance, please refer to documentation of each solver and to the basic modelling of each problem herein formulated. 

---

### 🔹 Centralized Example

**Case file:**

[`test/matpower_cases/pglib_opf_case118_ieee.m`](https://github.com/MateusKaiss/Lirion.jl/blob/master/test/matpower_cases/pglib_opf_case118_ieee.m)


#### Running centralized AC OPF

The first argument is the path to a MATPOWER case file on your local machine.

```julia
using Lirion

out = solve("pglib_opf_case118_ieee.m";
    model = ac(),
    algorithm = centralized(),
    outputs = false,
    plots = false,
)

out.objective
```

You may switch to any of the eleven available formulations by changing the `model` argument:

```julia
model = dc()      # DC linear approximation       (LP)
model = ac()      # AC polar (V, θ)               (NLP)
model = acr()     # AC rectangular (vr, vi)       (NLP)
model = ivr()     # Current-voltage rectangular   (NLP)
model = lpac()    # Linear-programming AC         (QCQP)
model = socwr()   # W-rectangular SOCP            (SOCP)
model = qc()      # Quadratic-convex relaxation   (QCQP)
model = sdpwrm()  # Semidefinite W-matrix         (SDP)
model = socbf()   # Branch-flow SOCP              (SOCP)
model = bfa()     # Linearized branch flow        (LP)
model = nfa()     # Network-flow approximation    (LP)
```

### 🔹 Distributed Examples

A distributed run uses the same `solve(...)` entry point; you choose the
**algorithm** (`admm` / `aladin`), the **coupling**, and the **sectioning**.

**Case file:**
[`test/matpower_cases/pglib_opf_case73_ieee_rts.m`](https://github.com/MateusKaiss/Lirion.jl/blob/master/test/matpower_cases/pglib_opf_case73_ieee_rts.m)

#### Run a distributed AC OPF with ADMM (hybrid coupling)

```julia
using Lirion

out = solve("pglib_opf_case73_ieee_rts.m";
    model = ac(),
    algorithm = admm(
        coupling = :hybrid,
        sectioning = :tie_line,
        max_iter = 350,

        # A single rho automatically fills the required penalty parameters
        rho = 1e3,

        eps_abs = 1e-7,
        eps_rel = 1e-7,

        # Adaptive penalty update
        adaptive_rho = true,
        mu_incr = 3.0,
        mu_decr = 3.0,
        tau_incr = 1.05,
        tau_decr = 0.98,
    ),
    outputs = false,
    plots = true,
)

out.diagnostics
```

ADMM supports **adaptive penalty (ρ) updates** following residual-balancing rules
(`mu_incr`, `mu_decr`, `tau_incr`, `tau_decr`). Each coupling has its own required ρ
parameters: a single `rho=...` fills them all, or set them individually (`rho_theta`,
`rho_v`, `rho_p`, `rho_q`, `rho_v2`, `rho_r`, `rho_t`), depending on the formulation.

#### Run a distributed AC OPF with ALADIN (state coupling)

```julia
using Lirion

out = solve("pglib_opf_case73_ieee_rts.m";
    model = ac(),
    algorithm = aladin(
        coupling = :state,
        sectioning = :tie_line,
        max_iter = 100,

        # Augmented-Lagrangian penalties
        rho = 1e3,   # local proximal penalty (subproblems)
        mu  = 1e4,   # coordinator consensus penalty (QP)

        eps_abs = 1e-4,
        eps_rel = 1e-4,

        # Optional variants
        stopping_criterion = :boyd,   # :boyd (2-norm) or :hfd (∞-norm)
        adaptive = false,             # when true, grows rho / mu across iterations
    ),
    outputs = false,
    plots = false,
)

out.diagnostics
```

ALADIN solves local NLPs and reconciles the zones through a **coordinator QP** built
from first/second-order sensitivities. Beyond `rho`/`mu` it exposes a few variants:

- **`proximal_mode`** — `:primal` (default) or `:pflow` (`:power`/`:hybrid`, DC/SOCWR only)
- **`stopping_criterion`** — `:boyd` (scaled 2-norm) or `:hfd` (∞-norm)
- **`adaptive`** — grow `rho` / `mu` via `rho_growth` / `mu_growth` when enabled

#### Switching the configuration

```julia
# Algorithm
algorithm = admm(...)        # ADMM
algorithm = aladin(...)      # ALADIN

# Coupling — which boundary quantities are exchanged
coupling = :state            # voltage state (θ / V / v²,R,T)
coupling = :power            # branch active/reactive flows
coupling = :hybrid           # both

# Network sectioning — how the grid is split into zones
sectioning = :tie_line       # zones share boundary buses
sectioning = :dummy_bus      # tie-lines split by dummy buses
sectioning = :split_bus      # boundary buses cloned across zones
```

#### Automatic partitioning

You can let Lirion assign the zones automatically with impedance-weighted **Metis**
graph partitioning (balanced zones, connectivity-repaired). This needs the optional
extension dependencies:

```julia
using Lirion, Metis, Graphs, SimpleWeightedGraphs

# Split a case into 4 zones → writes pglib_opf_case118_ieee_4zones.m
create_partition("pglib_opf_case118_ieee.m"; k = 4)
```

The resulting zoned case is then solved with any distributed algorithm/sectioning.
For `:split_bus`, use `create_partition_split_bus(...)`: unlike the others, it
requires every boundary bus to border **at most one** other zone (no bus on 3+
zones), a bilateral constraint the plain partitioner does not enforce.

## Results

`solve(...)` returns a `LirionResult` holding the solution as DataFrames (mirroring
the MATPOWER schema), plus summary fields:

```julia
out.objective   # optimal cost
out.status      # :optimal, :local_optimal, ...
out.bus         # per-bus voltages / angles
out.gen         # generator dispatch (pg, qg)
out.branch      # branch flows (pf, qf, pt, qt)
```

Distributed runs additionally expose the per-zone results and the boundary consensus:

```julia
out.zones         # per-zone (bus, gen, branch, cost)
out.consensus     # boundary-variable deviations across zones
out.diagnostics   # iteration history (residuals, ρ, ...)
```

Export everything to CSV (requires `using CSV, DataFrames`):

```julia
Lirion.to_csvs(out, "results/")
```

## Performance & Threading

Lirion has two very different workloads and each one has its own
optimal threading configuration. The golden rule is **never nest
parallelism**: either go *deep* (one solver, many BLAS threads) or
*wide* (many zones, one BLAS thread each).

### Centralized — leave the defaults

One solver instance; give BLAS all the physical cores so matrix
factorizations are fast. No environment variables needed — in particular,
do **not** set `OPENBLAS_NUM_THREADS=1` here, or factorizations will be slow.

### Distributed — pin BLAS and set Julia threads

The one setting that **must** be in place before launching Julia is
`JULIA_NUM_THREADS` (it cannot be changed at runtime). Lirion already scopes
BLAS to a single thread during the parallel zone solves, so pinning the BLAS
environment variables below is belt-and-suspenders — it just covers the phases
outside that scope.

**Linux (bash / zsh):**

```bash
export JULIA_NUM_THREADS=$(getconf _NPROCESSORS_ONLN)
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
julia --project
```

**macOS (bash / zsh):**

```bash
export JULIA_NUM_THREADS=$(sysctl -n hw.ncpu)
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
julia --project
```

**Windows (PowerShell):**

```powershell
$env:JULIA_NUM_THREADS   = $env:NUMBER_OF_PROCESSORS
$env:OPENBLAS_NUM_THREADS = "1"
$env:OMP_NUM_THREADS      = "1"
$env:MKL_NUM_THREADS      = "1"
julia --project
```

After loading, Lirion prints the active configuration:

```
[ Info: [Lirion] Threads: Julia=28 | BLAS=1
```

For the full rationale (why hyperthreads hurt BLAS, why Lirion
auto-scopes BLAS during distributed runs, verification tips), see
the **Performance & Threading** page in the documentation.

---

## Documentation

Full documentation — installation, centralized & distributed workflows, solver
configuration, examples, and API reference — is in preparation. The link will
appear here once it is published.

## Mathematical Foundations Website

Lirion is accompanied by a dedicated technical & instructive website, created by Andressa Mathias Pedroso and Mateus Kaiss, providing learning resources for anyone interested in expanding their knowledge of these topics.

Rather than just cataloguing the formulations, the site builds the theory **from
the ground up** — based on Mateus Kaiss's studies — deriving each model from first
principles, including:

- **Physical foundations** — from Maxwell's equations and circuit laws to the power-flow equations
- **Optimization theory** — convexity, KKT conditions, duality, and the relaxation hierarchy
- **The eleven OPF formulations** (DC, AC, ACR, IVR, LPAC, SOCWR, QC, SDPWRM, SOCBF, BFA, NFA) and how each is derived
- **Distributed methods** — ADMM and ALADIN derivations and penalty/update rules
- **Coupling-variable design** and network decomposition analysis

*Link coming soon.*


---

## Scientific Outputs & Deliverables

This software is part of an ongoing research program. The following works have been produced as part of this project's development:

1. **[White Paper]** - *Lirion.jl: A Distributed Framework for Power Systems*, Mateus Kaiss et al. (In Preparation).
2. *On the Impact of Network Sectioning and Coupling Variable Choices in ADMM-Based Distributed Optimal Power Flow*, Mateus Kaiss et al. (In Preparation).
3. **BSc Thesis** - *Desenvolvimento de Ferramenta para Determinação do Fluxo de Potência Ótimo em Redes Elétricas*, Mateus Kaiss, UFPR (2025). [Link to University Repository/PDF]

## Citing Lirion

If you find Lirion useful in your work, we kindly request that you cite it:

```bibtex
@software{kaiss2026lirion,
  author  = {Kaiss, Mateus and Gebbran, Daniel and Küster, Kristie Kaminski},
  title   = {{Lirion.jl}: A {Julia} Framework for Centralized and Distributed Optimal Power Flow},
  year    = {2026},
  version = {0.1.0},
  url     = {https://github.com/MateusKaiss/Lirion.jl}
}
```

A white paper describing Lirion is in preparation; this entry will be updated with the publication details once it is available.

## Acknowledgments

Lirion.jl is the result of a collective research effort under the coordination of the Federal University of Parana. The project was conceptualized by Mateus Kaiss and Daniel Gebbran, evolving towards a comprehensive open-source framework. The key collaborators of this work are:

- **Mateus Kaiss**: Researcher and Lead Developer. Responsible for R\&D, implementation, maintenance, and technical architecture of the project.

- **Daniel Gebbran, PhD**: Research Lead. Established the scientific workflow \& logic, and provides continuous high-level technical direction, including system architecture design, collaborative debugging, and feature implementation oversight. 

- **Prof. Dr. Kristie Kaminski Küster**: Academic Advisor. Provides critical guidance in power systems engineering, reviews the project's logic and implementation, and coordinates the research group's scientific strategy.

- **Andressa Mathias Pedroso**: Lead Web Developer. Responsible for the design and development of the official Lirion platform and online presence.

## Funding

The authors gratefully acknowledge Brazil's Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq) for its support through grant 447271/2024-5 — established in 2024 and executed from 2025 through present — led by Daniel Gebbran as the Principal Investigator. This support has been essential in providing the institutional framework, scientific scope, and computational resources essential for the conceptualization and development of Lirion.jl, enabling the continued development, scientific expansion, and open-source release of this work.

## Julia community

Lirion was developed with deep respect for the Julia open-source optimization ecosystem. In particular, **PowerModels.jl** has been an important source of inspiration throughout the developer’s academic training, ever since the start of his undergraduate studies.

Lirion does not aim to replicate PowerModels, but rather to explore complementary research directions, especially in distributed optimization and new methods for OPF in electrical power systems — driven by the creator’s strong interest in applied mathematics, optimization theory, and scientific programming.

**PowerModels.jl**:
https://github.com/lanl-ansi/PowerModels.jl

## License

Lirion.jl is released under the MIT License.

See the [LICENSE](LICENSE) file for details.
