<p align="center">
  <img src="docs/src/assets/logo.png" alt="Lirion.jl" 
  width="600"/>
</p>

<h3 align="center">
  Large-scale Integrated Research on Innovative Optimization of Networks

</h3>

A Julia package for solving Optimal Power Flow (OPF) problems using both centralized and distributed algorithms.

## Why Lirion?

Lirion was built to bridge implementation clarity and theoretical rigor in Optimal Power Flow research.

It provides:

- **Nine centralized OPF formulations**:
  - **DC** — linear active-power approximation (LP)
  - **AC** — full nonlinear AC in polar coordinates (NLP)
  - **ACR** — full nonlinear AC in rectangular coordinates (NLP)
  - **IVR** — current-voltage rectangular formulation (NLP)
  - **LPAC** — linear-programming AC approximation (Coffrin & Van Hentenryck, QCQP)
  - **SOCWR** — Jabr W-rectangular second-order cone relaxation (SOCP)
  - **SOCBF** — Farivar–Low branch-flow second-order cone relaxation (SOCP)
  - **BFA** — linearized branch-flow approximation / LinDistFlow (LP)
  - **NFA** — network-flow approximation, active-power only (LP)
- **Distributed OPF** via ADMM (DC / AC / SOCWR) (and upcoming, via ALADIN)
- Modular **coupling-variable** experimentation: `:state`, `:power`, `:hybrid`
- Flexible **partitioning** methods: `:tie_line`, `:dummy_bus`
- Smarter **sectioning** strategy, based on `:metis`
- **Adaptive penalty (ρ) updates** with per-block residual balancing
- Structured diagnostics and plots for research

---

## Installation

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

Below are examples showing the **centralized** and **distributed ADMM** workflows in Lirion.

The example cases below are based on the well-known **PGLib-OPF** benchmark library and are included in the repository for convenience.

If no solver is specified, Lirion automatically selects its default solver per formulation:

- **LP/QP formulations (DC, BFA, NFA)** → OSQP
- **NLP formulations (AC, ACR, IVR, LPAC)** → Ipopt
- **SOCP formulations (SOCWR, SOCBF)** → Clarabel

For detailed installation instruction and solver selection guidance, please refer to documentation of each solver and to the basic modelling of each problem herein formulated. 

(Upcoming) **Installation Guide** ...

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
    optimizer = Ipopt.Optimizer,
    outputs = false,
    plots = false,
)

out.result
```

You may switch to any of the nine available formulations by changing the `model` argument:

```julia
model = dc()      # DC linear approximation       (LP,   OSQP)
model = ac()      # AC polar (V, θ)               (NLP,  Ipopt)
model = acr()     # AC rectangular (vr, vi)       (NLP,  Ipopt)
model = ivr()     # Current-voltage rectangular   (NLP,  Ipopt)
model = lpac()    # Linear-programming AC         (QCQP, Ipopt)
model = socwr()   # Jabr W-rectangular SOCP       (SOCP, Clarabel)
model = socbf()   # Farivar–Low branch-flow SOCP  (SOCP, Clarabel)
model = bfa()     # Linearized branch flow        (LP,   OSQP)
model = nfa()     # Network-flow approximation    (LP,   OSQP)
```

### 🔹 Distributed Example (ADMM — DC / AC / SOCWR)
#### Run a Distributed AC OPF with ADMM (Hybrid Coupling)


**Case file:**

[`test/matpower_cases/pglib_opf_case73_ieee_rts.m`](https://github.com/MateusKaiss/Lirion.jl/blob/master/test/distributed/pglib_opf_case73_ieee_rts.m)


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

Lirion supports **adaptive penalty parameter updates** following residual balancing rules.
The parameters `mu_incr`, `mu_decr`, `tau_incr`, and `tau_decr` control how the ADMM penalty is adjusted based on primal and dual residuals.

The **coupling** argument selects which variables are exchanged between subproblems:

- `:state` — only state variables (θ, V, or W) are coupled
- `:power` — only branch active/reactive powers are coupled
- `:hybrid` — both state and power couplings (per-block ρ)

Each coupling has its own required ρ parameters. Passing a single `rho=...` keyword fills all required penalties at once, or they can be set individually (e.g., `rho_theta`, `rho_v`, `rho_p`, `rho_q`, `rho_v2`, `rho_r`, `rho_t`), depending on the formulation.

For a detailed explanation of the update rules and their theoretical motivation, see the documentation:

**ADMM Parameter Tuning Guide** (link to docs section)

## Performance & Threading

Lirion has two very different workloads and each one has its own
optimal threading configuration. The golden rule is **never nest
parallelism**: either go *deep* (one solver, many BLAS threads) or
*wide* (many zones, one BLAS thread each).

### Centralized — leave the defaults

One solver instance; give BLAS all the physical cores so matrix
factorizations are fast. No environment variables needed.

### Distributed — pin BLAS and set Julia threads

Set these **before** launching Julia.

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

## Mathematical Foundations Website

Lirion is accompanied by a dedicated technical \& instructive website, providing learning resources for any person who might be interested in expanding their knowledge in such topics.

The website contains deeper mathematical material underlying Lirion’s design, including:

- OPF formulations (DC, AC, ACR, IVR, LPAC, SOCWR, SOCBF, BFA, NFA)
- ADMM and ALADIN derivations and penalty update rules
- Coupling-variable design theory
- Structural and decomposition analysis

## Documentation

The full documentation is available online and includes:

- Installation guide  
- Modeling workflows (centralized and distributed)  
- Solver configuration  
- Examples and API reference  

The full documentation site is being prepared; the link will appear here once it is published.

For a detailed explanation of the update rules and their theoretical motivation, upcoming documentation will address it: **ADMM Parameter Tuning Guide** (Upcoming)

---

## Acknowledgments

Lirion.jl is the result of a collective research effort under the coordination of the Federal University of Parana. The project was conceptualized by Mateus Kaiss and Daniel Gebbran, evolving towards a comprehensive open-source framework. The key collaborators of this work are:

- **Mateus Kaiss**: Researcher and Lead Developer. Responsible for R\&D, implementation, maintenance, and technical architecture of the project.

- **Daniel Gebbran, PhD**: Research Lead. Established the scientific workflow \& logic, and provides continuous high-level technical direction, including system architecture design, collaborative debugging, and feature implementation oversight. 

- **Prof. Dr. Kristie Kaminski Küster**: Academic Advisor. Provides critical guidance in power systems engineering, reviews the project's logic and implementation, and coordinates the research group's scientific strategy.

- **Andressa Mathias Pedroso**: Lead Web Developer. Responsible for the design and development of the official Lirion platform and online presence.

## Scientific Outputs & Deliverables

This software is part of an ongoing research program. The following works have been produced as part of this project's development:

1. **[White Paper]** - *Lirion.jl: A Distributed Framework for Power Systems* (In Preparation).
2. *On the Impact of Network Sectioning and Coupling Variable Choices in ADMM-Based Distributed Optimal Power Flow* (In Preparation).
3. **BSc Thesis** - *Desenvolvimento de Ferramenta para Determinação do Fluxo de Potência Ótimo em Redes Elétricas*, Mateus Kaiss, UFPR (2025). [Link to University Repository/PDF]
4. **Official Documentation** — Coming soon.


## Funding

The authors gratefully acknowledge Brazil's Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq) for its support through grant 447271/2024-5 — established in 2024 and executed from 2025 through present — led by Daniel Gebbran as the Principal Investigator. This support has been essential in providing the institutional framework, scientific scope, and computational resources essential for the conceptualization and development of Lirion.jl, enabling the continued development, scientific expansion, and open-source release of this work.


## Julia community

Lirion was developed with deep respect for the Julia open-source optimization ecosystem. In particular, **PowerModels.jl** has been an important source of inspiration throughout all the developer’s academic formation.

Lirion does not aim to replicate PowerModels, but rather to explore complementary research directions, especially in distributed optimization and new methods for OPF in electrical power systems — driven by the creator’s strong interest in applied mathematics, optimization theory, and scientific programming.

**PowerModels.jl**:
https://github.com/lanl-ansi/PowerModels.jl


## License

Lirion.jl is released under the MIT License.

See the [LICENSE](LICENSE) file for details.
