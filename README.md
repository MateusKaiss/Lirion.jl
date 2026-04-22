<p align="center">
  <img src="docs/src/assets/logo.png" alt="Lirion.jl" 
  width="600"/>
</p>

<h3 align="center">
  Large-scale Integrated Research on Innovative Optimization of Networks

</h3>

<p align="center">
  <b>Created and maintained by <a href="https://github.com/MateusKaiss">Mateus Kaiss</a></b><br/>
</p>

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
- **Distributed OPF** via ADMM (DC / AC / SOCWR) and ALADIN
- Modular **coupling-variable** experimentation: `:state`, `:power`, `:hybrid`
- Flexible **partitioning/sectioning** strategies: `:tie_line`, `:dummy_bus`, `:metis`
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

Below are examples showing the **centralized**, **distributed ADMM**, and **distributed ALADIN** workflows in Lirion.

The example cases below are based on the well-known **PGLib-OPF** benchmark library and are included in the repository for convenience.

If no solver is specified, Lirion automatically selects its default solver per formulation:

- **LP/QP formulations (DC, BFA, NFA)** → OSQP
- **NLP formulations (AC, ACR, IVR, LPAC)** → Ipopt
- **SOCP formulations (SOCWR, SOCBF)** → Clarabel

For detailed installation instructions, solver selection guidance, and a complete list of supported solvers, see the documentation:

**Installation Guide**
https://your-docs-link-here.com/installation

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

---

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

## Documentation

The full documentation is available online and includes:

- Installation guide  
- Modeling workflows (centralized and distributed)  
- Solver configuration  
- Examples and API reference  

https://your-docs-link-here.com

---

## Mathematical Foundations Website

Lirion is accompanied by a dedicated technical website jointly created by Mateus Kaiss, creator of Lirion.jl, and Andressa Mathias Pedroso, Founder of Creatum.

The website contains deeper mathematical material underlying Lirion’s design, including:

- OPF formulations (DC, AC, ACR, IVR, LPAC, SOCWR, SOCBF, BFA, NFA)
- ADMM and ALADIN derivations and penalty update rules
- Coupling-variable design theory
- Structural and decomposition analysis

The material originates from structured research developed during undergraduate studies and is currently being expanded as part of ongoing graduate research.

https://your-website-link-here.com


## LirionLab.jl — Experimental Research Companion

**Lirion.jl is intentionally kept stable, minimal, and tightly scoped.** Its core formulations, API, and architecture are maintained by the principal author and evolve under a conservative roadmap to preserve reproducibility and scientific integrity.

To support ongoing research, student projects, and exploratory methodologies without compromising the stability of the core package, a separate companion repository is planned:

**LirionLab.jl** — an open research playground where collaborators, students, and researchers can prototype new decomposition schemes, coupling strategies, solver integrations, and OPF formulations on top of Lirion.jl.

The intended relationship is:

- **Lirion.jl** — stable core; curated formulations and algorithms; minimal external contributions; strict API stability.
- **LirionLab.jl** — experimental extensions; depends on Lirion.jl as a library; open to broader contributions; no API stability guarantees.

This separation ensures that Lirion.jl remains a dependable scientific tool, while LirionLab.jl provides a structured space for methodological exploration and academic collaboration.


## Acknowledgments

Lirion.jl originated from the principal developer’s undergraduate thesis in power systems and has since evolved into an ongoing research line in his Master's program (since March 2026). The package was conceived, designed, and developed by **Mateus Kaiss**, as its principal author and maintainer, with support from:

- **Dr. Daniel Gebbran Cons Bacilla Ferreira** for his supervision of the research involving distributed methodologies and for his significant role in strengthening the scientific contributions of this work.

- **Prof. Dr. Kristie Kaminski Küster** for her supervision and guidance in power systems research, as well as for reviewing the implemented code.

- **Andressa Mathias Pedroso** for leading and implementing the design and development of the official Lirion website.


Lirion was developed with deep respect for the Julia open-source optimization ecosystem. In particular, **PowerModels.jl** has been an important source of inspiration throughout all the developer’s academic formation.

Lirion does not aim to replicate PowerModels, but rather to explore complementary research directions, especially in distributed optimization and new methods for OPF in electrical power systems — driven by the creator’s strong interest in applied mathematics, optimization theory, and scientific programming.

**PowerModels.jl**:
https://github.com/lanl-ansi/PowerModels.jl


## License

Lirion.jl is released under the MIT License.

See the [LICENSE](LICENSE) file for details.
