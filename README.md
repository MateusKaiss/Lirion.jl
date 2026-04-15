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

## Documentation

The full documentation is available online and includes:

- Installation guide  
- Modeling workflows (centralized and distributed)  
- Solver configuration  
- Examples and API reference  

https://your-docs-link-here.com

---

## Acknowledgments

Lirion.jl originated from the undergraduate thesis of **Mateus Kaiss**, its principal developer, and has since evolved into a CNPq-funded research project. The key collaborators of this work are:

- **Dr. Daniel Gebbran Cons Bacilla Ferreira**, leader of the CNPq-funded research project, for his supervision of the research involving distributed methodologies and for his significant role in strengthening the scientific contributions of this work.

- **Prof. Dr. Kristie Kaminski Küster**, academic advisor of both the undergraduate thesis and the Master's research that emerged from this project, for her guidance in power systems research and for reviewing the implemented code.

- **Andressa Mathias Pedroso**, principal responsible for the design and development of the official Lirion website.

## Funding

The authors gratefully acknowledge CNPq for its support, which has been essential in enabling the continued development, scientific expansion, and open-source release of this work.

Lirion was developed with deep respect for the Julia open-source optimization ecosystem. In particular, **PowerModels.jl** has been an important source of inspiration throughout all the developer’s academic formation.

Lirion does not aim to replicate PowerModels, but rather to explore complementary research directions, especially in distributed optimization and new methods for OPF in electrical power systems — driven by the creator’s strong interest in applied mathematics, optimization theory, and scientific programming.

**PowerModels.jl**:
https://github.com/lanl-ansi/PowerModels.jl


## Mathematical Foundations Website

Lirion is accompanied by a dedicated technical website jointly created by Mateus Kaiss and Andressa Mathias Pedroso.

The website contains deeper mathematical material underlying Lirion’s design, including:

- OPF formulations (DC, AC, ACR, IVR, LPAC, SOCWR, SOCBF, BFA, NFA)
- ADMM and ALADIN derivations and penalty update rules
- Coupling-variable design theory
- Structural and decomposition analysis

The material originates from structured research developed during undergraduate studies and is currently being expanded as part of ongoing graduate research.

https://your-website-link-here.com

## License

Lirion.jl is released under the MIT License.

See the [LICENSE](LICENSE) file for details.
