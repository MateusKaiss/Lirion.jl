<p align="center">
  <img src="docs/src/assets/logo.png" alt="Lirion.jl" 
  width="600"/>
</p>

<h3 align="center">
  Large-scale Integrated Research on Innovative Optimization of Networks

</h3>

A Julia package for solving Optimal Power Flow (OPF) problems using both centralized and distributed algorithms.  

> **Current focus:** Advanced ALADIN implementations and structured experimentation with decomposition and coupling strategies. Additional distributed methods are under active development.


## Why Lirion?

Lirion was built to bridge implementation clarity and theoretical rigor in Optimal Power Flow research.

It provides:

- Centralized DC, AC, and SOCP OPF
- Distributed OPF via ADMM and ALADIN
- Modular coupling-variable experimentation
- Flexible partitioning strategies
- Structured diagnostics for research

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

Below are two examples, one **centralized** and one **distributed (ADMM)**, to get started with Lirion: 

The example cases below are based on the well-known **PGLib-OPF** benchmark library and are included in the repository for convenience.

If no solver is specified, Lirion automatically selects its default solver:

- **DC → OSQP**
- **AC → Ipopt**
- **SOCP → Clarabel**

For detailed installation instructions, solver selection guidance, and a complete list of supported solvers, see the documentation:

**Installation Guide**  
https://your-docs-link-here.com/installation

---

### 🔹 Centralized Example (AC / DC / SOCP)

**Case file:**

[`test/matpower_cases/pglib_opf_case118_ieee.m`](https://github.com/MateusKaiss/Lirion.jl/blob/master/test/matpower_cases/pglib_opf_case118_ieee.m)


#### Runing centralized AC OPF

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

You may also switch the formulation:

```julia
model = dc()     # DC OPF  → typically uses OSQP.Optimizer
model = socp()   # SOCP OPF → typically uses Clarabel.Optimizer
```

### 🔹 Distributed Example ( ADMM AC / ADMM DC / ADMM SOCP)
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

## Mathematical Foundations Website

Lirion is accompanied by a dedicated technical website jointly created by Mateus Kaiss, creator of Lirion.jl, and Andressa Mathias Pedroso, Founder of Creatum.

The website contains deeper mathematical material underlying Lirion’s design, including:

- OPF formulations (DC, AC, SOCP)  
- ADMM derivations and penalty update rules  
- Coupling-variable design theory  
- Structural and decomposition analysis  

The material originates from structured research developed during undergraduate studies and is currently being expanded as part of ongoing graduate research.

https://your-website-link-here.com


## Acknowledgments

Lirion.jl originated as part of the undergraduate thesis of its principal developer, in the field of power systems. The package was conceived, designed, and developed by **Mateus Kaiss**, as its principal author and maintainer, with support from:

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
