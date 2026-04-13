# Quick Start

This page shows the shortest path from **zero** to a working OPF solve in Lirion:
1) run a **centralized** OPF (single area), then  
2) run a **distributed** OPF with **ADMM** (multi-area).

---

## 1) Install and load Lirion

```julia
julia> import Pkg
julia> Pkg.add("Lirion")
julia> using Lirion
```

## 2) Download Example Cases

Many of the most widely used and benchmarked OPF test cases are available in the **PGLib-OPF library**, maintained by the Power Grid Library project.

You can explore and download additional cases directly from:

👉 **PGLib-OPF Repository**  
https://github.com/power-grid-lib/pglib-opf

---

For this Quick Start, however, we recommend downloading two specific cases directly from Lirion’s test suite to ensure compatibility with the examples below.

---

### Centralized Example Case (AC / DC / SOCWR)

Case file:

`pglib_opf_case118_ieee.m`

Download the raw file:

https://raw.githubusercontent.com/MateusKaiss/Lirion.jl/master/test/matpower_cases/pglib_opf_case118_ieee.m

---

### Distributed Example Case (ADMM-ready)

Case file:

`pglib_opf_case73_ieee_rts.m`

Download the raw file:

https://raw.githubusercontent.com/MateusKaiss/Lirion.jl/master/test/distributed/pglib_opf_case73_ieee_rts.m

---

These two cases are sufficient to reproduce the centralized and distributed examples presented in this Quick Start.

3) Centralized OPF (one-liner)

Lirion supports centralized DC, AC, and SOCWR OPF via:

dc(), ac(), socwr()

algorithm = centralized() (default)

Centralized AC OPF
julia> out = solve("cases/pglib_opf_case118_ieee.m";
           model = ac(),
           algorithm = centralized(),
           outputs = false,
           plots = false,
       );

julia> out.result

Centralized DC OPF
julia> out = solve("cases/pglib_opf_case118_ieee.m";
           model = dc(),
           algorithm = centralized(),
       );


If you don’t pass an explicit solver, Lirion uses its default solver selection:

DC → OSQP

AC → Ipopt

SOCWR → Clarabel

4) Distributed OPF with ADMM (minimal configuration)

Distributed runs use:

algorithm = admm(...)

a coupling choice: :state, :power, or :hybrid

a sectioning strategy (currently :tie_line)

ADMM (recommended starter): AC + hybrid coupling
julia> out = solve("cases/pglib_opf_case73_ieee_rts.m";
           model = ac(),
           algorithm = admm(
               coupling = :hybrid,
               sectioning = :tie_line,
               max_iter = 350,
               # a single rho fills the required rhos for the chosen model/coupling
               rho = 1e3,
               eps_abs = 1e-7,
               eps_rel = 1e-7,
               adaptive_rho = true,
           ),
           outputs = false,
           plots = true,
       );

julia> out.diagnostics

5) What to change next (the 4 knobs)

Once the examples above run, most experimentation in Lirion is just changing:

(1) OPF formulation

model = dc() / ac() / socwr()

(2) Centralized vs Distributed

algorithm = centralized()

algorithm = admm(...)

(3) Coupling choice (ADMM)

coupling = :state (boundary state consensus)

coupling = :power (tie-line power/flow consensus)

coupling = :hybrid (state + power)

(4) Penalty parameters (ρ)

You can pass either:

a single rho = ... (fills the required rhos for that configuration), or

the explicit rhos required by each model/coupling (advanced tuning).

Example (explicit tuning for AC + hybrid):

julia> alg = admm(
           coupling = :hybrid,
           sectioning = :tie_line,
           max_iter = 500,
           rho_theta = 1e3,
           rho_v     = 1e3,
           rho_p     = 1e2,
           rho_q     = 1e2,
           eps_abs   = 1e-7,
           eps_rel   = 1e-7,
           adaptive_rho = true,
       );
