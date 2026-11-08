# Third-Party Notices

Lirion.jl is distributed under the MIT License (see [`LICENSE`](LICENSE)).
This document lists third-party material that is either redistributed
within this repository or that Lirion.jl is partially derived from, along
with the corresponding attributions and license texts required by the
upstream projects.

---

## 1. InfrastructureModels.jl — adapted source code

- **Project:** InfrastructureModels.jl
- **Upstream:** https://github.com/lanl-ansi/InfrastructureModels.jl
- **Copyright:** (c) 2016, Triad National Security, LLC. All rights reserved.
- **License:** BSD-3-clause-style (LANL/Triad)
- **Where used in Lirion:** [`src/parsing/matlab_parser_core.jl`](src/parsing/matlab_parser_core.jl)
  contains routines adapted from the InfrastructureModels MATLAB parser.
  The file's header explicitly identifies it as a modified/adapted work.
- **Full license text:** [`licenses/InfrastructureModels_LICENSE.txt`](licenses/InfrastructureModels_LICENSE.txt)

---

## 2. PGLib-OPF — benchmark data

- **Project:** IEEE PES Power Grid Library — Optimal Power Flow (PGLib-OPF)
- **Upstream:** https://github.com/power-grid-lib/pglib-opf
- **Data License:** Creative Commons Attribution 4.0 International (CC-BY-4.0)
- **Software License:** MIT
- **Where used in Lirion:** the MATPOWER case files under
  [`test/matpower_cases/`](test/matpower_cases/) whose names start with
  `pglib_opf_*` are taken from PGLib-OPF and are used solely as
  benchmark/test inputs. Each case file preserves its original copyright
  header (Creative Commons attribution and original IEEE/PSAP/PEGASE/etc.
  source credits) inline.

---

## 3. PowerModels.jl — test-only dependency and inspiration

- **Project:** PowerModels.jl
- **Upstream:** https://github.com/lanl-ansi/PowerModels.jl
- **Copyright:** (c) 2016, Los Alamos National Security, LLC. All rights reserved.
- **License:** BSD-3-clause-style (LANL)
- **Relationship to Lirion:** PowerModels.jl is **not** redistributed with
  Lirion.jl. It appears only in `[extras]` / `[targets].test` of
  `Project.toml` and is fetched by `Pkg` exclusively when the user runs
  `Pkg.test("Lirion")` to compare results against PowerModels. No
  PowerModels source code is bundled in this repository. The
  acknowledgments section of the [README](README.md#acknowledgments)
  credits PowerModels.jl as a long-standing source of inspiration during
  the principal author's academic formation.

---

## 4. Other runtime dependencies

The remaining dependencies declared in [`Project.toml`](Project.toml)
(JuMP, Clarabel, DataFrames, Ipopt, MathOptInterface, OSQP, etc.) are
fetched by the Julia package manager from the General registry at install
time. Their source code is **not** redistributed inside this repository,
so their licenses do not need to be reproduced here. End users who install
Lirion.jl receive each of those packages directly from its upstream source
under its own license.
