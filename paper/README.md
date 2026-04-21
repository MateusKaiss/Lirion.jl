# Paper reproduction

This folder is **not** part of the Lirion package — `Pkg.add("Lirion")`
does not pull it in. It is bundled with the repository only so that
reviewers and researchers can reproduce the paper's centralized vs.
distributed OPF benchmarks locally.

## Contents

- `paper_results.jl` — main reproduction script
- `cases/<case>/std/*.m` — PGLib-OPF std-variant cases, with pre-computed
  Metis partitions (`*_tie.m` for TieLine / DummyBus, `*_splitbus.m`
  for SplitBus with the bilateral-zone repair)

## Cases and partitions

| Case | Buses | Zones | Tie-lines | Tie-lines (SplitBus) |
|---|---:|---:|---:|---:|
| case14_ieee | 14 | 2 | 3 | 3 |
| case30_ieee | 30 | 2 | 6 | 6 |
| case57_ieee | 57 | 2 | 8 | 8 |

## Usage

```julia
using Lirion
# Optional but recommended: also have one of Gurobi / Ipopt / Clarabel / OSQP loaded.
include(joinpath(pkgdir(Lirion), "paper", "paper_results.jl"))
```

The script:

1. Solves the centralized DC / AC / SOCWR OPF on each case as reference.
2. Runs the distributed ADMM across every combination of
   `{state, power, hybrid}` coupling × `{tie_line, dummy_bus, split_bus}`
   sectioning (18 combinations per case).
3. Prints the per-run gap, iteration count, and wall time, plus a
   summary table at the end.

## ρ and `max_iter` table

The reproduction uses the per-`(case, base_model, coupling)` ρ values
established by the paper's tuning study (see `paper_results.jl`:
`CASE_OVERRIDES`). In particular:

- `state` coupling on any base needs ρ = 1e5 on `case14_ieee` and
  ρ = 1e6 on `case30_ieee` / `case57_ieee` (smaller ρ leaves the
  primal consensus unclosed for tens of thousands of iterations).
- DC `power`-only on `case30_ieee` needs ρ = 1e8 with `max_iter = 10000`
  to escape the known ~1.2 % angle-drift plateau. AC/SOCWR `power`
  cannot absorb that ρ (solver ill-conditioning) and stay at the
  default ρ = 1e3.
- `hybrid` converges at the default ρ = 1e3 across all cases.

Every other combination uses `ρ = {state: 1e5, power: 1e3, hybrid: 1e3}`
with `max_iter = 5000` and `eps_abs = eps_rel = 1e-4`.
