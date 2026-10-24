# Reproducing Table I

This folder regenerates the distributed results table of the paper (ADMM and
ALADIN on case14, case24, case39, case73, case118 and case500, three
formulations, six coupling/sectioning rows each) from the exact parameters the
published numbers were obtained with, and checks the regenerated numbers
against the published ones.

One command per case, run from the repository root:

```sh
julia --project=. paper/reproduction/run.jl case14
```

It solves every cell of that case with the parameters in
`cases/case14/gold_params.toml`, builds the Tier-1 / Tier-2 crossing table from
the saved residual histories, and compares it with
`cases/case14/expected_crossing.csv`. The last line of the output reads

```
verify_crossing: 36 reproduce, 0 mismatch, 0 missing
```

when the case reproduces. Without a case argument every case is run; the
positional filters `case k formulation sectioning coupling` select a single
cell (`-` or absent = wildcard), for example
`julia --project=. paper/reproduction/run.jl case118 4 ac tie_line state`.

## What is in this folder

| path | what it is |
|---|---|
| `run.jl` | The driver. Reads the parameter file(s), solves each cell, writes one history CSV and one meta TOML per run, then calls the two scripts below. |
| `scripts/build_crossing_table.jl` | Post-hoc scorer. Reads the histories and metas, finds the Tier-1 and Tier-2 crossings and the verdict of every run, writes `results/table/crossing_table.csv` (and `.md`). |
| `scripts/verify_crossing.jl` | Compares a regenerated `crossing_table.csv` with one or more `expected_crossing.csv` files, row by row. |
| `scripts/repro_centralized_baseline.csv` | Centralized oracle objectives (one per case and formulation) used for the optimality gate. |
| `scripts/lib/` | Shared helpers used by the driver: paths, solver factories, CSV/TOML writers. |
| `cases/<case>/gold_params.toml` | The parameters of every cell of that case (ADMM and ALADIN blocks). Cells the paper prints as x carry `status = "non_convergent"` and are not solved. |
| `cases/<case>/expected_crossing.csv` | The published crossing table rows of that case: the file a reproduction has to match. |
| `cases/<case>/SCORECARD.md` | Human summary of those rows: Tier-1 (a), Tier-2 (b), gap and verdict per cell. |
| `expected_crossing_full.csv` | All 216 rows (six cases) in one file, identical to the per-case files concatenated. |
| `results/` | Created by `run.jl`: `results/runs/` (per-run histories and metas) and `results/table/` (the regenerated crossing table). Not versioned. |

## Prerequisites

- A checkout of this repository with its Julia environment instantiated:
  `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
- The partitioned cases under `paper/cases/partitioned/<case>/std/k<k>/case.m`
  (case14 k2, case24 k3, case39 k3, case73 k3, case118 k4, case500 k20). They
  ship with the repository.
- **HSL MA57.** The AC and SOCWRNL zone solvers run Ipopt with the MA57 linear
  solver, loaded from a licensed `libhsl.so`. Point `LIRION_HSLLIB` at the full
  HSL build (about 15.5 MB). A small placeholder `libhsl.so` of a few kilobytes
  will load but makes Ipopt fail with `Restoration_Failed`; do not rely on
  auto-detection, set the path explicitly.
- **Gurobi** with a valid license, for the DC zone subproblems and the DC
  centralized oracle. Gurobi is pinned to two threads in `scripts/lib/tuning_common.jl`
  (the case500 DC cells are thread-sensitive).
- Clarabel (ALADIN coordinator QP, all formulations) and Ipopt come with the
  Julia environment.

Environment, exported before running:

```sh
# required
export LIRION_HSLLIB=/path/to/libhsl.so      # full HSL build with MA57
export LIRION_BOUND_DETECT=primal_dual       # bound-activity detector used by every DC/AC ALADIN cell
export LIRION_STRICT_SOLVERS=1               # a solver fallback is an error, never silent
# recommended
export JULIA_NUM_THREADS=8                   # BLAS is pinned to 1 inside the driver
```

The driver sets `LIRION_BOUND_DETECT=primal_dual` and `LIRION_STRICT_SOLVERS=1`
itself if they are missing; exporting them keeps the provenance explicit.

Variables that must **not** be set (the driver aborts if they are, because
each one changes the trajectory away from the published configuration):
`LIRION_GAUGEFIX`, `LIRION_ZONE_WARMSTART`, `LIRION_ZONE_TOL_TIGHT`,
`LIRION_AC_PSD_MODE`. `LIRION_SOCWRNL_LIFTED_FLOW` is managed per cell by the
driver (set only for the three case500 SOCWRNL-ADMM cells whose block carries
`lifted_flow = true`); do not export it.

Optional overrides: `LIRION_REPRO_OUT` (output directory, default
`paper/reproduction/results/runs`), `LIRION_REPRO_GOLD` (an explicit parameter
TOML instead of `cases/<case>/gold_params.toml`), `LIRION_REPRO_NO_TABLE=1`
(skip the table and the comparison), `LIRION_REPRO_SMOKE=1` (allow a run
without MA57, for engineering checks only; the numbers will not match).

## What "reproduces" means

A row of the crossing table reproduces when these four columns equal the
values in `expected_crossing.csv`:

- `tier1_first_touch` — the paper's **a**, the first iteration at which the
  consensus mismatch is below 0.01 (the mismatch-only stop of the literature);
- `tier2_iter` — the paper's **b**, the certified iteration (primal residual
  and stationarity or dual residual under their bars, objective within 1% of
  the centralized solution);
- `tier2_first_touch_raw` — the raw Tier-2 first touch, which is the run's own
  stopping iteration;
- `verdict` — `pass`, or `non_convergent` for the cells the paper prints as x.

`scripts/verify_crossing.jl` does exactly this check and prints the objective
gap next to each row for information. Under a pinned toolchain (same Gurobi,
Clarabel and Ipopt versions, the full MA57, BLAS on one thread) the gaps are
bit-identical; otherwise they agree to solver tolerance. The paper prints
iteration counts of 10000 and above truncated to thousands (14763 appears as
14k); the CSV files keep the exact integers.

Cells that never certified during the reproduction have no parameter set to run.
The driver writes a placeholder meta for them (no solve), the crossing table
lists them as `non_convergent` with the reproduction note in `verdict_reason`, and
the comparison treats them like any other row. Their Tier-1 value in the paper
comes from the reproduction run, not from this driver.

A mismatch on a certified cell almost always points at the environment rather
than at the method: an AC toggle not applied, a wrong `bound_detect`, a leaked
or missing `LIRION_SOCWRNL_LIFTED_FLOW`, or the placeholder `libhsl.so`. Every
run's meta TOML records the resolved solvers and every toggle that was active
(`manifest` and `run_manifest`), which is where to look first.

## Running time

Wall time of the published run, single machine, one Julia process:

| case | zones | all 36 rows | longest single cell |
|---|---|---|---|
| case14 | 2 | about 7 min | 2 min (AC state/tie ALADIN) |
| case24 | 3 | about 9 min | 6 min (SOCWRNL power/tie ALADIN) |
| case39 | 3 | about 27 min | 19 min (AC hybrid/dummy ALADIN) |
| case73 | 3 | about 34 min | 14 min (AC power/tie ALADIN) |
| case118 | 4 | about 38 min | 9 min (AC power/tie ALADIN) |
| case500 | 20 | about 6 h | 2 h 45 min (SOCWRNL hybrid/dummy ALADIN) |

Run one heavy Julia process at a time; the AC and SOCWRNL zone solves are
compute-bound. For the full set use a detached session (`nohup`), and on a
thermally constrained machine pin the process (`taskset`,
`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`).

At case500 the zone solvers are sensitive to thread count and to floating-point
order, so an iteration count there may move by a few units between machines.
The certificates themselves (which cells pass) are stable; the paper prints
those counts truncated to thousands.

## Output files

Per run, in `results/runs/`:

- `<case>_k<k>_<form>_<sect>_<coup>_<algo>.csv` — one row per iteration:
  `iter, r_primal, eps_primal, r_dual, eps_dual, r_step_primal,
  eps_step_primal, lambda_step, sd, psd_delta`.
- `<case>_k<k>_<form>_<sect>_<coup>_<algo>.toml` — objective, iterations,
  runtime, the dual scaling `s_d`, `n_prox`, `atlambda_l1`, the full runtime
  manifest and the git commit of the source tree.

The scorer can be rerun on its own, and the comparison against any subset of
cases, without solving again:

```sh
julia --project=. paper/reproduction/scripts/build_crossing_table.jl \
    paper/reproduction/results/runs \
    paper/reproduction/scripts/repro_centralized_baseline.csv \
    paper/reproduction/results/table
julia paper/reproduction/scripts/verify_crossing.jl \
    paper/reproduction/results/table/crossing_table.csv \
    paper/reproduction/cases/case14/expected_crossing.csv
```
