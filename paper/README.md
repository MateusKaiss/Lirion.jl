# Lirion distributed-OPF paper — experiments

This folder is **not** part of the Lirion package (`Pkg.add("Lirion")`
does not pull it in). It bundles the reproducible experiment pipeline so
reviewers can rebuild the cases, partitions, tuning and tables locally.

## Design: partition vs. sectioning

In current Lirion, **sectioning is a runtime choice**, not a file. The
solver applies `TieLine` / `DummyBus` / `SplitBus` in memory from the
`bus.zone` column of the case:

```julia
solve(case; model = ac, algorithm = admm(sectioning = :dummy_bus, coupling = :hybrid))
# -> build_distributed_system(case, algorithm.sectioning)
#    -> prepare_distributed_partitions(case, sectioning)  # reads case.bus.zone
```

Therefore a partition is **one `case.m` carrying only the `zone`
column**, one per `(case, variant, k)`. There are **no** `_tie.m` /
`_dummybus.m` / `_splitbus.m` files — the sectionings to test are
recorded as metadata (`sectionings_planned`) in each `meta.toml`.

## PGLib-OPF data

Cases come from [PGLib-OPF](https://github.com/power-grid-lib/pglib-opf).
Root resolution order:

1. environment variable `PGLIB_OPF_ROOT`;
2. fallback `/home/mateus-kaiss/Downloads/papers/pglib-opf-master`.

Variants: `std` at `<root>/pglib_opf_<name>.m`, `api` at
`<root>/api/pglib_opf_<name>__api.m`, `sad` at
`<root>/sad/pglib_opf_<name>__sad.m`. The scripts **copy** (never
move/alter) the selected cases into `cases/raw/pglib/{std,api,sad}/`.

## Partitioning deps

`02_build_partitions.jl` needs the Lirion partition extension:

```julia
using Pkg; Pkg.add(["Metis", "Graphs", "SimpleWeightedGraphs"])
```

## Layout

```
cases/selected_cases.toml          # declarative case list
cases/raw/pglib/{std,api,sad}/     # copied PGLib .m + PGLIB_LICENSE + PGLIB_README
cases/partitioned/<case>/<variant>/k<k>/{case.m, meta.toml}
cases/partitioned/{_summary,_diagnostic}.toml
configs/partitioning.toml          # per-case k, seeds, refine
configs/selection_criteria.toml    # partition + tuning selection rules
configs/tuning_grid_{admm,aladin}.toml  # placeholders (not run yet)
scripts/lib/common.jl              # shared helpers
scripts/prepare/{01_select_cases,02_build_partitions}.jl
scripts/experiments/               # (future) solver / tuning drivers
scripts/tables/                    # (future) table + figure generation
results/logs/                      # per-script logs (git-ignored)
```

## Reproduce / audit (no solver)

```julia
using Lirion
using Metis, Graphs, SimpleWeightedGraphs
include(joinpath(pkgdir(Lirion), "paper", "scripts", "prepare", "01_select_cases.jl"))
include(joinpath(pkgdir(Lirion), "paper", "scripts", "prepare", "02_build_partitions.jl"))
```

Normal use reads the ready-made cases in `cases/partitioned/`; the
`prepare/` scripts are only the reconstruction/audit recipe.

## Partition status & main matrix

Each `meta.toml` carries `status`, `clean_for_main`, and
`optional_or_sensitivity`:

- `clean` (`clean_for_main = true`) — every zone connected **and** has a
  dispatchable generator. **The main matrix of the paper uses ONLY these.**
- `connected_no_gen` (`clean_for_main = false`,
  `optional_or_sensitivity = true`) — the partition is connected, but some
  zone has no dispatchable generator of its own and imports power through
  the tie-lines. This is a data property of the case (few dispatchable
  generators, e.g. case14/case30), not a partitioning defect. Kept as
  **optional / sensitivity** cases, out of the headline matrix.
- `not_clean` — some zone is not connected on the electrical graph.
  After the component-based connectivity repair
  (`Lirion.enforce_zone_connectivity!`) this count is **0**; a partition
  can only be `not_clean` if the underlying graph itself has more
  electrical islands than zones can absorb.

Aggregate `status_counts_by_case_k` are recorded in
`cases/partitioned/_summary.toml` (`_meta`).

## Data license & attribution

The cases in `cases/raw/` and `cases/partitioned/` derive from the
**IEEE PES Power Grid Library — Optimal Power Flow (PGLib-OPF)**,
<https://github.com/power-grid-lib/pglib-opf>. PGLib-OPF **data** is
licensed **CC-BY 4.0** (<http://creativecommons.org/licenses/by/4.0/>);
the PGLib software is MIT. The upstream `LICENSE` ships verbatim as
`cases/raw/pglib/PGLIB_LICENSE`, and each `.m` retains its per-case
copyright header. Partitioned cases are modifications (a Metis-assigned
`zone` column) of these CC-BY sources and remain CC-BY 4.0 with
attribution preserved.
