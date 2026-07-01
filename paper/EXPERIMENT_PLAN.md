# Experiment Plan — Lirion.jl distributed-OPF paper

Status: **design only.** No drivers, no solver runs, no tuning yet. This
document is the reviewed blueprint for `paper/scripts/experiments/`.

Ground rules for the whole stage:
- Never modify `src/`, `Project.toml`, or `Manifest.toml` for experiments.
- The partition is a single `case.m` (with a `zone` column) per
  `(case, variant, k)`; **sectioning is chosen at runtime**
  (`admm(sectioning=…)` / `aladin(sectioning=…)`), never materialized.
- Every run is recorded (including failures); nothing is silently dropped.

---

## 0. Conceptual pillars

1. **ADMM** — classical distributed consensus / decomposition.
2. **ALADIN** — distributed method better suited to **non-convex** problems.
3. **Consistency / sectioning formulations** compared in distributed OPF.
4. **PowerModels** — reference centralized power-flow formulation framework.
5. **PowerModelsADA** — reference distributed/ADA framework in Julia.

Thesis (not "we are faster than everyone"): Lirion.jl is a **verified,
traceable framework** to study *network decomposition (sectioning) ×
consistency coupling × formulation × distributed algorithm* in OPF, where:
- the **same** `bus→zone` partition drives different sectionings at runtime;
- centralized and distributed solutions are compared traceably;
- ADMM and ALADIN are evaluated in **one** experimental pipeline;
- the sectioning/consistency comparison generalizes **beyond DC** to AC and
  the SOCWR relaxation (extending Harris et al.; broader than PowerModelsADA's
  single decomposition).

---

## 1. Scientific objective

Produce a reproducible, verified results section showing how sectioning,
coupling, formulation, and algorithm interact in distributed OPF, using a
centralized baseline (itself cross-validated against PowerModels via
`LIRION_TEST_VS_PM`) as the oracle for gaps.

## 2. Hypotheses

- **H1 (convex vs non-convex):** in DC and SOCWR (convex), ADMM converges
  reliably to the centralized optimum (gap→0); in AC (non-convex), ADMM may
  stall/diverge while **ALADIN converges** (superlinear). *(pillars 1×2)*
- **H2 (shared partition):** the same `bus→zone` map under
  `tie_line`/`dummy_bus`/`split_bus` yields the same optimum where the
  sectioning is exact; ⚠️ with a **small, documented modeling gap** for
  `dummy_bus` in AC/SOCWR (line-charging redistribution — exact only in DC).
- **H3 (coupling):** coupling (`state/power/hybrid/pflow`) affects convergence
  rate/robustness, not the optimum.
- **H4 (scale):** convergence difficulty grows with `k` and the number of
  tie-lines.
- **H5 (generalization):** reproduces the DC-only sectioning comparison of
  Harris et al. and extends it to AC and SOCWR.

## 3. Tables & figures

- **T1 — Centralized baseline:** `case × variant × formulation` → status,
  objective, solve_time, solver.
- **T2 — Distributed matrix:** `(case, k, formulation, sectioning, coupling,
  algorithm)` → gap_pct, iterations, residuals, runtime, converged?.
- **T3 — ADMM vs ALADIN on AC:** where each converges (evidence for H1).
- **T4 — Sectioning comparison:** tie/dummy/split at fixed formulation/coupling
  → iterations, gap, tie-lines.
- **T5 — Scaling:** iterations/time vs `k` (and vs #buses).
- **Fig1 —** convergence curves (residual vs iteration), ADMM vs ALADIN.
- **Fig2 —** gap heatmap `sectioning × coupling` per formulation.
- **Fig3 —** scaling curve.
- **Fig4 (optional) —** DummyBus exactness: DC ~1e-12 vs small AC/SOCWR
  modeling gap (links to the existing Ybus equivalence test).

## 4. Execution order

**A** centralized baseline → **B** distributed smoke test → **C** ADMM/ALADIN
on one schema → **D** declared two-level tuning → **E** main matrix. Each stage
starts only after the previous one is validated; you review between stages.

## 5. Cases per stage

- **A (centralized):** distinct cases of the clean set — case57, case73,
  case118, case179, case300, case500, case588, case1354 (+ optional case24,
  case39). Start on `std`; expand to `api`/`sad` later. (Centralized is
  independent of `k`/sectioning.)
- **B (smoke):** `case57_ieee k2`, `case73_ieee_rts k2/k3`, `case118_ieee k2`.
- **C:** same small set, final schema, adding AC.
- **D (tuning):** two levels — see §8. D1 (base) on a representative subset;
  D2 (per-config) over the whole main matrix.
- **E (main matrix):** the **27 `clean_for_main = true`** combinations.
  The **14 `connected_no_gen`** stay **out of the headline** — optional /
  sensitivity only (connected partitions where some zone has no dispatchable
  generator and imports through tie-lines; e.g. case14, case30).

## 6. Formulations — order

**DC** and **SOCWR** first (convex, reliable, interpretable gap), then **AC**
(non-convex; where ADMM struggles and ALADIN shines). Others (QC/IVR/…) later.
⚠️ The distributed layer currently covers **DC / AC / SOCWR**.

## 7. Sectionings & couplings — order

- Sectioning: `tie_line` first (exact, simplest), then `dummy_bus`, `split_bus`.
- Coupling: `state` and `hybrid` first; then `power` and `pflow`
  (⚠️ `pflow` requires `coupling = :power` or `:hybrid`).

### Algorithm coverage
- **Both ADMM and ALADIN** are considered for the three main formulations
  **DC, AC, and SOCWR**, for **every** `(sectioning, coupling)` cell **that the
  API supports**. ALADIN is not restricted to AC — it is evaluated on the
  convex formulations too (where `rho = 0` reduces it to dual decomposition),
  so ADMM and ALADIN are compared on the same footing everywhere the
  combination is available.

### Support matrix (do not assume; discover in the smoke test)
- We do **not** assume that every
  `formulation × algorithm × sectioning × coupling` cell works. The smoke test
  (Stage B) builds an authoritative **support matrix**, saved as
  `artifacts/support_matrix.toml`, classifying each cell as one of:
  - `supported` — runs and converges as expected;
  - `not_implemented` — the API/code path does not exist
    (e.g. ⚠️ `SOCWR × dummy_bus × ALADIN`, intentionally not implemented);
  - `unstable_or_fails` — implemented but does not converge / errors under the
    declared settings (recorded with the failure reason);
  - `skipped_by_design` — deliberately excluded (e.g. `pflow` with
    `coupling = :state`, which the API forbids).
- Only `supported` cells enter the main matrix (Stage E). The others are
  reported explicitly (transparency), never silently omitted.

## 8. Tuning — two declared, reproducible levels

Calibration is **not** done only on small cases. There are two levels, both
with fixed grids and a fixed selection rule.

### Level 1 — Global/base tuning (defaults)
- **Purpose:** good **initial defaults** per `(formulation, algorithm,
  coupling)` — a starting point, not a final result.
- **How:** a **coarse** declared grid
  (`configs/tuning_grid_admm_base.toml`, `…_aladin_base.toml`) run on a
  representative subset (e.g. case57, case118, case300 at clean `k`).
- **Output:** `artifacts/tuning/defaults.toml` per `(formulation, algorithm,
  coupling)` — initial ρ/μ/growth/adaptive. Used to (a) center the Level-2
  search and (b) act as fallback when Level 2 finds nothing better.

### Level 2 — Per-configuration tuning (final paper results)
- **Scope:** for **each** applicable cell
  `case × variant × k × formulation × algorithm × sectioning × coupling`,
  run the **same** fixed grid (`configs/tuning_grid_admm.toml` /
  `…_aladin.toml`) — a **focused ladder** around the Level-1 default (not a
  blind cartesian product; see cost control).
- **Same search logic everywhere** — small, medium, and large cases alike, no
  manual exceptions.
- **Save every attempt** (including failures) to `results/tuning/` (local).
- **Automatic winner** by the declared rule below; feeds the main matrix.

### Selection rule (lexicographic, declared before running)
Rank attempts of a cell by:
1. **converged** — `primal_l2 < eps_abs` **and** `dual_l2 < eps_abs` within
   `max_iter`;
2. **gap** — `gap_pct < gap_target` (vs the same-formulation centralized);
3. **smallest residual** — combined primal+dual norm at the final point;
4. **fewest iterations**, then **lowest runtime** as the final tiebreak.

> Iterations are prioritized over runtime because iteration count is a cleaner,
> machine-independent measure of algorithm behavior (runtime varies with
> machine, solver, precompile, environment). Runtime is only the last tiebreak.

If **no** attempt converges, the "least bad" is selected by the same order,
flagged `convergence_flag = false`, and kept visible in the matrix.

### Per-case overrides
- Allowed only via `configs/tuning_overrides.toml`, and **each override
  requires** a `reason` + `author`/`date`. No documented reason → no override.
- An override appears in the result row as `tuning_level = "override"` with the
  `override_reason` (auditable).

### Cost control (declared, not cherry-picking)
The full cartesian product explodes, so to stay **declared + tractable**:
- Level 1 narrows the grid → Level 2 uses a **focused ladder** (e.g. 5–7 ρ
  values around the default × {adaptive on/off}), declared in the config.
- **Declared, optional early-stop:** the search may stop at the first parameter
  meeting (1)+(2), **but logs everything it tried** and the stop flag;
  `EXHAUSTIVE = true` forces the whole grid. The **same** policy applies to
  every cell.
- Any pruning/limit is logged (never silently truncated).

## 9. Gap vs centralized

`gap_pct = 100 · (obj_dist − obj_cen) / |obj_cen|`, always comparing the **same
formulation** (DC dist vs DC cen, AC vs AC, SOCWR vs SOCWR).
- **Exact** sectionings (`tie_line` always; `split_bus`; `dummy_bus` in DC): at
  convergence the gap is purely the convergence gap → should tend to 0.
- ⚠️ `dummy_bus` in **AC/SOCWR**: there is an **inherent modeling gap**
  (line-charging). To separate modeling from convergence, also use as oracle
  the **centralized optimum of the dummy-bus-augmented network** (not only the
  original). Record **both** numbers and label the gap's nature.
- The centralized reference itself is trustworthy because it is cross-validated
  against PowerModels (`LIRION_TEST_VS_PM`).

## 10. Convergence-failure handling

- `max_iter`, `eps_abs`, `eps_rel` declared a priori (in `configs/`).
- Every run produces a row with `convergence_flag` + `failure_reason`;
  **never dropped**.
- Distinguish: (a) **algorithm non-convergence** (hit `max_iter` without eps),
  (b) **subproblem solver failure** (infeasible/numerical),
  (c) **divergence**.
- `time_limit` per run; a partial result (best iterate) is recorded with flags.

## 11. Anti-cherry-picking policy

- Grids and criteria fixed **before** running (`configs/*`).
- The **same** search logic for every applicable cell (§8).
- **All** attempts saved (full tuning history), including failures.
- `parameter_set` logged per run; fixed seeds.
- Overrides only with a documented reason.
- The **entire** matrix is reported, failures included; same `eps/max_iter`
  across compared methods.

## 12. Result schema (single row, shared by all stages)

```
case, variant, k, formulation, algorithm, sectioning, coupling,
status, objective, centralized_objective, gap_pct,
iterations, primal_residual, dual_residual, runtime_s, solver,
parameter_set, convergence_flag, failure_reason,
tuning_level (base|per_config|override), grid_id, attempt_id,
is_selected (bool), n_attempts_in_cell, early_stopped (bool),
override_reason, case_path, timestamp
```

- Centralized rows set `algorithm = centralized` and leave
  sectioning/coupling/iterations/residuals empty.
- ADMM and ALADIN emit the **same** schema → directly comparable tables.

## 13. Directory layout

### `paper/scripts/experiments/`
```
experiments/
  lib/
    schema.jl     # result-row struct + CSV/TOML writers (one schema)
    io.jl         # load partitioned <case>/<variant>/k<k>/case.m + meta.toml;
                  # resolve solver per formulation
    runners.jl    # solve() wrappers: centralized / admm / aladin + gap calc
    metrics.jl    # residuals, gap, convergence
    tuning.jl     # shared search + selection engine (used by both levels)
  run_centralized.jl          # Stage A
  run_distributed_smoke.jl    # Stage B
  run_admm_grid.jl            # Stages C/E (ADMM) — consumes Level-2 winners
  run_aladin_grid.jl          # Stages C/E (ALADIN) — consumes Level-2 winners
  run_tuning.jl               # Stage D: --level base | --level per_config
```
Solver per formulation: DC → OSQP/Gurobi, AC → Ipopt, SOCWR → Clarabel/Gurobi.
`run_tuning.jl` has two modes but **one** search/selection engine (`lib/tuning.jl`),
so the logic is identical for base and per-config levels.

### `paper/configs/` (to add / adjust)
```
tuning_grid_admm_base.toml     tuning_grid_aladin_base.toml   # Level 1 (coarse)
tuning_grid_admm.toml          tuning_grid_aladin.toml        # Level 2 (focused)
tuning_overrides.toml          # overrides, reason required
selection_criteria.toml        # rewrite [protocol]: tuning_set = Level-1 subset;
                               # NO more "transfer params to large cases"
```

### `paper/artifacts/` — versioned (curated, reproducible)
```
artifacts/
  grids/               # copies of the declared grids used for a run
  support_matrix.toml  # smoke-test cell classification (supported/not_implemented/…)
  tuning/
    defaults.toml   # Level-1 defaults per (formulation, algorithm, coupling)
    winners.toml    # Level-2 winning config per cell + convergence/failure flags
  tables/           # final aggregated T1–T5
  figures/          # final Fig1–Fig4
```

### `paper/results/` — local scratch, git-ignored
```
results/
  centralized/   # baseline runs (raw)
  distributed/   # smoke + grid runs (raw, unified schema)
  tuning/        # FULL attempt history (large) — local only
  logs/          # per-run logs and per-iteration residual traces
```

### Commit vs ignore
- **Git receives** curated, reproducible artifacts under `paper/artifacts/`:
  declared grids, Level-1 defaults, Level-2 summary, winning config per cell,
  convergence/failure flags, final aggregated tables, final figures.
- **`paper/results/` stays local/ignored** (scratch: logs, large dumps, full raw
  tuning history). The root `.gitignore` already ignores `**/results/`; keeping
  curated output in `paper/artifacts/` avoids touching that ignore rule.

---

## Open items to confirm before writing drivers

1. **Dummy-bus AC/SOCWR oracle:** record the gap against **both** the original
   centralized and the dummy-bus-augmented centralized (recommended).
2. **Support matrix:** the smoke test (Stage B) produces
   `artifacts/support_matrix.toml`, classifying every
   `(formulation, algorithm, sectioning, coupling)` cell as `supported`,
   `not_implemented`, `unstable_or_fails`, or `skipped_by_design`. ADMM **and**
   ALADIN are covered for DC/AC/SOCWR wherever the API supports the cell. Only
   `supported` cells feed Stage E; the rest are reported, not omitted.
3. **`paper/artifacts/`** is the versioned home for curated outputs; drivers
   will be added only after this plan is approved.
