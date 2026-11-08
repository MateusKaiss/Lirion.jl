# case118_ieee (k = 4 zones) - scorecard

Summary of this case's rows in `expected_crossing.csv` (the file `run.jl` regenerates and checks). One row per cell and algorithm: 3 formulations x 6 coupling/sectioning rows x {ADMM, ALADIN} = 36 rows.

## How to read a row

- **Tier-1 (a)** - first iteration at which the consensus mismatch satisfies ||r||_2 <= 0.01 (rad / p.u.). This is the mismatch-only stop used in the literature. Column `tier1_first_touch`.
- **Tier-2 (b)** - the certified iteration: primal residual under its Boyd bar and the stationarity residual under its bar (ADMM: dual residual; ALADIN: scaled step residual), with the objective within 1% of the centralized solution. Column `tier2_iter`.
- **Verdict** - `pass` = certified (Tier-1 and Tier-2 both reached, gap gate satisfied); `not certified` = what the paper prints as x. Cells that never certified in the reproduction are not solved by `run.jl` (no parameters to run); they appear as placeholders with the reproduction note in `verdict_reason`.
- **Gap** - objective relative to the centralized oracle of the same formulation, in %. For `power` coupling the interzonal voltage law is relaxed, so that gap measures relaxation slack, not solver error; those rows are gated on convergence only.
- The paper prints iteration counts >= 10000 truncated to thousands (e.g. 14763 -> 14k). This file keeps the exact integers.
- ALADIN superscript: o = original Alg. 3 (Wolfe + fallback), a = Alg. 3 golden-section, m = Armijo, f = full step. The dagger (lifted) marks the SOCWRNL ADMM cells run with the lifted-flow penalty.

## DC

| coupling / sectioning | algorithm | Tier-1 (a) | Tier-2 (b) | gap % | iterations | verdict | notes |
|---|---|---|---|---|---|---|---|
| state / tie_line | ADMM | 21 | 15434 | +0.00428 | 15434 | pass |  |
| state / tie_line | ALADIN | 113 | 188 | -0.000184 | 349 | pass | glob. a |
| power / tie_line | ADMM | 21 | 173 | +0.00107 | 173 | pass | relaxation gap |
| power / tie_line | ALADIN | 110 | 696 | -0.015 | 776 | pass | glob. m; relaxation gap |
| power / dummy_bus | ADMM | 10 | 99 | -0.0264 | 99 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 85 | 980 | -0.0169 | 5241 | pass | glob. m; relaxation gap |
| hybrid / tie_line | ADMM | 106 | 418 | +0.000919 | 418 | pass |  |
| hybrid / tie_line | ALADIN | 82 | 99 | +1.51e-06 | 103 | pass | glob. o |
| hybrid / dummy_bus | ADMM | 489 | 3582 | -5.59e-05 | 3582 | pass |  |
| hybrid / dummy_bus | ALADIN | 204 | 275 | +9.23e-07 | 276 | pass | glob. o |
| hybrid / split_bus | ADMM | 328 | 4774 | -8.86e-05 | 4774 | pass |  |
| hybrid / split_bus | ALADIN | 52 | 69 | +6.32e-06 | 75 | pass | glob. o |

## AC

| coupling / sectioning | algorithm | Tier-1 (a) | Tier-2 (b) | gap % | iterations | verdict | notes |
|---|---|---|---|---|---|---|---|
| state / tie_line | ADMM | 33 | 2211 | +0.000715 | 2211 | pass |  |
| state / tie_line | ALADIN | 732 | x | - | - | not certified | no configuration certified in the reproduction; ADMM covers the cell; Tier-1 value from the reproduction run, not re-run here |
| power / tie_line | ADMM | 38 | 69 | -0.403 | 69 | pass | relaxation gap |
| power / tie_line | ALADIN | 554 | 1283 | -0.39 | 1413 | pass | glob. o; relaxation gap |
| power / dummy_bus | ADMM | 36 | 67 | -0.406 | 67 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 195 | 327 | -0.428 | 331 | pass | glob. a; relaxation gap |
| hybrid / tie_line | ADMM | 69 | 106 | +0.000673 | 106 | pass |  |
| hybrid / tie_line | ALADIN | 511 | x | - | - | not certified | reaches Tier-1 but never certifies Tier-2 on the frozen code; Tier-1 value from the reproduction run, not re-run here |
| hybrid / dummy_bus | ADMM | 79 | 133 | -0.00243 | 133 | pass |  |
| hybrid / dummy_bus | ALADIN | 235 | 699 | -1.83e-05 | 704 | pass | glob. a |
| hybrid / split_bus | ADMM | 115 | 186 | -0.00134 | 186 | pass |  |
| hybrid / split_bus | ALADIN | 98 | 237 | -7.97e-05 | 253 | pass | glob. a |

## SOCWRNL

| coupling / sectioning | algorithm | Tier-1 (a) | Tier-2 (b) | gap % | iterations | verdict | notes |
|---|---|---|---|---|---|---|---|
| state / tie_line | ADMM | 60 | 96 | -0.0039 | 96 | pass |  |
| state / tie_line | ALADIN | 2931 | 3197 | -0.0664 | 3304 | pass | glob. m |
| power / tie_line | ADMM | 38 | 45 | -0.000812 | 45 | pass | relaxation gap |
| power / tie_line | ALADIN | 4375 | 6396 | -0.0151 | 6396 | pass | glob. o; relaxation gap |
| power / dummy_bus | ADMM | 38 | 51 | -0.0288 | 51 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 361 | 447 | -0.0212 | 447 | pass | glob. o; relaxation gap |
| hybrid / tie_line | ADMM | 47 | 53 | -0.0262 | 53 | pass |  |
| hybrid / tie_line | ALADIN | 158 | 203 | +6.45e-05 | 204 | pass | glob. o |
| hybrid / dummy_bus | ADMM | 39 | 50 | -0.00735 | 50 | pass |  |
| hybrid / dummy_bus | ALADIN | 76 | 126 | -0.00189 | 126 | pass | glob. o |
| hybrid / split_bus | ADMM | 46 | 70 | +0.000766 | 70 | pass |  |
| hybrid / split_bus | ALADIN | 98 | 228 | +4.56e-05 | 242 | pass | glob. o |

## Totals

- ADMM certifies 18 of 18 cells.
- ALADIN certifies 16 of 18 cells.

A reproduction matches when, for every row, `tier1_first_touch`, `tier2_iter`, `tier2_first_touch_raw` and `verdict` equal the values in `expected_crossing.csv` (`scripts/verify_crossing.jl` does this check).
