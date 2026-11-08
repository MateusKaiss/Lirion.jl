# case24_ieee_rts (k = 3 zones) - scorecard

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
| state / tie_line | ADMM | 34 | 2317 | +6.25e-05 | 2317 | pass |  |
| state / tie_line | ALADIN | 187 | 191 | -0.0461 | 191 | pass | glob. o |
| power / tie_line | ADMM | 15 | 32 | +0.0101 | 32 | pass | relaxation gap |
| power / tie_line | ALADIN | 10 | 25 | +2.64e-06 | 30 | pass | glob. o; relaxation gap |
| power / dummy_bus | ADMM | 14 | 26 | -0.0187 | 26 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 24 | 33 | -2.95e-07 | 41 | pass | glob. a; relaxation gap |
| hybrid / tie_line | ADMM | 279 | 1416 | +0.00109 | 1416 | pass |  |
| hybrid / tie_line | ALADIN | 27 | 30 | +1.46e-05 | 31 | pass | glob. o |
| hybrid / dummy_bus | ADMM | 273 | 2011 | +0.000869 | 2011 | pass |  |
| hybrid / dummy_bus | ALADIN | 37 | 57 | +2.65e-06 | 72 | pass | glob. o |
| hybrid / split_bus | ADMM | 501 | 2830 | +0.000212 | 2830 | pass |  |
| hybrid / split_bus | ALADIN | 10 | 12 | +6.12e-06 | 13 | pass | glob. o |

## AC

| coupling / sectioning | algorithm | Tier-1 (a) | Tier-2 (b) | gap % | iterations | verdict | notes |
|---|---|---|---|---|---|---|---|
| state / tie_line | ADMM | 39 | 352 | +0.0254 | 352 | pass |  |
| state / tie_line | ALADIN | 1163 | x | - | - | not certified | no configuration certified in the reproduction; ADMM covers the cell; Tier-1 value from the reproduction run, not re-run here |
| power / tie_line | ADMM | 47 | 68 | -0.0103 | 68 | pass | relaxation gap |
| power / tie_line | ALADIN | 83 | x | - | - | not certified | relaxation gap; no configuration certified in the reproduction; ADMM covers the cell; Tier-1 value from the reproduction run, not re-run here |
| power / dummy_bus | ADMM | 47 | 64 | -0.0242 | 64 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 97 | 154 | -0.0263 | 176 | pass | glob. a; relaxation gap |
| hybrid / tie_line | ADMM | 96 | 140 | +0.00207 | 140 | pass |  |
| hybrid / tie_line | ALADIN | 175 | 415 | -0.00025 | 446 | pass | glob. o |
| hybrid / dummy_bus | ADMM | 66 | 122 | -0.000291 | 122 | pass |  |
| hybrid / dummy_bus | ALADIN | 270 | 378 | -0.0011 | 378 | pass | glob. a |
| hybrid / split_bus | ADMM | 93 | 170 | -0.00236 | 170 | pass |  |
| hybrid / split_bus | ALADIN | 18 | 72 | -0.00023 | 93 | pass | glob. a |

## SOCWRNL

| coupling / sectioning | algorithm | Tier-1 (a) | Tier-2 (b) | gap % | iterations | verdict | notes |
|---|---|---|---|---|---|---|---|
| state / tie_line | ADMM | 97 | 329 | -0.00478 | 329 | pass |  |
| state / tie_line | ALADIN | x | x | - | - | not certified | no configuration certified in the reproduction (gap gate fails); ADMM covers the cell |
| power / tie_line | ADMM | 47 | 69 | -0.00618 | 69 | pass | relaxation gap |
| power / tie_line | ALADIN | 22318 | 24218 | -0.0258 | 24218 | pass | glob. m; relaxation gap |
| power / dummy_bus | ADMM | 47 | 65 | -0.0203 | 65 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 130 | 167 | -0.0222 | 169 | pass | glob. o; relaxation gap |
| hybrid / tie_line | ADMM | 46 | 71 | +0.0029 | 71 | pass |  |
| hybrid / tie_line | ALADIN | 319 | 348 | -0.00038 | 369 | pass | glob. a |
| hybrid / dummy_bus | ADMM | 56 | 90 | +0.000107 | 90 | pass |  |
| hybrid / dummy_bus | ALADIN | 190 | 278 | -0.000193 | 280 | pass | glob. o |
| hybrid / split_bus | ADMM | 58 | 155 | -0.000286 | 155 | pass |  |
| hybrid / split_bus | ALADIN | 217 | 295 | +0.000145 | 331 | pass | glob. f |

## Totals

- ADMM certifies 18 of 18 cells.
- ALADIN certifies 15 of 18 cells.

A reproduction matches when, for every row, `tier1_first_touch`, `tier2_iter`, `tier2_first_touch_raw` and `verdict` equal the values in `expected_crossing.csv` (`scripts/verify_crossing.jl` does this check).
