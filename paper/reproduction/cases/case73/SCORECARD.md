# case73_ieee_rts (k = 3 zones) - scorecard

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
| state / tie_line | ADMM | 9 | 14763 | +7.6e-05 | 14763 | pass |  |
| state / tie_line | ALADIN | 18 | 27 | +3.25e-06 | 31 | pass | glob. o |
| power / tie_line | ADMM | 4 | 320 | +1.42e-05 | 320 | pass | relaxation gap |
| power / tie_line | ALADIN | 12 | 16 | +3.24e-07 | 18 | pass | glob. o; relaxation gap |
| power / dummy_bus | ADMM | 4 | 320 | +1.42e-05 | 320 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 18 | 18 | -8.57e-07 | 19 | pass | glob. o; relaxation gap |
| hybrid / tie_line | ADMM | 75 | 701 | -9.01e-05 | 701 | pass |  |
| hybrid / tie_line | ALADIN | 16 | 16 | +3.96e-05 | 16 | pass | glob. o |
| hybrid / dummy_bus | ADMM | 314 | 1733 | -4.67e-05 | 1733 | pass |  |
| hybrid / dummy_bus | ALADIN | 26 | 28 | +3.5e-05 | 28 | pass | glob. o |
| hybrid / split_bus | ADMM | 15 | 1258 | +5.29e-05 | 1258 | pass |  |
| hybrid / split_bus | ALADIN | 9 | 10 | +5.1e-05 | 10 | pass | glob. a |

## AC

| coupling / sectioning | algorithm | Tier-1 (a) | Tier-2 (b) | gap % | iterations | verdict | notes |
|---|---|---|---|---|---|---|---|
| state / tie_line | ADMM | 24 | 81 | -0.046 | 81 | pass |  |
| state / tie_line | ALADIN | 132 | x | - | - | not certified | no configuration certified in the reproduction; ADMM covers the cell; Tier-1 value from the reproduction run, not re-run here |
| power / tie_line | ADMM | 2 | 41 | -0.00299 | 41 | pass | relaxation gap |
| power / tie_line | ALADIN | 105 | 3419 | -0.00226 | 3430 | pass | glob. o; relaxation gap |
| power / dummy_bus | ADMM | 2 | 48 | -0.00282 | 48 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 91 | 161 | -0.00207 | 163 | pass | glob. o; relaxation gap |
| hybrid / tie_line | ADMM | 56 | 105 | -0.000442 | 105 | pass |  |
| hybrid / tie_line | ALADIN | 222 | 620 | -0.000233 | 686 | pass | glob. a |
| hybrid / dummy_bus | ADMM | 67 | 119 | -0.000403 | 119 | pass |  |
| hybrid / dummy_bus | ALADIN | 109 | 255 | -0.00167 | 255 | pass | glob. o |
| hybrid / split_bus | ADMM | 70 | 112 | -0.000198 | 112 | pass |  |
| hybrid / split_bus | ALADIN | 35 | 93 | -0.000231 | 132 | pass | glob. a |

## SOCWRNL

| coupling / sectioning | algorithm | Tier-1 (a) | Tier-2 (b) | gap % | iterations | verdict | notes |
|---|---|---|---|---|---|---|---|
| state / tie_line | ADMM | 36 | 60 | +0.000504 | 60 | pass |  |
| state / tie_line | ALADIN | 5136 | 5963 | -0.0517 | 5963 | pass | glob. o |
| power / tie_line | ADMM | 45 | 59 | -0.000195 | 59 | pass | relaxation gap |
| power / tie_line | ALADIN | 6103 | 7578 | -0.000226 | 7769 | pass | glob. o; relaxation gap |
| power / dummy_bus | ADMM | 36 | 41 | -0.000614 | 41 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 68 | 328 | -0.000337 | 332 | pass | glob. o; relaxation gap |
| hybrid / tie_line | ADMM | 44 | 53 | -0.000298 | 53 | pass |  |
| hybrid / tie_line | ALADIN | 113 | 235 | +8.3e-05 | 236 | pass | glob. a |
| hybrid / dummy_bus | ADMM | 38 | 61 | +2.51e-05 | 61 | pass |  |
| hybrid / dummy_bus | ALADIN | 97 | 125 | +0.000324 | 129 | pass | glob. o |
| hybrid / split_bus | ADMM | 36 | 58 | +1.45e-05 | 58 | pass |  |
| hybrid / split_bus | ALADIN | 67 | 142 | +1.32e-05 | 202 | pass | glob. o |

## Totals

- ADMM certifies 18 of 18 cells.
- ALADIN certifies 17 of 18 cells.

A reproduction matches when, for every row, `tier1_first_touch`, `tier2_iter`, `tier2_first_touch_raw` and `verdict` equal the values in `expected_crossing.csv` (`scripts/verify_crossing.jl` does this check).
