# case500_goc (k = 20 zones) - scorecard

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
| state / tie_line | ADMM | 204 | 5006 | -0.48 | 5006 | pass |  |
| state / tie_line | ALADIN | 1336 | 39878 | -2.26e-05 | 56930 | pass | glob. m |
| power / tie_line | ADMM | 48 | 134 | -0.152 | 134 | pass | relaxation gap |
| power / tie_line | ALADIN | 18258 | 22612 | -0.188 | 22612 | pass | glob. o; relaxation gap |
| power / dummy_bus | ADMM | 86 | 104 | -0.143 | 104 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 16164 | 19628 | -0.177 | 19628 | pass | glob. o; relaxation gap |
| hybrid / tie_line | ADMM | 8896 | 10458 | -0.0139 | 10458 | pass |  |
| hybrid / tie_line | ALADIN | 3545 | 3672 | +2.29e-07 | 3714 | pass | glob. o |
| hybrid / dummy_bus | ADMM | 11617 | 35179 | -0.00127 | 35179 | pass |  |
| hybrid / dummy_bus | ALADIN | x | x | - | - | not certified | stalls before the certificate (best attempt noted in expected_crossing.csv) |
| hybrid / split_bus | ADMM | 15001 | 40360 | +0.00141 | 40360 | pass |  |
| hybrid / split_bus | ALADIN | 135 | 240 | +2.36e-07 | 299 | pass | glob. o |

## AC

| coupling / sectioning | algorithm | Tier-1 (a) | Tier-2 (b) | gap % | iterations | verdict | notes |
|---|---|---|---|---|---|---|---|
| state / tie_line | ADMM | x | x | - | - | not certified | stalls before the certificate (best attempt noted in expected_crossing.csv) |
| state / tie_line | ALADIN | x | x | - | - | not certified | stalls before the certificate (best attempt noted in expected_crossing.csv) |
| power / tie_line | ADMM | 132 | 299 | -0.243 | 299 | pass | relaxation gap |
| power / tie_line | ALADIN | x | x | - | - | not certified | relaxation gap; stalls before the certificate (best attempt noted in expected_crossing.csv) |
| power / dummy_bus | ADMM | 4909 | 5649 | -0.348 | 5649 | pass | relaxation gap |
| power / dummy_bus | ALADIN | x | x | - | - | not certified | relaxation gap; stalls before the certificate (best attempt noted in expected_crossing.csv) |
| hybrid / tie_line | ADMM | 1539 | 1640 | -0.048 | 1640 | pass |  |
| hybrid / tie_line | ALADIN | 6273 | x | - | - | not certified | reaches Tier-1 but the stationarity bar is never met; Tier-1 value from the reproduction run, not re-run here |
| hybrid / dummy_bus | ADMM | 510 | 622 | -0.0942 | 622 | pass |  |
| hybrid / dummy_bus | ALADIN | x | x | - | - | not certified | stalls before the certificate (best attempt noted in expected_crossing.csv) |
| hybrid / split_bus | ADMM | 546 | 715 | -0.0928 | 715 | pass |  |
| hybrid / split_bus | ALADIN | 116 | x | - | - | not certified | the only certifying run used an off-grid variant; kept uncertified; Tier-1 value from the reproduction run, not re-run here |

## SOCWRNL

| coupling / sectioning | algorithm | Tier-1 (a) | Tier-2 (b) | gap % | iterations | verdict | notes |
|---|---|---|---|---|---|---|---|
| state / tie_line | ADMM | x | x | - | - | not certified | stalls before the certificate (best attempt noted in expected_crossing.csv) |
| state / tie_line | ALADIN | x | x | - | - | not certified | stalls before the certificate (best attempt noted in expected_crossing.csv) |
| power / tie_line | ADMM | 14717 | 16348 | -0.144 | 16348 | pass | relaxation gap |
| power / tie_line | ALADIN | x | x | - | - | not certified | relaxation gap; stalls before the certificate (best attempt noted in expected_crossing.csv) |
| power / dummy_bus | ADMM | 136 | 317 | -0.0789 | 317 | pass | lifted; relaxation gap |
| power / dummy_bus | ALADIN | 11k | x | - | - | not certified | relaxation gap; a late run certified but its extended confirmation was never completed; kept uncertified; Tier-1 value from the reproduction run, not re-run here |
| hybrid / tie_line | ADMM | 209 | 329 | +0.00272 | 329 | pass | lifted |
| hybrid / tie_line | ALADIN | 198 | x | - | - | not certified | reaches Tier-1 but the stationarity bar is never met; Tier-1 value from the reproduction run, not re-run here |
| hybrid / dummy_bus | ADMM | 260 | 386 | +0.000619 | 386 | pass | lifted |
| hybrid / dummy_bus | ALADIN | 11707 | 13046 | -0.0331 | 13046 | pass | glob. o |
| hybrid / split_bus | ADMM | 358 | 390 | -0.0679 | 390 | pass |  |
| hybrid / split_bus | ALADIN | 988 | 1308 | -0.00703 | 1308 | pass | glob. o |

## Totals

- ADMM certifies 16 of 18 cells.
- ALADIN certifies 7 of 18 cells.

A reproduction matches when, for every row, `tier1_first_touch`, `tier2_iter`, `tier2_first_touch_raw` and `verdict` equal the values in `expected_crossing.csv` (`scripts/verify_crossing.jl` does this check).
