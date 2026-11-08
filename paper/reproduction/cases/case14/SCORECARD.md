# case14_ieee (k = 2 zones) - scorecard

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
| state / tie_line | ADMM | 17 | 674 | -0.137 | 674 | pass |  |
| state / tie_line | ALADIN | 3 | 3 | -2.46e-05 | 3 | pass | glob. o |
| power / tie_line | ADMM | 12 | 33 | +0.00914 | 33 | pass | relaxation gap |
| power / tie_line | ALADIN | 3 | 3 | +2.02e-10 | 3 | pass | glob. m; relaxation gap |
| power / dummy_bus | ADMM | 10 | 28 | +0.00413 | 28 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 3 | 3 | -4.52e-07 | 3 | pass | glob. o; relaxation gap |
| hybrid / tie_line | ADMM | 22 | 83 | +0.00121 | 83 | pass |  |
| hybrid / tie_line | ALADIN | 3 | 3 | +2.19e-08 | 3 | pass | glob. m |
| hybrid / dummy_bus | ADMM | 26 | 187 | +0.000563 | 187 | pass |  |
| hybrid / dummy_bus | ALADIN | 3 | 3 | +3.17e-08 | 3 | pass | glob. m |
| hybrid / split_bus | ADMM | 18 | 220 | +0.00158 | 220 | pass |  |
| hybrid / split_bus | ALADIN | 3 | 3 | -5.87e-07 | 3 | pass | glob. m |

## AC

| coupling / sectioning | algorithm | Tier-1 (a) | Tier-2 (b) | gap % | iterations | verdict | notes |
|---|---|---|---|---|---|---|---|
| state / tie_line | ADMM | 17 | 3590 | -5.56e-05 | 3590 | pass |  |
| state / tie_line | ALADIN | 265 | 1246 | -0.000124 | 1921 | pass | glob. a |
| power / tie_line | ADMM | 12 | 578 | -0.00721 | 578 | pass | relaxation gap |
| power / tie_line | ALADIN | 24 | 105 | -0.00707 | 105 | pass | glob. a; relaxation gap |
| power / dummy_bus | ADMM | 10 | 4406 | -0.0435 | 4406 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 607 | 1100 | -0.0614 | 1100 | pass | glob. o; relaxation gap |
| hybrid / tie_line | ADMM | 25 | 115 | +0.000888 | 115 | pass |  |
| hybrid / tie_line | ALADIN | 88 | 215 | -9.13e-05 | 215 | pass | glob. a |
| hybrid / dummy_bus | ADMM | 29 | 222 | +7.35e-05 | 222 | pass |  |
| hybrid / dummy_bus | ALADIN | 194 | 610 | -0.000244 | 610 | pass | glob. a |
| hybrid / split_bus | ADMM | 34 | 249 | -0.000831 | 249 | pass |  |
| hybrid / split_bus | ALADIN | 8 | 20 | -9.14e-05 | 20 | pass | glob. a |

## SOCWRNL

| coupling / sectioning | algorithm | Tier-1 (a) | Tier-2 (b) | gap % | iterations | verdict | notes |
|---|---|---|---|---|---|---|---|
| state / tie_line | ADMM | 93 | 175 | +0.0406 | 175 | pass |  |
| state / tie_line | ALADIN | 467 | 625 | -0.145 | 656 | pass | glob. m |
| power / tie_line | ADMM | 32 | 69 | -0.0104 | 69 | pass | relaxation gap |
| power / tie_line | ALADIN | 1557 | 2003 | -0.0258 | 2010 | pass | glob. o; relaxation gap |
| power / dummy_bus | ADMM | 31 | 54 | -0.0397 | 54 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 152 | 206 | -0.0544 | 208 | pass | glob. o; relaxation gap |
| hybrid / tie_line | ADMM | 29 | 36 | -0.0534 | 36 | pass |  |
| hybrid / tie_line | ALADIN | 67 | 96 | -6.09e-05 | 96 | pass | glob. o |
| hybrid / dummy_bus | ADMM | 37 | 77 | +0.00429 | 77 | pass |  |
| hybrid / dummy_bus | ALADIN | 66 | 132 | +0.000232 | 132 | pass | glob. o |
| hybrid / split_bus | ADMM | 16 | 69 | +0.00216 | 69 | pass |  |
| hybrid / split_bus | ALADIN | 64 | 221 | +2.56e-05 | 235 | pass | glob. a |

## Totals

- ADMM certifies 18 of 18 cells.
- ALADIN certifies 18 of 18 cells.

A reproduction matches when, for every row, `tier1_first_touch`, `tier2_iter`, `tier2_first_touch_raw` and `verdict` equal the values in `expected_crossing.csv` (`scripts/verify_crossing.jl` does this check).
