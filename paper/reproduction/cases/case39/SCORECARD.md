# case39_epri (k = 3 zones) - scorecard

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
| state / tie_line | ADMM | 10 | 742 | -0.0227 | 742 | pass |  |
| state / tie_line | ALADIN | 297 | 308 | -0.00912 | 308 | pass | glob. o |
| power / tie_line | ADMM | 4 | 1011 | -0.756 | 1011 | pass | relaxation gap |
| power / tie_line | ALADIN | 115 | 124 | -0.752 | 126 | pass | glob. o; relaxation gap |
| power / dummy_bus | ADMM | 2 | 982 | -0.744 | 982 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 27 | 2137 | -0.752 | 6692 | pass | glob. m; relaxation gap |
| hybrid / tie_line | ADMM | 74 | 7633 | +0.00153 | 7633 | pass |  |
| hybrid / tie_line | ALADIN | 112 | 136 | +2.19e-06 | 139 | pass | glob. o |
| hybrid / dummy_bus | ADMM | 288 | 24738 | -0.00144 | 24738 | pass |  |
| hybrid / dummy_bus | ALADIN | 57 | 58 | -7.16e-06 | 59 | pass | glob. o |
| hybrid / split_bus | ADMM | 105 | 17988 | +0.00198 | 17988 | pass |  |
| hybrid / split_bus | ALADIN | 17 | 62 | +5.48e-06 | 64 | pass | glob. a |

## AC

| coupling / sectioning | algorithm | Tier-1 (a) | Tier-2 (b) | gap % | iterations | verdict | notes |
|---|---|---|---|---|---|---|---|
| state / tie_line | ADMM | 39 | 721 | +0.0203 | 721 | pass |  |
| state / tie_line | ALADIN | x | x | - | - | not certified | no configuration certified in the reproduction; ADMM covers the cell |
| power / tie_line | ADMM | 36 | 55 | -0.521 | 55 | pass | relaxation gap |
| power / tie_line | ALADIN | 458 | x | - | - | not certified | relaxation gap; no configuration certified in the reproduction; ADMM covers the cell; Tier-1 value from the reproduction run, not re-run here |
| power / dummy_bus | ADMM | 28 | 35 | -0.546 | 35 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 641 | 811 | -0.551 | 813 | pass | glob. o; relaxation gap |
| hybrid / tie_line | ADMM | 117 | 749 | +0.00125 | 749 | pass |  |
| hybrid / tie_line | ALADIN | 5786 | x | - | - | not certified | no configuration certified in the reproduction; ADMM covers the cell; Tier-1 value from the reproduction run, not re-run here |
| hybrid / dummy_bus | ADMM | 307 | 837 | -0.00203 | 837 | pass |  |
| hybrid / dummy_bus | ALADIN | 1523 | 5537 | -1.01e-05 | 5598 | pass | glob. a |
| hybrid / split_bus | ADMM | 197 | 966 | +0.00135 | 966 | pass |  |
| hybrid / split_bus | ALADIN | 192 | 972 | -3.64e-05 | 1102 | pass | glob. m |

## SOCWRNL

| coupling / sectioning | algorithm | Tier-1 (a) | Tier-2 (b) | gap % | iterations | verdict | notes |
|---|---|---|---|---|---|---|---|
| state / tie_line | ADMM | 317 | 346 | -0.582 | 346 | pass |  |
| state / tie_line | ALADIN | 1211 | 1703 | +0.00105 | 1941 | pass | glob. o |
| power / tie_line | ADMM | 29 | 36 | +0.00238 | 36 | pass | relaxation gap |
| power / tie_line | ALADIN | 809 | 986 | +0.000831 | 1017 | pass | glob. m; relaxation gap |
| power / dummy_bus | ADMM | 28 | 35 | +0.00309 | 35 | pass | relaxation gap |
| power / dummy_bus | ALADIN | 78 | 84 | -0.000371 | 85 | pass | glob. m; relaxation gap |
| hybrid / tie_line | ADMM | 47 | 91 | +3.8e-05 | 91 | pass |  |
| hybrid / tie_line | ALADIN | 452 | 851 | +7.89e-06 | 855 | pass | glob. o |
| hybrid / dummy_bus | ADMM | 50 | 300 | -0.00017 | 300 | pass |  |
| hybrid / dummy_bus | ALADIN | 439 | 611 | +0.000636 | 619 | pass | glob. o |
| hybrid / split_bus | ADMM | 25 | 159 | +6.09e-05 | 159 | pass |  |
| hybrid / split_bus | ALADIN | 41 | 193 | +1.28e-05 | 256 | pass | glob. o |

## Totals

- ADMM certifies 18 of 18 cells.
- ALADIN certifies 15 of 18 cells.

A reproduction matches when, for every row, `tier1_first_touch`, `tier2_iter`, `tier2_first_touch_raw` and `verdict` equal the values in `expected_crossing.csv` (`scripts/verify_crossing.jl` does this check).
