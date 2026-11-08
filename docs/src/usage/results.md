# Results

Every call to `solve` hands back a single `LirionResult`, no matter which model or algorithm you ran. It holds two things: a short summary of the solution (cost, status, timing) and the full per-element solution as `DataFrame`s whose columns follow the MATPOWER schema of the case you loaded.

The same type comes back from a centralized run and a distributed one. A centralized run leaves the distributed fields empty, and a distributed run fills them in. The two sections below split the fields along that line.

The centralized example uses `pglib_opf_case24_ieee_rts.m` from the [PGLib-OPF library](https://github.com/power-grid-lib/pglib-opf); the distributed one uses a three-zone partition of the same case, shipped with Lirion at `paper/cases/partitioned/case24_ieee_rts/std/k3/case.m`. If you have not set up the solvers yet, see the [Installation Guide](@ref).

## Centralized

A centralized solve populates the summary fields and the three solution tables.

```julia
using Lirion
result = solve("test/matpower_cases/pglib_opf_case24_ieee_rts.m"; model = ac())

result.objective    # total generation cost (Float64)
result.status       # :optimal, :local_optimal, :infeasible, ...
result.solve_time   # wall-clock seconds
```

`status` is read back from the solver. The interior-point models (AC, ACR, IVR, LPAC, QC, SOCWRNL on Ipopt) report `:local_optimal`, and the convex ones (DC, SOCWR, and the other relaxations) report `:optimal`. You may also see `:almost_optimal`, `:iteration_limit`, `:time_limit`, or `:infeasible`.

Four more summary fields record what was solved and how: `result.algorithm` (`:centralized`, or `:admm`/`:aladin` for a distributed run), `result.base_model` (the model, such as `:AC`, `:DC`, or `:SOCWR`), `result.coupling`, and `result.sectioning`. A centralized run leaves the last two at `:none`; a distributed run sets `coupling` to `:state`, `:power`, or `:hybrid` and `sectioning` to `:tieline`, `:dummybus`, or `:splitbus`. All four also go into `summary.csv`.

### Solution tables

`result.bus`, `result.gen`, and `result.branch` are `DataFrame`s. Which columns show up depends on what the formulation models, so a DC result has no voltage magnitude and no reactive power while an AC result has both.

`result.bus` has one row per bus.

| column | present | meaning |
| --- | --- | --- |
| `bus_i` | always | bus id |
| `zone` | always | zone label from the case |
| `pd` | always | active power demand |
| `qd` | always | reactive power demand |
| `va` | models with an angle | voltage angle |
| `vm` | models that track magnitude | voltage magnitude |
| `vm_sq` | SOCWR-style relaxations | squared voltage magnitude |

`result.gen` has one row per in-service generator. Only generators that are both flagged in service and present in the solved model appear.

| column | present | meaning |
| --- | --- | --- |
| `gen_id` | always | generator index |
| `gen_bus` | always | bus it sits on |
| `status` | always | in-service flag |
| `pg` | always | active dispatch |
| `qg` | models with reactive dispatch | reactive dispatch |
| `curtailment` | models with curtailment | curtailed power |

`result.branch` has one row per in-service branch, with flows reported at both ends.

| column | present | meaning |
| --- | --- | --- |
| `branch_id` | always | branch index |
| `from` | always | from-bus |
| `to` | always | to-bus |
| `status` | always | in-service flag |
| `pf` | always | active flow, from end |
| `pt` | always | active flow, to end |
| `qf` | models with reactive flow | reactive flow, from end |
| `qt` | models with reactive flow | reactive flow, to end |

Because these are plain `DataFrame`s, you work with them the usual way:

```julia
result.objective                           # the number you optimised
sum(result.gen.pg)                         # total active generation
result.gen[result.gen.pg .> 0.5, :]        # generators above 0.5 p.u.
result.bus[result.bus.vm .< 0.95, :]       # buses under 0.95 p.u.
```

To write everything to disk, `Lirion.to_csvs(result, "out/")` dumps a `summary.csv` plus `bus.csv`, `gen.csv`, and `branch.csv` (it needs `using CSV`).

## Decentralized

A distributed solve returns the same three top-level tables, reconciled into one network-wide view that lines up with a centralized result. On top of that it fills the fields that only make sense when the grid was split across zones. You can get a distributed result from any `admm(...)` or `aladin(...)` run.

```julia
result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = ac(),
    algorithm = admm(coupling = :hybrid, sectioning = :tie_line, max_iter = 350),
)

result.iter_done        # iterations run (Int; nothing for centralized)
length(result.zones)    # number of zones
result.consensus        # boundary-variable agreement
result.diagnostics      # convergence history
```

A distributed run's `status` is `:converged`, `:iteration_limit` (it reached `max_iter`), or, for ALADIN, `:stalled` (it stopped making progress short of convergence), taken from the algorithm's own termination rather than a solver.

The top-level `bus`, `gen`, and `branch` are the reconciled network. Shared buses have their voltage columns averaged across the zones that see them, and a tie-line split across zones is fused back into a single row. The merged branch table keys on `branch_id_global`, the branch's index in the original pre-partition case.

### Per-zone views

`result.zones` is a `Vector` of one `NamedTuple` per zone. Each entry carries the zone's own local solution:

| field | meaning |
| --- | --- |
| `id` | zone id |
| `bus` | that zone's bus table |
| `gen` | that zone's generator table |
| `branch` | that zone's branch table |
| `cost` | the zone's local generation cost |

The `bus`, `gen`, and `branch` tables inside a zone use the same columns as the centralized tables above (the zone branch table also carries `branch_id_global`). The total objective is the sum of the per-zone costs.

```julia
result.zones[1].cost     # zone 1's local generation cost
result.zones[1].bus      # zone 1's buses, before reconciliation
result.zones[1].gen
```

### Consensus

`result.consensus` is a `DataFrame` summarising how well neighbouring zones agree on their shared boundary variables. It is in long format with one row per boundary variable, element, and zone.

| column | meaning |
| --- | --- |
| `variable` | the boundary quantity, a `Symbol` such as `:theta`, `:v_mag`, `:v2`, `:p_flow`, `:q_flow`, `:p_inj`, or `:q_inj` (and `:p_flow_dummy` / `:q_flow_dummy` under `:dummy_bus`) |
| `bus_id` | bus the reading is on (`missing` for a branch quantity) |
| `branch_id` | branch the reading is on (`missing` for a bus quantity) |
| `zone` | the zone that reported this value |
| `value` | what this zone saw for the variable |
| `consensus` | the mean across all zones sharing the element |
| `deviation` | `value - consensus`, how far this zone is off the agreed value |

At convergence the `deviation` column sits near zero across the board. It is the quickest way to see where the zones have not quite agreed.

### Internals

A few more fields hold the objects behind a result, rarely needed day to day.

- `result.model` is the underlying `JuMP` model for a centralized run (`nothing` for a distributed one), and `result.system` is the parsed system: a `SystemData` for a centralized run, the partitioned system for a distributed one.
- `result.partitions` is the `Partitions` object with the per-zone models and the index maps between zone-local and global numbering (`nothing` for a centralized run).
- `result.diagnostics` is the iteration-by-iteration history. ADMM returns an `ADMMHistory` with `r_primal`, `r_dual`, `eps_primal`, `eps_dual`, the penalty track `rhos`, and a `meta` dictionary. ALADIN returns an `ALADINHistory`, which adds the step residuals `r_step_primal` and `eps_step_primal`, the dual update `lambda_step`, and the coordinator penalty `mus`. A centralized run carries an empty `NoDiagnostics`.

Each residual field is a vector with one entry per iteration, so you can plot the convergence directly:

```julia
result.diagnostics.r_primal     # primal residual per iteration
result.diagnostics.r_dual       # dual residual per iteration
```
