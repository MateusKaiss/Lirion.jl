# ADMM settings

`admm(...)` builds the ADMM algorithm you hand to `solve` through the `algorithm` keyword. It runs on the `dc()`, `ac()`, `socwr()`, and `socwrnl()` models. The call takes a handful of knobs and leaves the rest at their defaults:

```julia
admm(; coupling = :hybrid, sectioning = :tie_line, max_iter = 350, kwargs...)
```

A first run can keep every default and only set the penalty:

```julia
using Lirion

result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = ac(),
    algorithm = admm(coupling = :hybrid, sectioning = :tie_line, rho = 1e3),
)
```

Each zone solves its own subproblem with a default solver chosen for the model. To point every zone at a different solver, pass `optimizer` to `solve` — it is a `solve` keyword, not an `admm` one — exactly as you would for a centralized run.

## Coupling

`coupling` picks what neighbouring zones agree on at their shared boundary. It is one of three symbols, and the default is `:hybrid`.

- `:state` couples the boundary bus state each model carries: the angle on DC, the angle and magnitude on AC, and, on SOCWR and SOCWRNL, the W-space terms — the squared voltage magnitude plus the branch coupling terms `R` and `T`.
- `:power` couples the power crossing the boundary: active only on DC, active and reactive on AC, SOCWR, and SOCWRNL.
- `:hybrid` couples both the state and the power.

Any other value throws an `ArgumentError` that lists the three allowed symbols.

## Sectioning

`sectioning` picks how the network is cut into zones. The default is `:tie_line`.

- `:tie_line` cuts the branches that run between zones and shares the buses at their ends.
- `:dummy_bus` splits each inter-zone branch at an added dummy bus.
- `:split_bus` splits a shared bus and makes the zones that meet there agree on its injection.

Not every coupling works with every cut. The accepted combinations are:

| `coupling` | `:tie_line` | `:dummy_bus` | `:split_bus` |
| --- | --- | --- | --- |
| `:state` | yes | no | no |
| `:power` | yes | yes | no |
| `:hybrid` | yes | yes | yes |

`:split_bus` requires `:hybrid`, and `:dummy_bus` requires `:power` or `:hybrid`. Anything outside the table throws an `ArgumentError` that gives the reason and the allowed options, so a split can never run silently under-coupled.

## Penalty (`rho`)

ADMM pushes the zones into agreement with a quadratic penalty on their boundary mismatch, and each coupled quantity carries its own penalty weight. The simplest way to set them is the generic `rho`, which fills every weight the `(model, coupling)` pair needs with the same value:

```julia
admm(coupling = :hybrid, rho = 1e3)
```

For finer control, set the per-block weights by name. Which names apply depends on the model and the coupling:

| Model | `:state` | `:power` | `:hybrid` |
| --- | --- | --- | --- |
| `dc()` | `rho_theta` | `rho_p` | `rho_theta`, `rho_p` |
| `ac()` | `rho_theta`, `rho_v` | `rho_p`, `rho_q` | `rho_theta`, `rho_v`, `rho_p`, `rho_q` |
| `socwr()`, `socwrnl()` | `rho_v2`, `rho_r`, `rho_t` | `rho_p`, `rho_q` | `rho_v2`, `rho_p`, `rho_q`, `rho_r`, `rho_t` |

Pass either the generic `rho`, or every per-block weight for the combination. All weights must be positive.

```julia
result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = ac(),
    algorithm = admm(
        coupling = :hybrid,
        sectioning = :tie_line,
        rho_theta = 1e3,
        rho_v = 1e3,
        rho_p = 1e3,
        rho_q = 1e3,
        adaptive_rho = true,
    ),
)
```

Passing a weight that does not belong to the pair throws an error that prints the allowed set. The check runs when you solve, not when you build the algorithm:

```julia
# AC + :power accepts only rho_p and rho_q, so this errors on solve
# and the message lists the allowed names.
solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = ac(),
    algorithm = admm(coupling = :power, rho_theta = 1e3),
)
```

SOCWR (and SOCWRNL) work in W-space, so they reject `rho_theta` and `rho_v` outright and direct you to the W-space weights (`rho_v2`, `rho_r`, `rho_t`) instead.

## Adaptive penalty

A fixed `rho` is often enough, but the best value is problem-dependent and a poor guess shows up as slow convergence: set it too high and the zones agree almost at once while the shared value barely moves, so the cost crawls down; set it too low and each zone optimizes freely but the zones take forever to agree. Passing `adaptive_rho = true` (the default is `false`) lets Lirion retune the penalty as the run proceeds, which makes the outcome far less sensitive to the `rho` you start from. This is why the knob exists: it turns the one number you most have to guess into something the algorithm manages for you.

The rule is the standard ADMM residual-balancing scheme, applied to each penalty block on its own. Every iteration past a warm-up, Lirion compares that block's **primal** residual `r` (how much the zones still disagree) with its **dual** residual `s` (how much the consensus value moved between iterations):

- if the primal residual dominates (`r > mu_incr * s`), consensus is lagging, so it multiplies that block's `rho` by `tau_incr` to pull the copies together harder;
- if the dual residual dominates (`s > mu_decr * r`), the penalty is dragging the consensus around, so it multiplies `rho` by `tau_decr` to loosen it;
- otherwise it leaves `rho` where it is.

Two guards keep this from thrashing. It does nothing for the first `adapt_warmup` iterations, where the residuals are still transient and would only mislead it, and it freezes once both residuals fall within `adapt_freeze ×` their stopping tolerances, so the final stretch converges on a steady penalty instead of chasing a moving one. Because each coupled quantity is balanced from its own residuals, a hybrid run can end up with a different `rho` on the state blocks than on the power blocks, which is exactly what you want when they converge at different rates.

Every knob has a default you rarely need to touch:

| knob | default | what it does |
| --- | --- | --- |
| `adaptive_rho` | `false` | turn the retuning on |
| `adapt_warmup` | `5` | iterations to wait before the first adjustment |
| `mu_incr` | `30.0` | raise `rho` when the primal residual exceeds `mu_incr ×` the dual |
| `mu_decr` | `10.0` | lower `rho` when the dual residual exceeds `mu_decr ×` the primal |
| `tau_incr` | `1.2` | factor `rho` is multiplied by when raised |
| `tau_decr` | `0.5` | factor `rho` is multiplied by when lowered |
| `rho_floor` | `0.0` | lower bound `rho` is never dropped below (`0` disables the floor) |
| `adapt_freeze` | `100.0` | freeze the retuning once both residuals are within this multiple of their tolerances |

The common case is to switch it on and leave the rest alone; the starting `rho` is still the point the retuning departs from, so it is worth setting:

```julia
result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = ac(),
    algorithm = admm(
        coupling = :hybrid,
        sectioning = :tie_line,
        rho = 1e3,
        adaptive_rho = true,
    ),
)
```

If a run oscillates, make the retuning more cautious: widen the trigger ratios so it only reacts to a clear imbalance, shrink the steps so each change is small, and set a floor so `rho` cannot collapse toward zero and drop the consensus penalty entirely.

```julia
result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = ac(),
    algorithm = admm(
        coupling = :hybrid,
        rho = 1e3,
        adaptive_rho = true,
        mu_incr = 50.0,      # only react to a large primal/dual imbalance
        mu_decr = 20.0,
        tau_incr = 1.1,      # small, cautious steps up
        tau_decr = 0.7,      # and down
        rho_floor = 1e2,     # never let rho fall below 100
        adapt_warmup = 10,   # ignore the first 10 iterations
    ),
)
```

## Zone warm-start

By default every zone starts the first iteration cold, from `z = 0` and `y = 0`. With `zone_warmstart = true` (the default is `false`) Lirion first solves each zone once on its own base objective — no penalty term, no dual term — and seeds the initial consensus from those solutions. A starting point that already respects each zone's own economics can shave iterations off the run.

It is opt-in because it is not free: it adds one full round of zone solves before the loop begins, and whether that pays for itself depends on the case. If any zone's warm-start solve fails, Lirion prints a warning and falls back to the cold start, so a failed warm-start costs time, never the result. You can also turn it on without editing the call by setting the `LIRION_ZONE_WARMSTART` environment variable to `1`, `true`, `yes`, or `on`.

```julia
result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = ac(),
    algorithm = admm(
        coupling = :hybrid,
        sectioning = :tie_line,
        rho = 1e3,
        zone_warmstart = true,
    ),
)
```

## Stopping

The run stops when the primal and dual residuals both fall below their tolerances, or when it reaches `max_iter`, whichever comes first. `max_iter` defaults to `350`.

Those tolerances follow the standard absolute-and-relative stopping criterion of Boyd et al. [Boyd2011](@cite): each one combines an absolute floor with a relative part that scales with the problem size and the norm of the current iterate, so the test stays meaningful across cases of different size. You set it through `eps_abs` and `eps_rel` (both default to `1e-7`); there is no single raw-residual threshold to pass. Loosen them to stop sooner with a looser consensus; tighten them to demand closer agreement. The full formula and its derivation live on the [ADMM theory page](../../teoria/distribuido/admm.md).

```julia
result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = dc(),
    algorithm = admm(
        coupling = :power,
        sectioning = :dummy_bus,
        rho = 1e3,
        max_iter = 500,
        eps_abs = 1e-6,
        eps_rel = 1e-6,
    ),
)
```

If the run reaches `max_iter` first, it returns the last iterate. To check how far it got and whether the residuals actually converged, read `result.iter_done` and the convergence history described in the Decentralized section of [Results](@ref).
