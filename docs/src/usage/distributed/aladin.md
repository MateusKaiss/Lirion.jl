# ALADIN settings

ALADIN is the second distributed solver in Lirion. You reach it the same way as everything else, by passing `algorithm = aladin(...)` to `solve`. Each zone solves its own local OPF, a central step folds the local sensitivities into one coordination QP, and the two alternate until the boundary agrees. For how the network is cut and what the zones share, see [Running a distributed solve](@ref) and the [ALADIN](@ref) theory page.

`aladin()` runs on `dc()`, `ac()`, `socwr()`, and `socwrnl()`. With no extra arguments it couples on the shared state across tie lines and stops after at most 100 iterations.

```julia
using Lirion

result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = ac(),
    algorithm = aladin(),
)
```

That call is `aladin(coupling = :state, sectioning = :tie_line, max_iter = 100)` written out in full. The sections below cover each knob.

## Coupling and sectioning for ALADIN

`coupling` decides which boundary quantities the zones agree on, and `sectioning` decides how the network is split into zones. The pairs ALADIN accepts are:

| `coupling` | `:tie_line` | `:dummy_bus` | `:split_bus` |
| --- | --- | --- | --- |
| `:state` | yes | no | no |
| `:power` | yes | yes | no |
| `:hybrid` | yes | yes | yes |

The rule behind the table is short. `:dummy_bus` only means something when the zones couple on power, so it needs `:power` or `:hybrid`. `:split_bus` shares bus injections, which only the hybrid workspaces build, so it needs `:hybrid`. `:state` has nothing to agree on across a dummy or split bus, so it stays on `:tie_line`.

One combination is left out on purpose. With the conic `socwr()`, ALADIN does not support `:dummy_bus` under any coupling. The coordination QP on the conic model relaxes the dummy-bus link and the run would settle on the wrong point, so the builder stops you before it does. Use `admm()` if you need conic SOCWR with a dummy bus, keep ALADIN on `:tie_line` or `:split_bus`, or switch to `socwrnl()`, which carries the curvature the QP needs and is **not** subject to this restriction.

```julia
result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = ac(),
    algorithm = aladin(coupling = :hybrid, sectioning = :split_bus),
)
```

Anything outside the table throws an `ArgumentError` that lists what the model and coupling allow, so an invalid call fails up front instead of running quietly wrong.

## Iterations and the penalty

`max_iter` caps the number of coordination rounds and defaults to 100. ALADIN often converges in fewer rounds than ADMM, which is why its default sits below ADMM's 350, but a hard case can still want more.

`rho` is the **local** penalty each zone puts on drifting away from the last agreed boundary point while it solves its own subproblem. It is a single positive number and defaults to `1e3`. A larger `rho` holds the zones closer to the current consensus each round, which steadies a run that is bouncing around, at the cost of smaller steps when a bigger one would do. The *central* coordination step has a separate penalty, `mu`, covered under [The coordination step](@ref) below.

```julia
result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = ac(),
    algorithm = aladin(rho = 1e4, max_iter = 200),
)
```

ALADIN's `rho` is one scalar, not the per-block `rho_theta` / `rho_p` family you set on ADMM. Those names are an ADMM feature and are rejected here. For the per-block penalties, see [ADMM settings](@ref).

Set `adaptive = true` to let both penalties grow as the run goes. Each round multiplies `rho` by `rho_growth` (default `1.3`, capped at `rho_upper`, default `1e6`) and `mu` by `mu_growth` (default `2.0`, capped at `mu_upper`, default `1e6`). The growth is unconditional, not a reaction to the residuals, which lets a stubborn case keep tightening. It is off by default, and with it off `rho` and `mu` stay fixed.

```julia
result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = ac(),
    algorithm = aladin(adaptive = true, max_iter = 200),
)
```

## Proximal mode

`proximal_mode` sets what the local `rho` penalty is written on. The default, `:primal`, penalizes the boundary variables directly. `:pflow` instead penalizes the power flows across the boundary, usually the better choice when the zones couple on power.

Because `:pflow` works on flows, it needs a coupling that actually shares them, so it is only allowed with `:power` or `:hybrid`; asking for it under `:state` throws. Under `:pflow`, `sigma_pflow` and `sigma_qflow` scale the active- and reactive-flow penalty terms relative to `rho`. Both default to `10.0` and must be positive. Under the default `:primal` mode they do nothing, and on DC `sigma_qflow` is inert in any case, since DC carries no reactive power.

```julia
result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = ac(),
    algorithm = aladin(
        coupling = :hybrid,
        sectioning = :tie_line,
        proximal_mode = :pflow,
        sigma_pflow = 20.0,
        stopping_criterion = :hfd,
    ),
)
```

## Globalization

Each coordination round proposes a step, and `globalization` decides how much of it to take. The update mixes three pieces: the local move `y − x`, the coordinator's correction `Δy`, and the dual move toward the QP multipliers. The globalization sets how far the step goes along each one. It affects whether a hard case converges more than any other setting, so Lirion offers four strategies. The default is `:armijo_l2`.

- **`:armijo_l2`** (default) — an ℓ2 Armijo backtracking line search on the augmented-Lagrangian merit. It takes the local and dual moves and backtracks the coordinator correction until it gives sufficient decrease. It is the most forgiving strategy, and the only one available on every model.
- **`:full_step`** — take the whole step undamped (the plain [Houska2016](@cite) Algorithm 2 update). It is the fastest per round when the problem is well-behaved, but with no line search to catch a bad step it can diverge on a hard case.
- **`:alg3`** — an ℓ1 exact-penalty merit (Algorithm 3 of [Houska2016](@cite)) with the paper's three-branch step (full step, local-only, or a golden-section line search on the dual), which buys stronger global-convergence guarantees than a plain line search.
- **`:alg3_original`** — the same, but the dual branch uses a Wolfe line search with a gradient fallback, and the exact-penalty weight can be pinned with `alg3_lam_bar` instead of chosen automatically. This is the paper-literal variant, kept for reproducibility.

`:full_step`, `:alg3`, and `:alg3_original` are implemented for `dc()`, `ac()`, `socwr()`, and `socwrnl()`; on any other model only `:armijo_l2` is available. An unknown value throws an `ArgumentError` that lists the four.

```julia
# full, undamped steps, fast when the case is well-conditioned
result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = dc(),
    algorithm = aladin(globalization = :full_step),
)
```

### Step coefficients

`alpha1`, `alpha2`, `alpha3` (all default `1.0`, all must be positive) are the step lengths of the three pieces of the update: `alpha1` scales the local move `y − x`, `alpha2` the coordinator correction `Δy`, and `alpha3` the dual move. How much each matters depends on the globalization. `:full_step` takes all three as given. `:armijo_l2` sets `alpha2` by line search (your value is its starting step) and applies `alpha1` and `alpha3` as configured, though tuning them there is rarely useful, since the line search validates the step assuming `alpha1 = 1`. `:alg3` and `:alg3_original` set their own step lengths each round and ignore all three. Leave them at `1.0` unless you are deliberately hand-tuning a `:full_step` run.

`alg3_lam_bar` (default `0.0`) is the fixed weight on the ℓ1 constraint-violation term in the Algorithm 3 merit, and it only has an effect under `:alg3_original`. At the default `0.0`, Lirion chooses the weight automatically; a positive value pins it for the whole run instead. The automatic default is the right choice unless you are reproducing the paper's constant-weight setting.

## The coordination step

Where `rho` penalizes each zone locally, `mu` (default `1e4`) is the penalty the **central** coordination QP puts on the predicted boundary mismatch, the augmented-Lagrangian term `λ·s + ½·mu·s²` on the coupling slack `s`. A larger `mu` makes the coordinator insist harder on closing the mismatch in a single step. Like `rho`, it grows under `adaptive = true` (by `mu_growth`, capped at `mu_upper`); otherwise it stays at `1e4`.

`qp_optimizer` sets the solver behind that coordination QP and defaults to `Clarabel`. It is a **separate** slot from the zone solver: `optimizer` (passed to `solve`) changes how each zone's local OPF is solved, exactly as in a centralized run, while `qp_optimizer` changes only the central QP. You rarely need to touch either.

```julia
result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = ac(),
    algorithm = aladin(mu = 1e5, max_iter = 200),
)
```

To solve that QP with a different solver, pass a factory as `qp_optimizer`. Any solver that handles the convex QP works; here it is pointed at OSQP instead of the default Clarabel:

```julia
using JuMP, OSQP

result = solve("paper/cases/partitioned/case24_ieee_rts/std/k3/case.m";
    model = ac(),
    algorithm = aladin(qp_optimizer = OSQP.Optimizer),
)
```

## Stopping criterion

ALADIN stops only when **both** residuals of its KKT check fall under tolerance: the **consensus mismatch** (the boundary copies agree) and the **stationarity residual** `ρ‖Σ(y−x)‖₁`, the ℓ1 size of the proximal step, which measures how close the iterate is to a KKT point. Algorithm 2 of [Houska2016](@cite) requires both. Testing the stationarity residual, not the mismatch alone, is what sets this apart from a mismatch-only stop.

The stationarity bar is `n_prox · eps_abs`, the ℓ1 analogue of Boyd's dual bound. It scales with the problem's dimension and is the same under either criterion. `stopping_criterion` changes only how the **consensus** residual is measured:

- **`:boyd`** (the default) — the consensus bar is Boyd's own primal test, `√n · eps_abs + eps_rel · ‖·‖` ([Boyd2011](@cite), §3.3.1), so it also grows with the problem size (`√n`).
- **`:hfd`** — the ALADIN paper's consensus bar, a flat `eps_abs` measured in the max-norm, which does not grow with the network.

So `:boyd` borrows Boyd's primal test for the consensus side; the stationarity side is the ℓ1 bar either way, not Boyd's ℓ2 dual residual. Both `eps_abs` and `eps_rel` default to `1e-4`. Tighten them for a closer consensus, loosen them to stop sooner.

## Advanced scaling (model-specific)

Two more sets of knobs exist, and both are best left alone.

`sigma_theta_shared`, `sigma_pg`, and `sigma_other` scale the per-class weights of each zone's proximal term on **DC** — the shared boundary angle, the generator setpoints, and everything else, respectively (defaults `10.0`, `1.0`, `1.0`). They are wired into the DC workspaces only; on AC, SOCWR, and SOCWRNL the per-class weights are fixed inside the workspace, and passing a non-default value there throws rather than being silently ignored.

`active_set_tol` is fixed at `1e-6`. ALADIN detects active box bounds with that hard-coded tolerance, so the knob exists only to reject a non-default value loudly instead of pretending to honour it.

## Reading convergence

An ALADIN run returns the usual `LirionResult`. `result.status` and `result.iter_done` tell you whether it converged and in how many rounds, and `result.diagnostics` holds an `ALADINHistory` with the per-round residuals. The Decentralized section of [Results](@ref) lists every field.
