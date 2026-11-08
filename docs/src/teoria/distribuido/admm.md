# ADMM

The alternating direction method of multipliers (**ADMM**) is the first of Lirion's two distributed algorithms. It solves the consensus problem set up on the [sectioning page](seccionamento_acoplamento.md#Consistency-constraints) by alternating three cheap steps: every region minimizes its own objective plus a quadratic penalty on its boundary copies, the copies of each shared quantity are averaged, and a dual price per copy is raised by whatever mismatch remains. It asks the regions for nothing but their solutions. Nothing else crosses the boundary, no gradient and no Hessian, which is why it runs on every model Lirion ships, the conic SOCWR included. [ALADIN](aladin.md) solves the same regional subproblems but coordinates them with second-order information and takes far fewer iterations near the solution; it needs that information from the subproblem, which the conic SOCWR supplies only in part, so under ALADIN the conic model stalls near its solver's accuracy. The trade is iteration count against generality. ADMM is slower to converge and indifferent to what the subproblem is. The form implemented is the general consensus ADMM of Boyd et al. [Boyd2011](@cite), with one penalty per coupled quantity.

!!! note "What ADMM guarantees"
    ADMM drives the copies of every shared quantity to a common value up to a tolerance; the copies agree to ``\epsilon``, never exactly. For a convex model (DC, SOCWR, SOCWRNL) the consensus problem is equivalent to the centralized one and ADMM converges to its optimum [Boyd2011](@cite), so the distributed cost meets the centralized oracle. For the non-convex AC model the same iteration is a local method with no convergence guarantee. It finds a point where the regions agree and each is locally optimal, which is what the centralized AC solve offers too.

## Derivation

### Step 1: the consensus problem

Partition the buses into regions ``r = 1, \dots, R`` and let region ``r`` own the variables ``w^{(r)}`` of its local OPF, with cost ``f_r`` and feasible set ``\mathcal{X}_r`` (its own balance, flow, and bound constraints, ghosts included). Write ``x_s^{(r)}`` for region ``r``'s copy of the shared boundary quantity ``s``, a variable or an expression of ``w^{(r)}``, and ``\mathcal{C}_s`` for the set of regions that hold a copy of ``s``. The distributed OPF is the sum of the regional problems tied by the consistency constraints,

```math
\min_{w,\,z} \;\sum_{r=1}^{R} f_r\!\left(w^{(r)}\right)
\quad\text{s.t.}\quad
w^{(r)} \in \mathcal{X}_r, \qquad
x_s^{(r)} = z_s \;\; \text{for every } s \text{ and every } r \in \mathcal{C}_s .
```

The only thing that couples the regions is the last line, and it couples each copy to a single consensus value ``z_s`` rather than to the other copies. The motivation is privacy, not speed, as the [sectioning page](seccionamento_acoplamento.md) explains; ``z`` is the one object the parties agree to share.

### Step 2: the augmented Lagrangian

Dual decomposition would attach a multiplier ``y_s^{(r)}`` to each consistency constraint and let the regions minimize the Lagrangian separately. That alone converges slowly and only under strong assumptions. The method of multipliers adds a quadratic penalty with weight ``\rho > 0`` on the same constraint,

```math
L_\rho(w, z, y) = \sum_{r} f_r\!\left(w^{(r)}\right)
+ \sum_{s}\sum_{r \in \mathcal{C}_s}
\left[\, y_s^{(r)}\left(x_s^{(r)} - z_s\right)
+ \frac{\rho}{2}\left(x_s^{(r)} - z_s\right)^2 \right] .
```

The penalty makes the problem strongly convex in the mismatch and robust to a poor dual estimate. It also destroys separability, because ``z_s`` now sits inside every region's square. ADMM recovers it by minimizing over ``w`` and ``z`` in alternation instead of jointly, which is the whole idea of the method [Boyd2011](@cite).

!!! note "Justification"
    Lirion keeps the **unscaled** form above, with ``y`` the true multiplier of the consistency constraint. Boyd also writes a scaled form with ``u = y/\rho``, which is shorter. The unscaled form pays with a quadratic term that has to be rebuilt in each region's objective whenever ``\rho`` changes, and gains that the duals survive a change of ``\rho`` untouched, which matters for the adaptive penalty of Step 7.

### Step 3: the region subproblem

Fix ``z`` and ``y`` at their current values ``z^k`` and ``y^k``. The augmented Lagrangian then splits by region, and region ``r`` minimizes its own share in parallel with the others,

```math
w^{(r),k+1} = \arg\min_{w^{(r)} \in \mathcal{X}_r}
\; f_r\!\left(w^{(r)}\right)
+ \sum_{s :\, r \in \mathcal{C}_s}
\left[ \frac{\rho}{2}\left(x_s^{(r)}\right)^2
+ \left(y_s^{(r),k} - \rho\, z_s^k\right) x_s^{(r)} \right] .
```

This is the penalty of Step 2 expanded, with the constant ``\tfrac{\rho}{2}(z_s^k)^2 - y_s^{(r),k} z_s^k`` dropped. The quadratic part depends only on ``\rho`` and is built once; the linear part carries the current consensus and price and is the only thing rewritten between iterations. Each region solves a full OPF of its own model with its own solver (Ipopt for AC and SOCWRNL, OSQP for DC, Clarabel for the conic SOCWR). A solve that does not return an acceptable status aborts the run rather than feed a bad point to the next step.

### Step 4: consensus by averaging

Fix the new copies and minimize the augmented Lagrangian over ``z``. Each ``z_s`` appears in the squares of its copies only, so the minimizer is the mean of the copies shifted by the mean of their prices,

```math
z_s^{k+1} = \frac{1}{|\mathcal{C}_s|}\sum_{r \in \mathcal{C}_s}
\left( x_s^{(r),k+1} + \frac{1}{\rho}\, y_s^{(r),k} \right)
= \frac{1}{|\mathcal{C}_s|}\sum_{r \in \mathcal{C}_s} x_s^{(r),k+1} .
```

The second equality holds because the prices of one shared quantity always sum to zero. They start at zero, and the dual step below adds ``\rho`` times a set of deviations from their own mean, whose sum is zero. Lirion therefore averages the copies directly. The consensus step knows nothing of the network. It is an average, not a power flow, and that is the sense in which ADMM is first-order: boundary information travels one hop per iteration.

Two shared quantities need a different average. The AC angle is averaged as a circular mean, ``\bar\theta = \operatorname{atan2}(\sum \sin\theta^{(r)}, \sum \cos\theta^{(r)})``, so copies on either side of ``\pm\pi`` do not average to the wrong branch. The split-bus injections are not copies of one value but terms of a sum that must vanish; their consensus is the projection of the copies onto the zero-sum hyperplane, each copy minus the group mean, which is what the general consensus form prescribes when the agreement set is a subspace rather than a point.

### Step 5: dual ascent

With the new copies and consensus in hand, every price is raised by ``\rho`` times the mismatch of its copy,

```math
y_s^{(r),k+1} = y_s^{(r),k} + \rho\left( x_s^{(r),k+1} - z_s^{k+1} \right) .
```

A copy above the consensus accumulates a positive price and is pushed down in the next subproblem through the linear term of Step 3; a copy below is pushed up. On a two-region tie-line cut with a single shared angle ``\theta_h``, the two prices are equal and opposite at every iteration, ``y_h^{(1)} = -y_h^{(2)}``, so the dual is the price one region pays the other for the boundary angle. The step length is ``\rho`` itself. With that choice the consensus step satisfies its own dual feasibility condition exactly, and whatever dual infeasibility remains sits on the region side, where the dual residual of the next step measures it [Boyd2011](@cite).

### Step 6: residuals and the stopping rule

Two quantities measure progress. The **primal residual** is the stacked mismatch between copies and consensus; the **dual residual** is ``\rho`` times the stacked change of the consensus between iterations, replicated once per copy,

```math
r^{k+1} = \Big\| \big( x_s^{(r),k+1} - z_s^{k+1} \big)_{s,\, r \in \mathcal{C}_s} \Big\|_2 , \qquad
s^{k+1} = \rho\, \Big\| \big( z_s^{k+1} - z_s^{k} \big)_{s,\, r \in \mathcal{C}_s} \Big\|_2 .
```

The first is the violation of the consistency constraint. The second is the violation of dual feasibility, the stationarity condition of the subproblem evaluated at the new consensus, and it vanishes only when the consensus stops moving. Both must be small; a run with tiny ``r`` and large ``s`` has copies that agree on a value still sliding downhill.

Small is relative to the problem. Following Boyd, each tolerance combines an absolute floor scaled by the number ``n`` of copies with a relative part scaled by the size of the iterate,

```math
\epsilon^{\mathrm{pri}} = \sqrt{n}\,\epsilon^{\mathrm{abs}} + \epsilon^{\mathrm{rel}} \max\!\big(\|x\|_2, \|z\|_2\big), \qquad
\epsilon^{\mathrm{dual}} = \sqrt{n}\,\epsilon^{\mathrm{abs}} + \epsilon^{\mathrm{rel}}\,\|y\|_2 ,
```

and the iteration stops when ``r^{k+1} < \epsilon^{\mathrm{pri}}`` and ``s^{k+1} < \epsilon^{\mathrm{dual}}`` together, or at `max_iter` (default 350), whichever comes first. Lirion sets both ``\epsilon^{\mathrm{abs}}`` and ``\epsilon^{\mathrm{rel}}`` to ``10^{-7}`` by default. That is a demanding target. The residuals cannot fall below the accuracy of the regional solvers themselves, so when a solver's own tolerance is looser than the target (OSQP runs at ``10^{-6}``, the conic solver at ``10^{-4}``) the residuals floor above it and the run ends at `max_iter` with the last iterate and the full residual history. The knobs `eps_abs`, `eps_rel`, and `max_iter` are on the [settings page](../../usage/distributed/admm.md#Stopping).

### Step 7: one penalty per block, and its adaptation

A single ``\rho`` cannot serve a hybrid coupling. A boundary angle is a few hundredths of a radian and a boundary flow a few per unit, so one weight is either too stiff for the state or too loose for the power. Lirion gives every coupled quantity its own penalty, ``\rho_\theta``, ``\rho_v``, ``\rho_p``, ``\rho_q`` on AC and ``\rho_{v^2}``, ``\rho_r``, ``\rho_t`` on the W-space models, and computes the residuals of Step 6 per block as well, combining them in quadrature for the stop test. The price is more knobs; the generic `rho` fills them all with one value when that is good enough.

The value of ``\rho`` sets the balance between the two residuals. A large ``\rho`` forces agreement and lets the consensus crawl; a small one lets each region optimize freely and defers agreement. The right value is problem-dependent. Residual balancing adjusts it from the residuals themselves [Boyd2011](@cite),

```math
\rho_b \leftarrow
\begin{cases}
\tau^{\mathrm{incr}}\,\rho_b & \text{if } r_b > \mu^{\mathrm{incr}}\, s_b ,\\[2pt]
\max\!\big(\tau^{\mathrm{decr}}\,\rho_b,\ \rho^{\mathrm{floor}}\big) & \text{if } s_b > \mu^{\mathrm{decr}}\, r_b ,\\[2pt]
\rho_b & \text{otherwise},
\end{cases}
```

for each block ``b`` on its own residual pair. When the copies lag the consensus the penalty tightens; when the consensus is being dragged around the penalty loosens. The rule is off by default (`adaptive_rho = false`) and, when on, it waits out a warm-up of 5 iterations, runs after the stop test with the residuals just computed, and freezes once both residuals are within ``100\times`` their tolerances, so the final approach happens on a fixed penalty. Because the duals are unscaled multipliers, they are left alone when ``\rho`` changes; only the quadratic term of Step 3 is rebuilt. The defaults are ``\mu^{\mathrm{incr}} = 30``, ``\mu^{\mathrm{decr}} = 10``, ``\tau^{\mathrm{incr}} = 1.2``, ``\tau^{\mathrm{decr}} = 0.5``, and no floor; the [settings page](../../usage/distributed/admm.md#Adaptive-penalty) discusses when to change them.

---

## ADMM iteration

The constructor `admm()` builds the configuration for `solve`; the loop below lives in `src/distributed_algorithms/admm/`. Everything is in per unit, with ``\theta_{ref} = 0`` pinned in each region.

### Sets

- ``r \in \{1, \dots, R\}``: regions, from the `bus.zone` column.
- ``b \in \mathcal{B}``: coupling blocks, one per coupled quantity, as fixed by the model and the `coupling`: ``\{\theta\}`` on DC `:state`, ``\{\theta, p\}`` on DC `:hybrid`, ``\{\theta, |V|, p, q\}`` on AC `:hybrid`, and so on.
- ``s \in \mathcal{S}_b``: shared boundary quantities of block ``b`` (a shared bus state, an inter-zonal branch flow, a dummy-bus flow, or a split-bus injection).
- ``\mathcal{C}_s \subseteq \{1, \dots, R\}``: the regions holding a copy of ``s``; ``|\mathcal{C}_s| \ge 2``, and more than two at a bus shared by three or more regions.
- ``n_b = \sum_{s \in \mathcal{S}_b} |\mathcal{C}_s|``: the number of copies in block ``b``.
- ``k``: the iteration counter.

### Parameters

- ``\rho_b > 0``: penalty of block ``b`` (`rho`, or `rho_theta`, `rho_v`, `rho_p`, `rho_q`, `rho_v2`, `rho_r`, `rho_t`).
- ``\epsilon^{\mathrm{abs}}, \epsilon^{\mathrm{rel}} \ge 0``: stopping tolerances (`eps_abs`, `eps_rel`, both ``10^{-7}``).
- ``K``: iteration cap (`max_iter`, 350).
- ``\mu^{\mathrm{incr}}, \mu^{\mathrm{decr}}, \tau^{\mathrm{incr}}, \tau^{\mathrm{decr}}, \rho^{\mathrm{floor}}``, warm-up ``k_w`` and freeze factor ``\phi``: the adaptive-penalty knobs (``30, 10, 1.2, 0.5, 0, 5, 100``).

### Variables

- ``w^{(r)}``: region ``r``'s local OPF variables; ``f_r`` its cost; ``\mathcal{X}_r`` its constraints.
- ``x_s^{(r)}``: region ``r``'s copy of ``s``, a variable or an expression of ``w^{(r)}``.
- ``z_s``: the consensus value of ``s``.
- ``y_s^{(r)}``: the dual price of copy ``(s, r)``.
- ``x_b, z_b, y_b``: the block-``b`` vectors of copies, consensus values replicated once per copy, and prices, each in ``\mathbb{R}^{n_b}``.

### Initialization

```math
z_s^{0} = 0, \qquad y_s^{(r),0} = 0 \qquad \text{for all } s,\ r \in \mathcal{C}_s . \tag{1}
```

The run starts cold. With `zone_warmstart = true` each region is first solved once on ``f_r`` alone and ``z^0`` is taken from the average of those solutions, except for the angle, which is gauge-dependent and stays at zero; ``y^0`` is zero in either case.

### Region subproblem

For every region ``r``, in parallel,

```math
w^{(r),k+1} \in \arg\min_{w^{(r)} \in \mathcal{X}_r}
\; f_r\!\left(w^{(r)}\right)
+ \sum_{b \in \mathcal{B}} \sum_{\substack{s \in \mathcal{S}_b \\ r \in \mathcal{C}_s}}
\left[ \frac{\rho_b}{2}\left(x_s^{(r)}\right)^2
+ \left(y_s^{(r),k} - \rho_b\, z_s^k\right) x_s^{(r)} \right],
\qquad r = 1, \dots, R . \tag{2}
```

The copies ``x_s^{(r),k+1}`` are read off the solution. A region whose solver reports neither an optimal (or locally optimal) nor a reduced-accuracy status stops the run.

### Consensus

```math
z_s^{k+1} = \frac{1}{|\mathcal{C}_s|} \sum_{r \in \mathcal{C}_s} x_s^{(r),k+1},
\qquad s \in \mathcal{S}_b,\ b \in \mathcal{B} . \tag{3}
```

Two exceptions. For the AC angle the mean is circular,

```math
\bar\theta_s^{\,k+1} = \operatorname{atan2}\!\Big( \sum_{r \in \mathcal{C}_s} \sin\theta_s^{(r),k+1},\ \sum_{r \in \mathcal{C}_s} \cos\theta_s^{(r),k+1} \Big), \tag{4}
```

and for a split-bus injection group ``g`` with copies ``p_{\mathrm{inj},g}^{(r)}``, ``r \in \mathcal{C}_g``, the consensus is one value per copy, the copy minus the group mean,

```math
z_{g}^{(r),k+1} = p_{\mathrm{inj},g}^{(r),k+1} - \frac{1}{|\mathcal{C}_g|}\sum_{r' \in \mathcal{C}_g} p_{\mathrm{inj},g}^{(r'),k+1},
\qquad r \in \mathcal{C}_g , \tag{5}
```

which is the projection of the copies onto ``\sum_r p_{\mathrm{inj},g}^{(r)} = 0``. The previous consensus ``z^k`` is kept for the dual residual.

### Dual update

```math
y_s^{(r),k+1} = y_s^{(r),k} + \rho_b \left( x_s^{(r),k+1} - z_s^{k+1} \right),
\qquad s \in \mathcal{S}_b,\ r \in \mathcal{C}_s,\ b \in \mathcal{B} . \tag{6}
```

For the AC angle the difference inside the parentheses is wrapped to ``(-\pi, \pi]``; for a split-bus injection ``z_s`` is the per-copy value of (5), so the update adds ``\rho_b`` times the group mean to every copy's price.

### Residuals

Per block, with the block vectors of the previous section,

```math
r_b^{k+1} = \left\| x_b^{k+1} - z_b^{k+1} \right\|_2 , \qquad
s_b^{k+1} = \rho_b \left\| z_b^{k+1} - z_b^{k} \right\|_2 ,
\qquad b \in \mathcal{B}, \tag{7}
```

and over all blocks in quadrature,

```math
r^{k+1} = \Big( \sum_{b \in \mathcal{B}} \big(r_b^{k+1}\big)^2 \Big)^{1/2}, \qquad
s^{k+1} = \Big( \sum_{b \in \mathcal{B}} \big(s_b^{k+1}\big)^2 \Big)^{1/2} . \tag{8}
```

The angle entries of ``r_b`` and ``s_b`` on AC are wrapped differences, as in (6).

### Tolerances and stop

Per block,

```math
\epsilon_b^{\mathrm{pri}} = \sqrt{n_b}\,\epsilon^{\mathrm{abs}} + \epsilon^{\mathrm{rel}} \max\!\big( \|x_b^{k+1}\|_2,\ \|z_b^{k+1}\|_2 \big), \qquad
\epsilon_b^{\mathrm{dual}} = \sqrt{n_b}\,\epsilon^{\mathrm{abs}} + \epsilon^{\mathrm{rel}}\, \|y_b^{k+1}\|_2 ,
\qquad b \in \mathcal{B}, \tag{9}
```

combined over blocks in quadrature exactly as in (8) to give ``\epsilon^{\mathrm{pri}}`` and ``\epsilon^{\mathrm{dual}}``. The iteration terminates at the first ``k`` with

```math
r^{k+1} < \epsilon^{\mathrm{pri}} \quad\text{and}\quad s^{k+1} < \epsilon^{\mathrm{dual}}, \tag{10}
```

or when ``k + 1 = K``, in which case the last iterate is returned.

### Penalty update

Only with `adaptive_rho = true`, only for ``k + 1 > k_w``, and skipped whenever ``r^{k+1} < \phi\,\epsilon^{\mathrm{pri}}`` and ``s^{k+1} < \phi\,\epsilon^{\mathrm{dual}}`` (the freeze). For each block, on its own residual pair from (7),

```math
\rho_b \leftarrow
\begin{cases}
\tau^{\mathrm{incr}}\,\rho_b & r_b^{k+1} > \mu^{\mathrm{incr}}\, s_b^{k+1} ,\\[2pt]
\max\!\big(\tau^{\mathrm{decr}}\,\rho_b,\ \rho^{\mathrm{floor}}\big) & s_b^{k+1} > \mu^{\mathrm{decr}}\, r_b^{k+1} ,\\[2pt]
\rho_b & \text{otherwise},
\end{cases}
\qquad b \in \mathcal{B}, \tag{11}
```

with the floor ignored when ``\rho^{\mathrm{floor}} = 0``. A changed ``\rho_b`` rebuilds the quadratic term of (2) in every region before the next iteration; the prices ``y`` are not touched.

### The ADMM step

1. Initialize ``z^0 = 0``, ``y^0 = 0`` by (1), or warm-start ``z^0`` from the unpenalized regional solutions.
2. Solve the regional subproblems (2) in parallel and read off the copies.
3. Average the copies into the consensus (3), using (4) for the AC angle and (5) for split-bus injections.
4. Raise the prices by (6).
5. Compute the residuals (7)–(8) and the tolerances (9); stop if (10) holds or the cap ``K`` is reached.
6. If the adaptive penalty is on and not frozen, update each ``\rho_b`` by (11) and rebuild the penalty terms.
7. Return to 2.

## What it captures / what it loses

!!! info "Philosophy"
    **Captures:** the full regional model. Every region solves its own exact or relaxed OPF with its own solver, in parallel, and the coordination needs only the solutions, so the method runs unchanged on DC, AC, SOCWR, and SOCWRNL, on any sectioning, and on buses shared by any number of regions. Convergence to the centralized optimum is guaranteed for the convex models, and the per-block penalties let state and power converge at their own rates.

    **Loses:** curvature across the boundary. The coordination step is an average, so the regions learn about each other one hop per iteration and the rate depends on ``\rho`` and on the conditioning of the cut; a poor ``\rho`` costs hundreds of iterations, and the adaptive rule narrows that risk without removing it. The consensus is reached to a tolerance, never exactly, so the boundary quantities read from two regions differ by up to ``\epsilon``. On the non-convex AC model the method is local and offers no guarantee.

## Lirion vs. the original

!!! warning "Lirion vs. the original"
    Lirion implements the general-form consensus ADMM of [Boyd2011](@cite) in its unscaled form, with the following choices of its own.

    **One penalty per coupled quantity** (``\rho_\theta, \rho_v, \rho_p, \rho_q, \rho_{v^2}, \rho_r, \rho_t``) instead of one global ``\rho``, with the residuals and the adaptive rule evaluated per block. The tolerance bars are grouped by coupling class (state, flow, injection) and combined in quadrature, which bounds Boyd's single stacked test rather than reproducing it in value.

    **The consensus is the plain average** of the copies, which equals Boyd's price-shifted average because the prices start at zero and sum to zero afterwards. The AC angle uses a circular mean and wrapped differences, which Boyd's Euclidean setting does not cover, and a split-bus injection group uses a zero-sum projection in place of an average.

    **Residual balancing** follows Boyd's rule with asymmetric defaults, ``\mu^{\mathrm{incr}} = 30``, ``\mu^{\mathrm{decr}} = 10``, ``\tau^{\mathrm{incr}} = 1.2``, ``\tau^{\mathrm{decr}} = 0.5``, in place of his ``\mu = 10``, ``\tau = 2``, and adds a warm-up, a freeze near convergence, and an optional floor. It is off by default.

    **The stopping rule** is Boyd's absolute-plus-relative test with ``\epsilon^{\mathrm{abs}} = \epsilon^{\mathrm{rel}} = 10^{-7}``, tighter than the default accuracy of some regional solvers, so a run can end at `max_iter` with its residuals floored at the solver's precision rather than at the target.

    On the conic SOCWR under `:power` coupling, and only when the regional solver is Gurobi, a small Tikhonov term (``10^{-6}`` on ``v^2``, ``10^{-8}`` on the flows) is added to the subproblem for numerical stability; with the default conic solver it is absent.

## See also

- Neighbours: [Decomposition: sectioning and coupling](seccionamento_acoplamento.md) (the consistency constraints this page enforces, and what each model shares) · [ALADIN](aladin.md) (the second-order alternative, and why it needs more from the subproblem).
- Foundations: [Foundations](@ref).
- Usage: [ADMM settings](../../usage/distributed/admm.md) (penalties, adaptive rule, warm-start, stopping knobs) · [Results](@ref) (the residual history).
- Code: constructor `admm()`; source `src/distributed_algorithms/admm/` (`run.jl` for the loop, `core/consensus.jl`, `core/dual_update.jl`, `core/residuals.jl`, `core/tolerances.jl`, `core/rho_update.jl`).
- Original: [Boyd2011](@cite), the consensus form and the stopping and residual-balancing rules.
