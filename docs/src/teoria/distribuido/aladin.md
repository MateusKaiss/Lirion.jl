# ALADIN

ALADIN, the Augmented Lagrangian Alternating Direction Inexact Newton method of
[Houska2016](@cite), is the second distributed solver in Lirion. Like ADMM it lets
each region solve its own OPF and reconciles the boundary copies afterwards, but
the reconciliation is different in kind. ADMM averages the copies and raises a
price, asking each region for nothing but its solution. ALADIN asks each region for its local **sensitivities**,
the gradient of its cost, an approximation of the Hessian of its Lagrangian, and
the Jacobian of its active constraints, and a central coordinator solves one
equality-constrained QP that moves all the regions at once. The step is
Newton-like, so close to the solution it takes far fewer iterations than ADMM.
The trade is information: the coordinator must see curvature and active sets, not
just boundary values, and a model that supplies only part of that curvature, the
conic SOCWR, is coordinated poorly, its step stalling near the conic solver's
accuracy. This page
derives the iteration; the knobs, `rho`, `mu`, the globalizations and the two
stopping criteria, are on the [ALADIN settings](../../usage/distributed/aladin.md)
page.

!!! note "What the ALADIN certificate guarantees"
    ALADIN stops only when two residuals are small at the same time: the
    **consensus mismatch** between the regional copies, and the **stationarity
    residual** ``\rho\lVert\Sigma(y-x)\rVert_1``, the size of the proximal step.
    [Houska2016](@cite) shows that the pair implies the local solutions satisfy
    the first-order KKT conditions of the full problem up to the tolerance, with
    the coordinator multiplier ``\lambda`` as the dual of the coupling. The stop
    therefore certifies an approximate KKT point, not merely agreement at the
    boundary. The ADMM stop also tests two residuals, but its second one is the
    motion of the consensus, ``\rho\lVert z^{k+1} - z^{k}\rVert``, a proxy for
    dual feasibility; it never evaluates stationarity at the regional solutions.

## Derivation

### The structured problem

Partition the network into regions ``r = 1, \dots, R`` as on the
[sectioning page](seccionamento_acoplamento.md), and let ``x^{(r)}`` collect the
variables region ``r`` owns, its ghost copies included. Each region has a cost
``f_r`` and a feasible set ``\mathcal{X}_r`` (its own balance, flow, bound, and
angle constraints). The consistency constraints couple the regions through their
copies, and nothing else does. Written as one problem,

```math
\min_{x^{(1)},\dots,x^{(R)}} \;\sum_{r=1}^{R} f_r\big(x^{(r)}\big)
\quad\text{s.t.}\quad
x^{(r)} \in \mathcal{X}_r, \qquad
\sum_{r=1}^{R} A_r\, x^{(r)} = 0 .
```

The matrix ``A_r`` picks out, with a sign, the boundary quantities region ``r``
holds a copy of. For a tie-line cut with one shared angle ``\theta_h`` between
two regions, ``A_1 x^{(1)} = \theta_h^{(1)}`` and ``A_2 x^{(2)} = -\theta_h^{(2)}``,
so the coupling row reads ``\theta_h^{(1)} - \theta_h^{(2)} = 0``. Under `:power`
the row is a difference of flows, which in the DC model are themselves linear in
the angles; under `:split_bus` it is a sum of injections. In every case the
coupling is a set of affine rows, which is the form [Houska2016](@cite) assumes.
Every coupling row carries one multiplier ``\lambda_j``, the price of its
mismatch, and ``\lambda`` is shared by all the regions that enter the row.

### Step 1: the local proximal subproblem

Each region cannot see the others, so it is given two things to stand in for
them: the current price vector ``\lambda`` and its own previous iterate
``x^{(r)}``, which acts as a **proximal centre**. The region then solves

```math
y^{(r)} \in \arg\min_{y \in \mathcal{X}_r}\;
f_r(y) + \lambda^{\top} A_r\, y
+ \frac{\rho}{2}\,\lVert y - x^{(r)} \rVert_{\Sigma_r}^{2} .
```

The linear term charges the region the price of its boundary copies. The
quadratic term, with ``\lVert v \rVert_{\Sigma}^2 = v^{\top}\Sigma v``, pulls the
solution toward the last agreed point, with the penalty ``\rho`` (`rho`, default
``10^3``) setting how hard. The weight matrix ``\Sigma_r`` is diagonal and groups
the variables into classes: on DC the shared boundary angles carry a weight of
10 and the generator setpoints and everything else a weight of 1, so the proximal
term holds the boundary closer than the interior. The subproblem keeps all of the
region's own constraints, so ``y^{(r)}`` is always locally feasible; what it
lacks is agreement with the neighbours.

!!! note "Why a proximal term and not an ADMM penalty"
    ADMM penalizes the distance to the *consensus* value. ALADIN penalizes the
    distance to the region's *own* previous point and leaves the consensus to
    the coordinator. The two coincide for a particular choice of the matrices
    (``\Sigma_r = A_r^{\top}A_r``, see section 5 of [Houska2016](@cite)), which is
    how the paper shows that one ALADIN iteration costs the same as one ADMM
    iteration when the coordinator is given the ADMM-equivalent curvature. The
    gain comes from giving it better curvature.

On the first iteration the proximal centre is the solver's start values;
afterwards Lirion resets it to the first local solution, ``x^{(r)} \leftarrow y^{(r)}``,
so the proximal step is zero at ``k = 1`` and convergence is only certified from
the second iteration on.

### Step 2: the sensitivities each region exports

After its solve, region ``r`` reports three objects about the point ``y^{(r)}``.

```math
g_r = \nabla f_r\big(y^{(r)}\big), \qquad
B_r \approx \nabla^2\Big[f_r + \kappa_r^{\top} h_r\Big]\big(y^{(r)}\big), \qquad
C_r = \frac{\partial h_r^{\mathrm{act}}}{\partial y}\Big|_{y^{(r)}} .
```

The vector ``g_r`` is the cost gradient. The matrix ``B_r`` is a symmetric,
positive-semidefinite approximation of the Hessian of the region's Lagrangian,
with ``\kappa_r`` the local multipliers of the constraints ``h_r \le 0``. The
matrix ``C_r`` is the Jacobian of the constraints that are **active** at
``y^{(r)}``: every equality row, every inequality whose multiplier is non-zero
beyond a tolerance, and every variable sitting on a bound (within ``10^{-6}`` of
it). Linearly dependent rows are dropped so the active Jacobian has full row rank.
This is the inexact-Newton half of the name: ``B_r`` and ``C_r`` need only
approximate the true matrices for the method to converge, and how well they do
sets the local rate.

What ``B_r`` contains depends on the model. On DC the cost is the only curvature,
so ``B_r`` is the diagonal of quadratic cost coefficients on the generator columns
with a small floor elsewhere. On AC and on SOCWRNL the region is a smooth NLP, and
Lirion assembles the exact Lagrangian Hessian ``\nabla^2 f_r + \sum_j \kappa_j
\nabla^2 h_j`` from the solver's evaluator, adds the proximal curvature
``\rho\Sigma_r``, and shifts the diagonal by the smallest ``\delta I`` that makes
the result positive semidefinite. On the conic SOCWR the cones have no smooth
Hessian to read, so ``B_r`` falls back to the DC-style diagonal plus
``\rho\Sigma_r``, and the cones contribute only their active-face normals to ``C_r``. That
missing curvature is the limitation described on the
[sectioning page](seccionamento_acoplamento.md#SOCWR,-and-its-limitations), and
the reason SOCWRNL exists.

### Step 3: the coordinator QP

The coordinator now has, for every region, a locally feasible point and a
quadratic model of the region's behaviour around it. It looks for the correction
``\Delta y^{(r)}`` of each point that closes the coupling while keeping the summed
quadratic models as small as possible, which is a QP,

```math
\begin{aligned}
\min_{\Delta y,\, s}\quad
& \sum_{r=1}^{R} \Big[\tfrac12\, \Delta y^{(r)\top} B_r\, \Delta y^{(r)}
+ g_r^{\top} \Delta y^{(r)}\Big]
+ \lambda^{\top} s + \tfrac{\mu}{2}\, \lVert s \rVert_2^2 \\
\text{s.t.}\quad
& \sum_{r=1}^{R} A_r \big(y^{(r)} + \Delta y^{(r)}\big) = s
\quad\big|\ \lambda_{\mathrm{QP}}, \\
& C_r\, \Delta y^{(r)} = 0, \qquad r = 1, \dots, R .
\end{aligned}
```

The first constraint is the coupling, written on the corrected points and
relaxed by a **slack** ``s``: the coordinator may leave a residual mismatch, but
pays ``\lambda^{\top}s + \tfrac{\mu}{2}\lVert s\rVert^2`` for it. The penalty
``\mu`` (`mu`, default ``10^4``) is the coordinator's own augmented-Lagrangian
weight, separate from the regions' ``\rho``; as ``\mu \to \infty`` the slack
closes and the row becomes a hard equality. The second constraint keeps each
correction inside the tangent space of the region's active set, so the QP does
not propose moving a generator past a bound that the local solve already found
binding; active bounds enter as ``\Delta y_k = 0``. The multiplier
``\lambda_{\mathrm{QP}}`` of the coupling rows is the coordinator's new estimate
of the price. In the two-region tie-line example the whole coupling block is the
single row ``\theta_h^{(1)} + \Delta\theta_h^{(1)} - \theta_h^{(2)} -
\Delta\theta_h^{(2)} = s``, and ``\lambda_{\mathrm{QP}}`` is one number.

The QP is convex when every ``B_r`` is positive semidefinite, which the shift
``\delta I`` guarantees, and it is structured: a block-diagonal Hessian tied only
by the coupling rows, with the corrections and the slack as unknowns.
Lirion hands it to an interior-point solver (`Clarabel` by default, `qp_optimizer`
to change it), because the QP becomes badly conditioned on large networks and a
first-order QP solver stops reaching tolerance there.

!!! note "Why the slack"
    Without ``s`` the coupling would be a hard equality and ``\lambda_{\mathrm{QP}}``
    its exact multiplier; the paper's section 4 shows the two forms are the same
    QP, with ``\lambda_{\mathrm{QP}} = \lambda + \mu\, s`` at the optimum. The
    slack version is better conditioned when ``\mu`` is large, because the large
    entries land on the diagonal, and it is always feasible: ``\Delta y = 0``,
    ``s = \sum_r A_r y^{(r)}`` is a feasible point whatever the active sets.

### Step 4: the primal and dual update

Each region now holds three points: its previous iterate ``x^{(r)}``, its local
solution ``y^{(r)}``, and the coordinator's correction ``\Delta y^{(r)}``. The
update blends them, with a separate step length on each move,

```math
x^{(r)} \leftarrow x^{(r)} + \alpha_1\big(y^{(r)} - x^{(r)}\big) + \alpha_2\, \Delta y^{(r)},
\qquad
\lambda \leftarrow \lambda + \alpha_3\big(\lambda_{\mathrm{QP}} - \lambda\big).
```

With ``\alpha_1 = \alpha_2 = \alpha_3 = 1`` this is the plain update of
Algorithm 2 in [Houska2016](@cite): the new iterate is the corrected local
solution ``y^{(r)} + \Delta y^{(r)}``, and the price jumps to the coordinator's
estimate. The new ``x^{(r)}`` becomes the proximal centre of the next local
solve, and the new ``\lambda`` its price, so the loop closes.

### Step 5: globalization

The full step is what gives the fast local rate, but far from the solution it
can overshoot, and Newton-type methods need a safeguard against that. ALADIN's
safeguard is the choice of ``\alpha_1, \alpha_2, \alpha_3``. Lirion's default,
`:armijo_l2`, keeps ``\alpha_1`` and ``\alpha_3`` at their configured values and
backtracks ``\alpha_2`` on the coordinator's own augmented Lagrangian,

```math
\phi(\alpha) = \sum_{r} f_r\big(y^{(r)} + \alpha\,\Delta y^{(r)}\big)
+ \lambda^{\top} c(\alpha) + \tfrac{\mu}{2}\, \lVert c(\alpha) \rVert_2^2 ,
\qquad
c(\alpha) = \sum_r A_r\big(y^{(r)} + \alpha\,\Delta y^{(r)}\big),
```

halving ``\alpha`` from 1 until ``\phi(\alpha) \le \phi(0) + 10^{-4}\,\alpha\,\phi'(0)``.
On DC and SOCWR a search that reaches the floor ``10^{-3}`` returns it; on AC it
instead falls back to Algorithm 3's dual-only step (``\alpha_1 = \alpha_2 = 0``,
a golden-section search on ``\alpha_3``). The merit is the one the QP was built
to decrease, so the test is cheap and local. The paper's own strategy, Algorithm 3, is an
exact ``\ell_1`` penalty with three branches (full step, local step only, or a
dual-only line search on ``\alpha_3``) and carries the global convergence proof;
Lirion offers it as `:alg3` and `:alg3_original`, and the undamped `:full_step`
for well-behaved cases. Which one to use, and what the step coefficients do under
each, is on the [ALADIN settings](../../usage/distributed/aladin.md#Globalization)
page.

### Step 6: the stopping rule

The iteration stops when the local solutions agree at the boundary *and* the
proximal step has died out. Both are measured on the fresh local solutions
``y^{(r)}``, before the update,

```math
r^{\mathrm{cons}} = \Big\lVert \sum_{r} A_r\, y^{(r)} \Big\rVert ,
\qquad
r^{\mathrm{stat}} = \rho\, \Big\lVert \big(\Sigma_r\,(y^{(r)} - x^{(r)})\big)_{r=1}^{R} \Big\rVert_1 ,
```

and the run is certified when ``r^{\mathrm{cons}} < \epsilon^{\mathrm{cons}}``
and ``r^{\mathrm{stat}} < \epsilon^{\mathrm{stat}}`` on the same iteration.
The first residual is the consensus mismatch. The second is the gradient of the
proximal term at ``y^{(r)}``, which by the subproblem's own optimality condition
equals ``\nabla f_r + A_r^{\top}\lambda + C_r^{*\top}\kappa_r`` up to sign, exact
only for an exact local solve and otherwise up to the regional solver's own KKT
tolerance: a small ``r^{\mathrm{stat}}`` means the current price ``\lambda``
already makes each region's solution stationary for the coupled problem. Dropping it would
leave a method that stops as soon as the copies agree, whatever the price, and
that is not ALADIN's stop; step 2 of Algorithm 2 in [Houska2016](@cite) tests
both.

The stationarity bar is ``\epsilon^{\mathrm{stat}} = n_{\mathrm{prox}}\,\epsilon^{\mathrm{abs}}``,
with ``n_{\mathrm{prox}}`` the number of variables carrying a proximal weight,
the ``\ell_1`` analogue of Boyd's dimension-scaled bound. The consensus bar
depends on `stopping_criterion`. Under the default `:boyd`,
``r^{\mathrm{cons}}`` is the 2-norm and the bar is Boyd's primal test
``\sqrt{n}\,\epsilon^{\mathrm{abs}} + \epsilon^{\mathrm{rel}} \max(\lVert y_I\rVert, \lVert y_J\rVert)``
[Boyd2011](@cite), with ``n`` the number of coupling rows and ``y_I, y_J`` the
two sides of those rows (the copies on one side of each cut and on the other);
under `:hfd` the bar is a flat ``\epsilon^{\mathrm{abs}}`` as in the paper,
measured in the max-norm where the paper uses the 1-norm. Both tolerances default
to ``10^{-4}``. The paper notes that the test is not scale invariant, and neither
is Lirion's. A third quantity, the dual step ``\mu\lVert\lambda^{+} - \lambda\rVert``,
is recorded for diagnostics but is not part of the stop test.

---

## ALADIN iteration

The constructor `aladin()` builds one local model per region, one coordinator
QP, and the loop below. It runs on `dc()`, `ac()`, `socwr()`, and `socwrnl()`.

### Sets

- ``\mathcal{R} = \{1, \dots, R\}``: regions (the `bus.zone` values).
- ``\mathcal{K}``: coupling rows, one per shared quantity (a boundary-bus state
  per copy pair, a boundary flow, a split-bus injection group), as fixed by the
  `coupling` and `sectioning` choice.
- ``\mathcal{A}_r \subseteq \{1, \dots, n_r\}``: columns of region ``r`` on an
  active bound at the current local solution.

### Parameters

- ``A_r \in \mathbb{R}^{|\mathcal{K}| \times n_r}``: signed selection of region
  ``r``'s boundary copies into the coupling rows.
- ``\Sigma_r \succ 0``: diagonal proximal weights by variable class
  (`sigma_theta_shared`, `sigma_pg`, `sigma_other` on DC; fixed per class on the
  other models).
- ``\rho > 0``: local proximal penalty (`rho`, default ``10^3``).
- ``\mu > 0``: coordinator slack penalty (`mu`, default ``10^4``).
- ``\epsilon^{\mathrm{abs}}, \epsilon^{\mathrm{rel}}``: tolerances (`eps_abs`,
  `eps_rel`, both ``10^{-4}``); ``\alpha_1, \alpha_2, \alpha_3``: step
  coefficients (all ``1``).

### Variables

- ``x^{(r)} \in \mathbb{R}^{n_r}``: current iterate of region ``r``, the proximal
  centre.
- ``y^{(r)} \in \mathbb{R}^{n_r}``: local solution of region ``r``.
- ``\Delta y^{(r)} \in \mathbb{R}^{n_r}``: coordinator correction for region ``r``.
- ``s \in \mathbb{R}^{|\mathcal{K}|}``: coupling slack.
- ``\lambda, \lambda_{\mathrm{QP}} \in \mathbb{R}^{|\mathcal{K}|}``: coupling
  multipliers, current and coordinator-proposed.

### Local subproblem

```math
y^{(r)} \in \arg\min_{y \in \mathcal{X}_r}\;
f_r(y) + \lambda^{\top} A_r\, y
+ \frac{\rho}{2}\,(y - x^{(r)})^{\top} \Sigma_r\, (y - x^{(r)}),
\qquad r \in \mathcal{R}. \tag{1}
```

### Sensitivities

```math
g_r = \nabla f_r\big(y^{(r)}\big), \qquad
B_r \approx \nabla^2\big[f_r + \kappa_r^{\top} h_r\big]\big(y^{(r)}\big) \succeq 0, \qquad
C_r = \nabla h_r^{\mathrm{act}}\big(y^{(r)}\big)^{\top},
\qquad r \in \mathcal{R}. \tag{2}
```

### Coordinator QP

```math
\min_{\Delta y,\, s}\;
\sum_{r \in \mathcal{R}} \Big[\tfrac12\, \Delta y^{(r)\top} B_r\, \Delta y^{(r)}
+ g_r^{\top} \Delta y^{(r)}\Big]
+ \lambda^{\top} s + \tfrac{\mu}{2}\, s^{\top} s \tag{3}
```

```math
\sum_{r \in \mathcal{R}} A_r \big(y^{(r)} + \Delta y^{(r)}\big) = s
\quad\big|\ \lambda_{\mathrm{QP}} \tag{4}
```

```math
C_r\, \Delta y^{(r)} = 0, \qquad
\Delta y^{(r)}_i = 0,\ i \in \mathcal{A}_r, \qquad r \in \mathcal{R}. \tag{5}
```

### Step update

```math
x^{(r)} \leftarrow x^{(r)} + \alpha_1\big(y^{(r)} - x^{(r)}\big) + \alpha_2\, \Delta y^{(r)},
\qquad r \in \mathcal{R}, \tag{6}
```

```math
\lambda \leftarrow \lambda + \alpha_3\big(\lambda_{\mathrm{QP}} - \lambda\big). \tag{7}
```

### Residuals and tolerances

```math
r^{\mathrm{cons}} = \Big\lVert \sum_{r \in \mathcal{R}} A_r\, y^{(r)} \Big\rVert_{2\ \text{or}\ \infty},
\qquad
r^{\mathrm{stat}} = \rho \sum_{r \in \mathcal{R}} \big\lVert \Sigma_r\,(y^{(r)} - x^{(r)}) \big\rVert_1 . \tag{8}
```

```math
\epsilon^{\mathrm{cons}} =
\begin{cases}
\sqrt{|\mathcal{K}|}\,\epsilon^{\mathrm{abs}} + \epsilon^{\mathrm{rel}} \max\big(\lVert y_I \rVert_2, \lVert y_J \rVert_2\big) & \texttt{:boyd} \\[2pt]
\epsilon^{\mathrm{abs}} & \texttt{:hfd}
\end{cases},
\qquad
\epsilon^{\mathrm{stat}} = n_{\mathrm{prox}}\,\epsilon^{\mathrm{abs}}. \tag{9}
```

Under `:boyd` the norm in (8) is the 2-norm; under `:hfd` it is the max-norm.
``y_I`` and ``y_J`` stack the copies on the two sides of the coupling rows, and
``n_{\mathrm{prox}}`` is the number of columns with a non-zero entry in some
``\Sigma_r``.

### The ALADIN step

Starting from ``\lambda = 0`` and the first local solve as the proximal centre,
repeat for ``k = 1, 2, \dots`` until `max_iter` (default 100):

1. Every region solves (1) in parallel at the current ``(x^{(r)}, \lambda)`` and
   exports (2).
2. Solve the QP (3)–(5) for ``\Delta y``, ``s``, and ``\lambda_{\mathrm{QP}}``.
3. Compute (8) and (9). If ``r^{\mathrm{cons}} < \epsilon^{\mathrm{cons}}`` and
   ``r^{\mathrm{stat}} < \epsilon^{\mathrm{stat}}``, and ``k \ge 2``, stop:
   ``y`` is the solution and ``\lambda`` its coupling dual.
4. Choose ``\alpha_1, \alpha_2, \alpha_3`` by the configured `globalization`
   (default: Armijo backtracking on ``\alpha_2``).
5. Apply (6) and (7). With `adaptive = true`, grow ``\rho`` and ``\mu`` by their
   factors, capped at their upper bounds. Return to 1.

## What it captures / what it loses

!!! info "Philosophy"
    **Captures:** the *curvature* of each region and its *active set*. Because
    the coordinator sees ``B_r`` and ``C_r``, its step is a Newton step on the
    coupled KKT system. Houska's quadratic rate needs exact Hessians, full steps,
    and ``\mu \to \infty`` at a regular KKT point; with the fixed approximations
    and the proximal term Lirion ships, the guaranteed rate is linear, in
    practice much faster than ADMM's. The stop certifies stationarity as well as
    agreement, so a converged run is an approximate KKT point of the whole problem.

    **Loses:** *locality of information*. ADMM exchanges boundary values only;
    ALADIN exports a gradient, a Hessian approximation, and an active Jacobian
    from every region to one coordinator, and solves a central QP each
    iteration. It also needs that curvature to exist: a region whose model is
    conic has no Lagrangian Hessian to export, which is why the conic SOCWR
    converges under ADMM and only stalls under ALADIN. And a Newton-type step
    must be globalized, which adds a line search and its own tuning.

## Lirion vs. the original

!!! warning "Lirion vs. the original"
    Lirion follows Algorithm 2 of [Houska2016](@cite) in structure, with these
    departures.

    **The matrix ``B_r``.** The paper prescribes ``B_r \approx \nabla^2 [f_r +
    \kappa_r^{\top} h_r]``, with the proximal term living in the local
    subproblem only. By default Lirion adds the proximal curvature
    ``\rho\Sigma_r`` to ``B_r`` on AC, SOCWR, and SOCWRNL, plus a ``10^{-6}``
    Tikhonov floor and whatever ``\delta I`` shift is needed for positive
    semidefiniteness; on DC ``B_r`` is the cost diagonal over the same ``10^{-6}``
    floor, and because the DC constraints are linear this is the exact Lagrangian
    Hessian. Adding
    ``\rho\Sigma_r`` damps the correction and slows the local rate at large
    ``\rho``, and it is the choice the shipped campaigns were run with.

    **The conic SOCWR.** The paper assumes twice-differentiable ``f_r`` and
    ``h_r``. The cones of the conic SOCWR are not, so its ``B_r`` carries no cone
    curvature and its ``C_r`` only the active cone-face normals; the coordinator step then
    makes little progress and the run settles near the conic solver's accuracy.
    ALADIN rejects `:dummy_bus` on `socwr()` for this reason, and SOCWRNL, which
    writes the cones as smooth inequalities, is the supported route to the Jabr
    relaxation under ALADIN.

    **The gradient.** The paper's step 3 uses the modified gradient ``g_r =
    \nabla f_r + (C_r^{*} - C_r)^{\top}\kappa_r``, correcting for inequalities
    left out of the active set. Lirion uses the plain cost gradient by default;
    the correction is an opt-in on the smooth NLP models.

    **The stop.** The paper uses one tolerance and the ``\ell_1`` norm on both
    residuals. Lirion keeps the ``\ell_1`` stationarity residual with a
    dimension-scaled bar ``n_{\mathrm{prox}}\,\epsilon^{\mathrm{abs}}``, and
    measures the consensus residual either by Boyd's ``\sqrt{n}``-scaled 2-norm
    test (`:boyd`, the default) or by a flat bar in the paper's style (`:hfd`),
    measured in the max-norm rather than the paper's 1-norm. Both options test
    consensus and stationarity.

    **Globalization.** The paper's Algorithm 3 is an exact ``\ell_1`` penalty
    with three branches. Lirion's default is an Armijo backtracking on
    ``\alpha_2`` alone, over the coordinator's augmented Lagrangian; Algorithm 3
    is available as `:alg3` and, with the paper's Wolfe branch plus a gradient
    fallback when the QP dual direction is not an ascent direction and a penalty
    weight fixed for the run, `:alg3_original`.

    **Active set.** Active bounds are detected by primal distance (``10^{-6}``)
    and active inequalities by the size of their local multiplier, then
    dependent rows are removed. The paper leaves the choice of ``C_r \approx
    C_r^{*}`` open.

## See also

- Neighbours: [ADMM](admm.md) (first-order, boundary-values-only exchange,
  carries the conic SOCWR) · [Decomposition: sectioning and coupling](seccionamento_acoplamento.md)
  (what the rows of ``A_r`` contain for each cut and model).
- Foundations: [Foundations](@ref).
- Usage: [ALADIN settings](../../usage/distributed/aladin.md) for `rho`, `mu`,
  `proximal_mode`, `globalization`, and `stopping_criterion`.
- Code: constructor `aladin()`; source `src/distributed_algorithms/aladin/`
  (`core/coordinator_qp.jl`, `core/residuals.jl`, `core/tolerances.jl`, and one
  folder per model).
- Original: [Houska2016](@cite).
