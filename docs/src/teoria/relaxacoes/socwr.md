# SOCWR (Jabr conic relaxation)

The SOCWR model is a **convex relaxation** of the exact AC power flow, written
in **W-space**. Every model up to this point ([AC Polar](@ref),
[AC Rectangular (ACR)](@ref), [IVR](@ref)) was exact, and paid for it with
non-convexity: no guarantee of a global optimum, and sensitivity to the
starting point. That difficulty is intrinsic and does not vanish by changing
coordinates. SOCWR trades exactness for convexity. It keeps the ACR
branch-flow expressions verbatim but promotes the four W blocks ``W_{ii}``,
``W_{jj}``, ``W_r``, ``W_i`` (the code identifiers `wii`, `wjj`, `wr`, `wi`)
from *expressions of the voltages* to *free decision variables*, tied together
by a single rotated second-order cone. Its optimum is therefore a lower bound
on the true AC optimum. The relaxation is due to Jabr [Jabr2006](@cite).

!!! note "Relaxation or approximation?"
    SOCWR is a **relaxation**, not an approximation: it enlarges the feasible
    region (the binding equality ``=`` becomes an inequality ``\le``), so with
    more freedom a cost-minimizing solution can only fall or stay equal. Its
    optimum is a **guaranteed lower bound** on the true AC optimum: a floor on
    the cost, never an estimate that might land above the true value. (Contrast the linear *approximations*, whose optima
    are neither a floor nor a ceiling.)

## Derivation

### The idea: relax

Relaxing means solving a different problem, easier and convex, that is close
enough to the original to be useful. Concretely, we **loosen the constraints**,
enlarging the feasible region. Because the objective minimizes cost, giving the
problem more freedom can only lower the optimum or leave it unchanged, so the
relaxed optimum is a **lower bound** of the true one.

### Promoting the blocks to variables

In the [AC Rectangular (ACR)](@ref) model, the branch-flow formulas are
**linear in the W blocks**. SOCWR reuses exactly those same four blocks. What
changes is not the flow but the *status* of the blocks: they become variables
rather than expressions. The branch parameters ``g, b`` (series
conductance and susceptance), ``g^{fr}, b^{fr}, g^{to}, b^{to}`` (end shunts),
the complex-tap parts ``t_r, t_i`` and ``tm = |T|^2 = \tau^2`` are exactly as
in the ACR / branch model.

### The flows (identical to the ACR)

The flow formulas are **identical to the ACR**. The only difference is that in
the ACR ``W_{ii}``, ``W_{jj}``, ``W_r`` and ``W_i`` are expressions of
``v^{\mathrm{re}}`` and ``v^{\mathrm{im}}``, whereas in SOCWR they are free
variables (tied by the cone).

**Active power, *from* side:**

```math
p_{ij} = \frac{g + g^{fr}}{tm}\,W_{ii}
+ \frac{-g\,t_r + b\,t_i}{tm}\,W_r
+ \frac{-b\,t_r - g\,t_i}{tm}\,W_i .
```

**Reactive power, *from* side:**

```math
q_{ij} = -\,\frac{b + b^{fr}}{tm}\,W_{ii}
- \frac{-b\,t_r - g\,t_i}{tm}\,W_r
+ \frac{-g\,t_r + b\,t_i}{tm}\,W_i .
```

**Active power, *to* side:**

```math
p_{ji} = (g + g^{to})\,W_{jj}
+ \frac{-g\,t_r - b\,t_i}{tm}\,W_r
+ \frac{b\,t_r - g\,t_i}{tm}\,W_i .
```
**Reactive power, *to* side:**

```math
q_{ji} = -(b + b^{to})\,W_{jj}
- \frac{-b\,t_r + g\,t_i}{tm}\,W_r
- \frac{-g\,t_r - b\,t_i}{tm}\,W_i .
```

### The blocks come from voltages, and can lie

In the ACR, the four blocks are **quadratic** (hence non-convex) expressions of
the rectangular voltages:

```math
W_{ii} = (v^{\mathrm{re}}_i)^2 + (v^{\mathrm{im}}_i)^2, \qquad
W_{jj} = (v^{\mathrm{re}}_j)^2 + (v^{\mathrm{im}}_j)^2,
```

```math
W_r = v^{\mathrm{re}}_i\,v^{\mathrm{re}}_j + v^{\mathrm{im}}_i\,v^{\mathrm{im}}_j, \qquad
W_i = v^{\mathrm{im}}_i\,v^{\mathrm{re}}_j - v^{\mathrm{re}}_i\,v^{\mathrm{im}}_j .
```

The flows are *constant* ``\times`` *block*, so SOCWR's proposal is to forget the
voltages ``v^{\mathrm{re}}, v^{\mathrm{im}}`` and keep only the blocks as free
variables.

!!! note "The loose variables lie"
    If ``W_{ii}``, ``W_{jj}``, ``W_r`` and ``W_i`` are left free, the optimizer
    picks combinations that correspond to *no voltage at all*. We need a
    constraint forcing the blocks to be consistent, so they come from genuine
    voltages.

### The Jabr tie

The blocks come from real voltages, with magnitudes ``v_i, v_j`` and angle
difference ``\theta_{ij} = \theta_i - \theta_j``:

```math
W_{ii} = v_i^2, \quad W_{jj} = v_j^2, \quad
W_r = v_i v_j \cos\theta_{ij}, \quad W_i = v_i v_j \sin\theta_{ij} .
```

These are the rectangular blocks of the previous section, rewritten with
``v^{\mathrm{re}} = v\cos\theta`` and ``v^{\mathrm{im}} = v\sin\theta``.

Forming ``W_r^2 + W_i^2``:

```math
W_r^2 + W_i^2 = \left(v_i v_j \cos\theta_{ij}\right)^2 + \left(v_i v_j \sin\theta_{ij}\right)^2
= v_i^2\,v_j^2\,\underbrace{\left(\cos^2\theta_{ij} + \sin^2\theta_{ij}\right)}_{=\,1} ,
```

so the four blocks are bound by the exact **Jabr tie**:

```math
W_r^2 + W_i^2 = W_{ii}\,W_{jj} .
```

This tie is the rank-one consistency of the voltage-product matrix ``W = V V^{*}``.
With ``W_{ii}, W_{jj}`` on the diagonal and ``W_r + j\,W_i`` off the diagonal,
``W_r^2 + W_i^2 = W_{ii} W_{jj}`` is its zero-determinant condition. SOCWR keeps
this per bus pair. The SDP relaxation instead keeps the whole matrix positive
semidefinite.

### The relaxation: `=` becomes `≤`

The tie is an equality between quadratic quantities, so it is non-convex. SOCWR
replaces the equality by an inequality,

```math
W_r^2 + W_i^2 \le W_{ii}\,W_{jj} ,
```

which turns the non-convex set into a convex one, a **rotated second-order cone**.
A rotated cone is the set of ``(u, v, x)`` with ``u, v \ge 0`` and a vector ``x``
obeying ``\|x\|^2 \le 2uv``. Taking ``u = W_{ii}/\sqrt{2}``, ``v = W_{jj}/\sqrt{2}`` and
``x = (W_r, W_i)``, the condition ``\|x\|^2 \le 2uv`` is exactly the tie above, since
the two ``\sqrt{2}`` cancel the factor 2.

!!! note "Shell versus solid"
    Geometrically, the equality describes a *shell* (a curved surface); the
    inequality ``\le`` describes its *solid interior*. The shell is non-convex;
    the solid is convex. Trading the surface for the solid is exactly what buys
    convexity.

### What the relaxation loses

Writing both sides back in physical terms shows what the inequality gives up:

```math
\underbrace{W_r^2 + W_i^2}_{(v_i v_j)^2\,(\cos^2\theta_{ij} + \sin^2\theta_{ij})}
\;\le\;
\underbrace{W_{ii}\,W_{jj}}_{(v_i v_j)^2} .
```

Points with ``W_r^2 + W_i^2 < W_{ii} W_{jj}`` satisfy the cone but correspond to
``\cos^2\theta_{ij} + \sin^2\theta_{ij} < 1``, which is physically impossible.
Geometrically, such a point has been pulled *inside the unit circle*. The model
starts to admit impossible solutions, whose blocks, translated back into voltage
and angle, fall outside the permitted ranges. The LNC cuts below tighten the
feasible region, removing some of these spurious points.

## Constraints

### Thermal limit (conic form)

The apparent-power limit is imposed directly as a second-order cone at each end
of the branch:

```math
\left(\mathrm{rate}_a,\; p_{ij},\; q_{ij}\right) \in \mathcal{Q}_3,
\qquad
\left(\mathrm{rate}_a,\; p_{ji},\; q_{ji}\right) \in \mathcal{Q}_3,
```

which is simply

```math
\sqrt{p_{ij}^2 + q_{ij}^2} \le \mathrm{rate}_a, \qquad
\sqrt{p_{ji}^2 + q_{ji}^2} \le \mathrm{rate}_a .
```

### Angle limits

In W-space the aperture limits are expressed through the blocks, using
``\tan\theta_{ij} = W_i / W_r``:

```math
\tan\theta^{\min}\,W_r \le W_i \le \tan\theta^{\max}\,W_r .
```

## LNC cuts (Lifted Nonlinear Cuts)

The conic tie ``W_r^2 + W_i^2 \le W_{ii} W_{jj}`` is **loose**. It uses only the
relation among the four blocks and ignores the voltage and angle limits. A
point can satisfy the cone and still have a magnitude or an angle outside its
physical range. The **LNC cuts** inject that missing information.

!!! note "Derived here"
    The LNC cuts are derived in full on this page. The [QC (Quadratic Convex)](@ref)
    relaxation reuses the very same two cuts and cross-references this
    derivation rather than repeating it.

### McCormick envelopes

McCormick is a general technique for products of variables ``a \cdot b``, which
are non-convex. Instead of imposing the exact product ``z = ab`` (a curved
surface, hard), it replaces that equality by four linear inequalities that box
the product inside a convex region. One works with the box (easy, linear)
instead of the curve (hard, non-convex), gaining convexity and paying with a
gap between the box and the true product.

Each of the four bounds comes from a **product of factors of known sign**,
expanded term by term: two lower bounds (factors on the same side) and two
upper bounds (factors on opposite sides).

**Lower bound 1.** Both factors at their minimums. Since ``a \ge a^{\min}`` and
``b \ge b^{\min}``, both factors are ``\ge 0``:

```math
(a - a^{\min})(b - b^{\min}) \ge 0 \;\Longrightarrow\;
ab \ge a^{\min} b + b^{\min} a - a^{\min} b^{\min} .
```

**Lower bound 2.** Both factors at their maximums. Since ``a \le a^{\max}`` and
``b \le b^{\max}``, the factors ``(a^{\max}-a)`` and ``(b^{\max}-b)`` are
``\ge 0``:

```math
(a^{\max} - a)(b^{\max} - b) \ge 0 \;\Longrightarrow\;
ab \ge a^{\max} b + b^{\max} a - a^{\max} b^{\max} .
```

**Upper bound 1.** Opposite sides (``a`` at its minimum, ``b`` at its maximum).
One factor ``\ge 0``, the other ``\le 0``, so the product is ``\le 0``:

```math
(a - a^{\min})(b - b^{\max}) \le 0 \;\Longrightarrow\;
ab \le a^{\min} b + b^{\max} a - a^{\min} b^{\max} .
```
**Upper bound 2.** Opposite sides (``a`` at its maximum, ``b`` at its minimum):

```math
(a - a^{\max})(b - b^{\min}) \le 0 \;\Longrightarrow\;
ab \le a^{\max} b + b^{\min} a - a^{\max} b^{\min} .
```

In SOCWR we use only the **two lower bounds**. The two upper bounds enter in the
[QC (Quadratic Convex)](@ref) relaxation.

### Building the cut

The cut is assembled in six steps.

#### Step 1: build the cut from the blocks

We would like to use ``\theta_{ij} \in [\theta^{\min}, \theta^{\max}]`` directly,
but the cut has to be written in the relaxation's blocks ``W_r, W_i``. A cosine
and a sine form a *direction* in the plane, so projecting the blocks onto a
reference direction ``\alpha`` (a dot product) recovers an angle:

```math
\cos\alpha\,\cos\theta_{ij} + \sin\alpha\,\sin\theta_{ij} = \cos(\theta_{ij} - \alpha),
```

```math
\cos\alpha\,W_r + \sin\alpha\,W_i = v_i v_j \cos(\theta_{ij} - \alpha).
```

#### Step 2: which reference `α`?

We are free to choose ``\alpha``, and we want the **strongest** cut. Because the
cosine will become the floor of the final inequality, we want the largest
possible cosine. As ``\theta_{ij}`` ranges over ``[\theta^{\min},
\theta^{\max}]``, the argument ``\theta_{ij} - \alpha`` ranges over
``[\theta^{\min} - \alpha,\ \theta^{\max} - \alpha]``, and the cosine is smaller
the farther the argument sits from zero. The worst case (smallest cosine) is the
end farthest from zero. To make the bound as good as possible, we minimize that
maximum distance.

!!! note "Why the center"
    The choice that makes the band most symmetric about zero is its **center**.
    Shifting by the center turns the band into ``[-d, +d]``, so the maximum
    distance to zero is as small as possible and the worst-case ``\cos d`` is as
    large as possible:

    ```math
    \Phi = \frac{\theta^{\max} + \theta^{\min}}{2} \;(\text{center}), \qquad
    d = \frac{\theta^{\max} - \theta^{\min}}{2} \;(\text{half-width}).
    ```

    Saying ``\theta_{ij} \in [\theta^{\min}, \theta^{\max}]`` is the same as
    ``|\theta_{ij} - \Phi| \le d``.

#### Step 3: the angle bound

With ``|\theta_{ij} - \Phi| \le d`` the argument lies in ``[-d, +d]``, where the
cosine is minimal at the ends:

```math
\cos(\theta_{ij} - \Phi) \ge \cos d .
```

Combining with Step 1:

```math
\cos\Phi\,W_r + \sin\Phi\,W_i = v_i v_j \cos(\theta_{ij} - \Phi)
\;\ge\; v_i v_j \cos d .
```

!!! note "Why a floor (≥)"
    SOCWR's slack is *downward*: the blocks can shrink
    (``\cos^2 + \sin^2 < 1``). The useful cut is therefore a **floor**. It forces
    the projection of the blocks to be at least ``\cos d \cdot v_i v_j``.

#### Step 4: the magnitude ↔ magnitude² bridge

The right-hand side carries ``v_i v_j`` (nonlinear). To make the cut linear in
the blocks, we must replace ``v_i v_j`` by something in ``W_{ii}, W_{jj}``.
First we link each magnitude ``v_i`` to its square ``W_{ii}``, using
opposite-sign factors of the *same* variable (since ``v_i \in [v_i^{\min},
v_i^{\max}]``):

```math
(v_i - v_i^{\min})(v_i - v_i^{\max}) \le 0 .
```

Expanding, with ``W_{ii} = v_i^2`` and ``s_f = v_i^{\min} + v_i^{\max}``:

```math
W_{ii} - s_f\,v_i + v_i^{\min} v_i^{\max} \le 0
\;\Longrightarrow\;
v_i \ge \frac{W_{ii} + v_i^{\min} v_i^{\max}}{s_f} .
```

Likewise for ``v_j``, with ``W_{jj} = v_j^2`` and ``s_t = v_j^{\min} + v_j^{\max}``:

```math
(v_j - v_j^{\min})(v_j - v_j^{\max}) \le 0
\;\Longrightarrow\;
v_j \ge \frac{W_{jj} + v_j^{\min} v_j^{\max}}{s_t} .
```

The sums ``s_f`` and ``s_t`` are the *from* and *to* bound-sums of the branch.

#### Step 5: McCormick of `v_i v_j` (the two origins)

We need a *lower* bound of ``v_i v_j`` (to chain with the ``\ge``). The two lower
McCormick bounds come from same-side factors, which is where the two cuts
originate:

**Origin of Cut 1.** Upper factors ``(v_i^{\max} - v_i)(v_j^{\max} - v_j) \ge 0``:

```math
v_i v_j \ge v_i^{\max} v_j + v_j^{\max} v_i - v_i^{\max} v_j^{\max} .
```

**Origin of Cut 2.** Lower factors ``(v_i - v_i^{\min})(v_j - v_j^{\min}) \ge 0``:

```math
v_i v_j \ge v_i^{\min} v_j + v_j^{\min} v_i - v_i^{\min} v_j^{\min} .
```

#### Step 6: assembling the cut

For each origin, multiply by ``s_f s_t`` (clearing the bridge denominators),
substitute the Step 4 bridges, and collect the constants. From the origin of
Cut 1 this yields a lower bound of ``v_i v_j`` that is linear in ``W_{ii},
W_{jj}``:

```math
s_f s_t\,v_i v_j \ge v_j^{\max} s_t\,W_{ii} + v_i^{\max} s_f\,W_{jj}
+ v_i^{\max} v_j^{\max}\left(v_i^{\min} v_j^{\min} - v_i^{\max} v_j^{\max}\right).
```

Chaining with Step 3 gives **Cut 1**:

```math
s_f s_t\left(\cos\Phi\,W_r + \sin\Phi\,W_i\right)
- v_j^{\max}\cos d\; s_t\,W_{ii}
- v_i^{\max}\cos d\; s_f\,W_{jj}
\;\ge\;
v_i^{\max} v_j^{\max}\cos d\left(v_i^{\min} v_j^{\min} - v_i^{\max} v_j^{\max}\right).
```

### Cut 2

The same calculation, starting from the origin of Cut 2 (multiply by
``s_f s_t``, substitute the bridges, collect the constants), gives a second
lower bound of ``v_i v_j`` linear in ``W_{ii}, W_{jj}``:

```math
s_f s_t\,v_i v_j \ge v_j^{\min} s_t\,W_{ii} + v_i^{\min} s_f\,W_{jj}
- v_i^{\min} v_j^{\min}\left(v_i^{\min} v_j^{\min} - v_i^{\max} v_j^{\max}\right),
```

and, chained with Step 3, **Cut 2**:

```math
s_f s_t\left(\cos\Phi\,W_r + \sin\Phi\,W_i\right)
- v_j^{\min}\cos d\; s_t\,W_{ii}
- v_i^{\min}\cos d\; s_f\,W_{jj}
\;\ge\;
-\,v_i^{\min} v_j^{\min}\cos d\left(v_i^{\min} v_j^{\min} - v_i^{\max} v_j^{\max}\right).
```

## SOCWR Formulation

The SOCWR constructor (`socwr()`, model type `SOCWR`) assembles the standard
skeleton (buses, branches, generators, storage, active/reactive balance,
objective) and dispatches the branch block to the conic W-space relaxation.
Without storage the problem is a pure **second-order cone program**, which Lirion
solves with `Clarabel` by default. The storage block adds a charge/discharge
binary (Eq. 9), so a storage-bearing instance is a **mixed-integer SOCP** and
needs a MISOCP-capable solver. The benchmark cases carry no storage.

### Sets

- ``\mathcal{N}``: set of buses.
- ``\mathcal{N}^{nv} \subseteq \mathcal{N}``: non-virtual buses, on which the power balance is enforced (`non_virtual_bus_ids`).
- ``\mathcal{L}``: set of active branches (`br.active`, in-service lines and
  transformers).
- ``\mathcal{P}``: set of distinct ordered bus pairs ``(\mathrm{fr},\mathrm{to})``
  carrying at least one active branch. Parallel branches share one pair.
- ``\mathcal{G}``: set of existing generators.
- ``\mathcal{G}^c \subseteq \mathcal{G}``: curtailable generators (`plant_type == "R"`).
- ``\mathcal{S}``: set of storage units.
- ``\mathcal{T} = \{1,\dots,T\}``: time periods.

### Indices and incidence

- ``b \in \mathcal{N}``: bus index; ``\ell \in \mathcal{L}``: branch index;
  ``p \in \mathcal{P}``: bus-pair index; ``g``: generator; ``s``: storage;
  ``t``: time.
- ``\mathrm{fr}(\ell),\ \mathrm{to}(\ell)``: *from* and *to* buses of branch
  ``\ell`` (`br.from`, `br.to`). For a branch ``\ell`` we abbreviate
  ``i = \mathrm{fr}(\ell)``, ``j = \mathrm{to}(\ell)``.
- ``p(\ell) \in \mathcal{P}``: the bus pair of branch ``\ell``;
  ``i(p),\ j(p)``: the two buses of pair ``p``.
- ``B(s)``: connection bus of storage ``s``.

### Parameters

**Branches.**

- ``g_\ell,\ b_\ell``: series conductance / susceptance,
  ``g = r/(r^2+x^2)``, ``b = -x/(r^2+x^2)``.
- ``g^{fr}_\ell,\ b^{fr}_\ell``: *from*-side shunt (code sets ``g^{fr}=0``,
  ``b^{fr} = b^{sh}_\ell/2``).
- ``g^{to}_\ell,\ b^{to}_\ell``: *to*-side shunt (``g^{to}=0``,
  ``b^{to} = b^{sh}_\ell/2``).
- ``t_r^\ell,\ t_i^\ell``: real / imaginary parts of the complex tap,
  ``t_r = \tau\cos\theta_{shift}``, ``t_i = \tau\sin\theta_{shift}``
  (`tap_real`, `tap_imag`).
- ``tm_\ell = |T_\ell|^2 = \tau_\ell^2``, with ``1/tm_\ell`` stored as
  ``\texttt{inv\_tap}_\ell^2 = (1/\tau_\ell)^2``.
- ``\mathrm{rate}_a^\ell``: apparent-power (thermal) limit (`rate_a`).
- ``\theta^{\min}_\ell,\ \theta^{\max}_\ell``: angle-difference limits in radians
  (`ang_min`, `ang_max`); either may be `NaN` (absent).

**Buses.**

- ``V_b^{\min},\ V_b^{\max}``: voltage-magnitude limits.
- ``g^{sh}_b,\ b^{sh}_b``: bus shunt conductance / susceptance
  (`shunt_conductance`, `shunt_susceptance`).
- ``P^d_{b,t},\ Q^d_{b,t}``: active / reactive demand.

**Generators.** As in the exact models: static limits ``P_g^{\min},P_g^{\max}``,
``Q_g^{\min},Q_g^{\max}``; per-period envelopes
``\overline{P}^{ts}_{g,t},\underline{P}^{ts}_{g,t}``,
``\overline{Q}^{ts}_{g,t},\underline{Q}^{ts}_{g,t}``; cost coefficients
``c_{g,t,k}``; curtailment cost ``c^{curt}_g``.

**Storage.** ``E_s^{\min},E_s^{\max}``, initial energy ``E^0_s``; charge /
discharge ramps ``C^{\max}_s,D^{\max}_s``; efficiencies ``\eta^c_s,\eta^d_s``;
converter ``r_s,x_s``; fixed losses ``P^{loss}_s,Q^{loss}_s``; inverter limit
``\overline{S}_s``; reactive limits ``Q_s^{\min},Q_s^{\max}``; the global voltage
floor ``V^{\min}_{\mathrm{glob}}`` (`model.ext[:v_min]`, default ``0.9``) used in
the current bound.

### Decision variables

The relaxation lives in **W-space**: the bus squares ``W_{bb}`` and the per-pair
products ``W^{r}_p, W^{i}_p`` are decision variables, tied by the Jabr cone and
carrying the box bounds derived below. For a branch ``\ell`` with
``i=\mathrm{fr}(\ell)``, ``j=\mathrm{to}(\ell)``, ``p=p(\ell)``, the four blocks
of the derivation are the views ``W_{ii}=W_{ii,t}``, ``W_{jj}=W_{jj,t}``,
``W_r=W^{r}_{p,t}``, ``W_i=W^{i}_{p,t}``.

```math
\theta_{b,t} \in \mathbb{R}\;\text{: auxiliary voltage angle at bus } b. \tag{1}
```

```math
W_{bb,t} \in \big[(V_b^{\min})^2,\ (V_b^{\max})^2\big]
\;\text{: squared voltage magnitude (}\texttt{v\_magnitude\_sq}\text{).} \tag{2}
```

```math
W^{r}_{p,t} \in [R^{\min}_p,\ R^{\max}_p]
\;\text{: real voltage product of pair } p\ (\texttt{R}). \tag{3}
```

```math
W^{i}_{p,t} \in [T^{\min}_p,\ T^{\max}_p]
\;\text{: imaginary voltage product of pair } p\ (\texttt{T\_var}). \tag{4}
```

```math
p^g_{g,t} \in [P_g^{\min}, P_g^{\max}], \qquad
Q^g_{g,t} \in [Q_g^{\min}, Q_g^{\max}]. \tag{5}
```

```math
p^c_{g,t} \in [0, \overline{P}^{ts}_{g,t}], \qquad g\in\mathcal{G}^c. \tag{6}
```

```math
E_{s,t} \in [E_s^{\min}, E_s^{\max}], \quad
c_{s,t} \in [0, C^{\max}_s], \quad
d_{s,t} \in [0, D^{\max}_s]. \tag{7}
```

```math
P^s_{s,t} \in [-\overline{S}_s, \overline{S}_s], \quad
Q^s_{s,t} \in [Q_s^{\min}, Q_s^{\max}], \quad
Q^{sc}_{s,t} \in [Q_s^{\min}, Q_s^{\max}]. \tag{8}
```

```math
z_{s,t} \in \{0,1\}\;\text{(charge/discharge indicator)}, \qquad
\mathrm{ccm}_{s,t} \in \Big[0,\ \big(\overline{S}_s/V^{\min}_{\mathrm{glob}}\big)^2\Big]
\;\text{(current magnitude}^2,\ \texttt{ccms}). \tag{9}
```

### Product bounds

The pair products are not fully free. Bounds (3)–(4) box each one, computed per
pair from the bus voltage-magnitude limits and the pair's angle window (the
PowerModels rule `ref_calc_voltage_product_bounds`). Across the parallel branches
of a pair ``p`` the window is the tightest of them,
``\theta^{\min}_p = \max_\ell \theta^{\min}_\ell`` and
``\theta^{\max}_p = \min_\ell \theta^{\max}_\ell``. With ``i=i(p)``, ``j=j(p)``
and both limits finite, the bounds split on the sign of the window. The reasoning
is monotonic. ``R = v_i v_j\cos\theta`` is largest when both magnitudes sit at
their maximum and ``|\theta|`` is smallest, and smallest at the opposite corner.
The imaginary product ``T = v_i v_j\sin\theta`` follows ``\sin``, which grows with
``\theta``. Because ``\cos`` is even and ``\sin`` is odd, the corner that gives each extreme
depends on the sign of the window, which is why there are three cases.

When ``\theta^{\min}_p \ge 0`` (the window sits at or above zero):

```math
\begin{aligned}
R^{\max}_p &= V_i^{\max} V_j^{\max}\cos\theta^{\min}_p, \\
R^{\min}_p &= V_i^{\min} V_j^{\min}\cos\theta^{\max}_p, \\
T^{\max}_p &= V_i^{\max} V_j^{\max}\sin\theta^{\max}_p, \\
T^{\min}_p &= V_i^{\min} V_j^{\min}\sin\theta^{\min}_p.
\end{aligned}
```

When ``\theta^{\max}_p \le 0`` (at or below zero):

```math
\begin{aligned}
R^{\max}_p &= V_i^{\max} V_j^{\max}\cos\theta^{\max}_p, \\
R^{\min}_p &= V_i^{\min} V_j^{\min}\cos\theta^{\min}_p, \\
T^{\max}_p &= V_i^{\min} V_j^{\min}\sin\theta^{\max}_p, \\
T^{\min}_p &= V_i^{\max} V_j^{\max}\sin\theta^{\min}_p.
\end{aligned}
```

Otherwise the window straddles zero (``\theta^{\min}_p < 0 < \theta^{\max}_p``):
``R^{\max}_p`` reaches the full magnitude product, and the others follow the
extreme magnitudes,

```math
\begin{aligned}
R^{\max}_p &= V_i^{\max} V_j^{\max}, \\
R^{\min}_p &= V_i^{\min} V_j^{\min}\,\min\!\big(\cos\theta^{\min}_p,\ \cos\theta^{\max}_p\big), \\
T^{\max}_p &= V_i^{\max} V_j^{\max}\sin\theta^{\max}_p, \\
T^{\min}_p &= V_i^{\max} V_j^{\max}\sin\theta^{\min}_p.
\end{aligned}
```

When a pair's angle limits are not both finite, no product box is imposed and both
``W^{r}_p`` and ``W^{i}_p`` stay free.

### Objective

Identical to the exact models. The fuel polynomial is evaluated on the power in
physical units (``S^{base}\,p^g``), and curtailment is priced linearly.

```math
F = \sum_{t\in\mathcal{T}} \sum_{g\in\mathcal{G}} \sum_{k=1}^{m_{g,t}}
c_{g,t,k}\left(S^{base}\,p^g_{g,t}\right)^{m_{g,t} - k}. \tag{10}
```

```math
F^{curt} = \sum_{g\in\mathcal{G}^c} \sum_{t\in\mathcal{T}}
c^{curt}_g\,S^{base}\,p^c_{g,t}, \qquad \min\; F + F^{curt}. \tag{11}
```

### Branch flows

The flows are **linear in the blocks** and identical to the ACR (the only change
is the status of the blocks). For ``\ell\in\mathcal{L}``, ``t\in\mathcal{T}``,
with ``i=\mathrm{fr}(\ell)``, ``j=\mathrm{to}(\ell)`` and the per-pair products
``W_r=W^{r}_{p(\ell),t}``, ``W_i=W^{i}_{p(\ell),t}``
(`p_flow_branches`, `q_flow_branches`, `p_reverse_flow_branches`,
`q_reverse_flow_branches`):

```math
p_{ij} = \frac{g_\ell + g^{fr}_\ell}{tm_\ell}\,W_{ii,t}
+ \frac{-g_\ell t_r^\ell + b_\ell t_i^\ell}{tm_\ell}\,W_r
+ \frac{-b_\ell t_r^\ell - g_\ell t_i^\ell}{tm_\ell}\,W_i. \tag{12}
```

```math
q_{ij} = -\frac{b_\ell + b^{fr}_\ell}{tm_\ell}\,W_{ii,t}
- \frac{-b_\ell t_r^\ell - g_\ell t_i^\ell}{tm_\ell}\,W_r
+ \frac{-g_\ell t_r^\ell + b_\ell t_i^\ell}{tm_\ell}\,W_i. \tag{13}
```

```math
p_{ji} = (g_\ell + g^{to}_\ell)\,W_{jj,t}
+ \frac{-g_\ell t_r^\ell - b_\ell t_i^\ell}{tm_\ell}\,W_r
+ \frac{b_\ell t_r^\ell - g_\ell t_i^\ell}{tm_\ell}\,W_i. \tag{14}
```

```math
q_{ji} = -(b_\ell + b^{to}_\ell)\,W_{jj,t}
- \frac{-b_\ell t_r^\ell + g_\ell t_i^\ell}{tm_\ell}\,W_r
- \frac{-g_\ell t_r^\ell - b_\ell t_i^\ell}{tm_\ell}\,W_i. \tag{15}
```

The *to*-side self terms ``(g+g^{to})W_{jj}`` and ``-(b+b^{to})W_{jj}`` carry
**no** ``1/tm`` factor, exactly as in the code.

### Power balance

On each non-virtual bus ``b\in\mathcal{N}^{nv}``, ``t\in\mathcal{T}``, the bus
shunt enters through the squared magnitude ``W_{bb,t}`` (storage injects with a
minus sign, ``P^s>0`` absorbs from the grid):

```math
\sum_{g\in\mathcal{G}(b)} p^g_{g,t}
- \sum_{s\in\mathcal{S}(b)} P^s_{s,t}
- P^d_{b,t}
- \sum_{\ell:\,\mathrm{fr}(\ell)=b} p_{ij}
- \sum_{\ell:\,\mathrm{to}(\ell)=b} p_{ji}
- g^{sh}_b\,W_{bb,t} = 0,
\qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{16}
```

```math
\sum_{g\in\mathcal{G}(b)} Q^g_{g,t}
- \sum_{s\in\mathcal{S}(b)} Q^s_{s,t}
- Q^d_{b,t}
- \sum_{\ell:\,\mathrm{fr}(\ell)=b} q_{ij}
- \sum_{\ell:\,\mathrm{to}(\ell)=b} q_{ji}
+ b^{sh}_b\,W_{bb,t} = 0,
\qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{17}
```

### Jabr cone

One rotated second-order cone **per bus pair**, not per branch. The code groups
parallel branches and emits the cone once, on the shared pair products, through
the representative branch of each pair.

```math
\Big(\tfrac{W_{i(p)\,i(p),t}}{\sqrt{2}},\ \tfrac{W_{j(p)\,j(p),t}}{\sqrt{2}},\
W^{r}_{p,t},\ W^{i}_{p,t}\Big) \in \mathcal{Q}^{\mathrm{rot}}_4,
\qquad p\in\mathcal{P},\ t\in\mathcal{T}, \tag{18}
```

which is exactly the relaxed Jabr tie, since
``2\cdot\tfrac{W_{i(p)i(p),t}}{\sqrt{2}}\cdot\tfrac{W_{j(p)j(p),t}}{\sqrt{2}} = W_{i(p)i(p),t}\,W_{j(p)j(p),t}``:

```math
(W^{r}_{p,t})^2 + (W^{i}_{p,t})^2 \le W_{i(p)\,i(p),t}\,W_{j(p)\,j(p),t}.
```

### LNC cuts (Lifted Nonlinear Cuts)

Added **per active branch** ``\ell`` only when both ``\theta^{\min}_\ell`` and
``\theta^{\max}_\ell`` are finite, both lie in ``[-\tfrac{\pi}{2},\tfrac{\pi}{2}]``,
and ``\theta^{\min}_\ell < \theta^{\max}_\ell``. With
``\Phi_\ell = \tfrac{\theta^{\max}_\ell+\theta^{\min}_\ell}{2}``,
``d_\ell = \tfrac{\theta^{\max}_\ell-\theta^{\min}_\ell}{2}``,
``s_f = V_i^{\min}+V_i^{\max}``, ``s_t = V_j^{\min}+V_j^{\max}``
(``i=\mathrm{fr}(\ell)``, ``j=\mathrm{to}(\ell)``, ``W_r=W^{r}_{p(\ell),t}``,
``W_i=W^{i}_{p(\ell),t}``):

```math
s_f s_t\big(\cos\Phi_\ell\,W_r + \sin\Phi_\ell\,W_i\big)
- V_j^{\max}\cos d_\ell\, s_t\,W_{ii,t}
- V_i^{\max}\cos d_\ell\, s_f\,W_{jj,t}
\ge V_i^{\max} V_j^{\max}\cos d_\ell\big(V_i^{\min} V_j^{\min} - V_i^{\max} V_j^{\max}\big). \tag{19}
```

```math
s_f s_t\big(\cos\Phi_\ell\,W_r + \sin\Phi_\ell\,W_i\big)
- V_j^{\min}\cos d_\ell\, s_t\,W_{ii,t}
- V_i^{\min}\cos d_\ell\, s_f\,W_{jj,t}
\ge -\,V_i^{\min} V_j^{\min}\cos d_\ell\big(V_i^{\min} V_j^{\min} - V_i^{\max} V_j^{\max}\big). \tag{20}
```

### Branch thermal limits (SOC)

One second-order cone at each end of every active branch:

```math
\big(\mathrm{rate}_a^\ell,\ p_{ij},\ q_{ij}\big) \in \mathcal{Q}_3, \qquad
\big(\mathrm{rate}_a^\ell,\ p_{ji},\ q_{ji}\big) \in \mathcal{Q}_3, \tag{21}
```

i.e. ``\sqrt{p_{ij}^2 + q_{ij}^2} \le \mathrm{rate}_a^\ell`` and
``\sqrt{p_{ji}^2 + q_{ji}^2} \le \mathrm{rate}_a^\ell``.

### Angle limits

The code carries an auxiliary bus angle ``\theta`` (kept so the model can state
the raw aperture, pin the reference bus, and interface with the distributed
coupling). Its raw difference is bounded directly, each side added only when
finite:

```math
\theta^{\min}_\ell \le \theta_{\mathrm{fr}(\ell),t} - \theta_{\mathrm{to}(\ell),t}
\le \theta^{\max}_\ell. \tag{22}
```

The window that the blocks actually see is enforced by the W-space tangent cuts,
with the limits clipped to ``\pm(\tfrac{\pi}{2}-10^{-6})`` before taking the
tangent:

```math
\tan\!\big(\theta^{\min}_\ell\big)\,W^{r}_{p(\ell),t}
\le W^{i}_{p(\ell),t} \le
\tan\!\big(\theta^{\max}_\ell\big)\,W^{r}_{p(\ell),t}. \tag{23}
```

No constraint ties ``\theta`` to the blocks, so (22) and the reference (24) shape
only this decoupled angle and do not tighten the relaxation. The relaxation itself
stays entirely in W-space, with its angle window set by the tangent cuts (23)
alone.

### Reference bus

```math
\theta_{\mathrm{ref},t} = 0, \qquad t\in\mathcal{T}. \tag{24}
```

### Generation constraints

Per-period envelopes (added only when strictly tighter than the scalar bound),
and renewable curtailment:

```math
p^g_{g,t} \le \overline{P}^{ts}_{g,t},\ g\in\mathcal{G}\setminus\mathcal{G}^c; \qquad
p^g_{g,t} \ge \underline{P}^{ts}_{g,t},\ g\in\mathcal{G}. \tag{25}
```

```math
\underline{Q}^{ts}_{g,t} \le Q^g_{g,t} \le \overline{Q}^{ts}_{g,t},
\qquad g\in\mathcal{G}. \tag{26}
```

```math
p^g_{g,t} + p^c_{g,t} = \overline{P}^{ts}_{g,t},
\qquad g\in\mathcal{G}^c. \tag{27}
```

### Storage constraints

Supported for SOCWR (`add_socwr_storage_constraints!`). The converter loss
``r|S|^2/V^2`` of the exact models ([AC Polar](@ref), Eqs. 34–35) is lifted the
same way the Jabr cone lifts the voltage products. A variable
``\mathrm{ccm}_{s,t} \ge |S|^2/W_{B(s)}`` (the rotated cone (30)) carries the
current magnitude squared, so the loss becomes linear, ``r_s\,\mathrm{ccm}``. Its
upper bound ``(\overline{S}_s/V^{\min}_{\mathrm{glob}})^2`` is ``|S|^2/V^2`` at the
worst case, maximum apparent power over minimum voltage. Charge/discharge
exclusivity is enforced with a **binary indicator** ``z_{s,t}`` rather than the
bilinear complementarity of the exact models.

**Energy balance** (``s\in\mathcal{S}``):

```math
E_{s,1} = E^0_s + \eta^c_s\,c_{s,1} - \frac{d_{s,1}}{\eta^d_s}; \qquad
E_{s,t} = E_{s,t-1} + \eta^c_s\,c_{s,t} - \frac{d_{s,t}}{\eta^d_s},\ t\ge 2. \tag{28}
```

**Charge / discharge exclusivity** (big-M on the indicator):

```math
c_{s,t} \le z_{s,t}\,C^{\max}_s, \qquad
d_{s,t} \le (1 - z_{s,t})\,D^{\max}_s. \tag{29}
```

**Current magnitude (rotated cone)** at the connection bus:

```math
(P^s_{s,t})^2 + (Q^s_{s,t})^2 \le W_{B(s),t}\,\mathrm{ccm}_{s,t}. \tag{30}
```

**Converter losses** (active and reactive, linear in the lifted current):

```math
P^s_{s,t} + \big(d_{s,t} - c_{s,t}\big) = P^{loss}_s + r_s\,\mathrm{ccm}_{s,t}. \tag{31}
```

```math
Q^s_{s,t} = Q^{sc}_{s,t} + Q^{loss}_s + x_s\,\mathrm{ccm}_{s,t}. \tag{32}
```

### Complete model

```math
\min_{\substack{\theta,\,W_{bb},\,W^{r},\,W^{i},\,p^g,\,Q^g,\,p^c,\\
E,\,c,\,d,\,P^s,\,Q^s,\,Q^{sc},\,z,\,\mathrm{ccm}}}
\quad F + F^{curt} \qquad (\text{Eq. } 11)
```

subject to:

| Block | Equations |
|---|---|
| Flows | (12)–(15) |
| Balance | (16), (17) |
| Jabr cone | (18) |
| LNC cuts | (19), (20) |
| Thermal (SOC) | (21) |
| Angle / reference | (22)–(24) |
| Generation | (25)–(27) |
| Storage | (28)–(32) |
| Bounds | (1)–(9) |

## What it captures / what it loses

!!! info "Philosophy"
    **Captures:** the *local* relation among the four blocks of each branch: the
    conic tie, plus the voltage/angle information re-injected by the LNC cuts.

    **Loses:** *cycle consistency* (KVL). The tie is written per bus pair. It does
    not coordinate the angles around a loop. In meshed networks the angles summed
    over a cycle may fail to close at zero, because each pair sets its
    ``\cos``/``\sin`` only to settle its own local account.

## Lirion vs. the original

!!! warning "Lirion vs. the original"
    Jabr's conic-programming formulation [Jabr2006](@cite) targets **radial
    distribution** load flow, where the relaxation is known to be exact under
    mild conditions. Lirion applies the same W-space second-order cone relaxation
    to **general (meshed) transmission** networks, where the per-pair tie no
    longer captures cycle consistency (see above) and a gap can open. The **LNC
    cuts** tighten it. Lirion also ships a smooth-NLP variant,
    `socwrnl`. It keeps the same relaxation but writes the Jabr tie as
    ``W_r^2 + W_i^2 \le W_{ii} W_{jj}`` and each thermal limit as
    ``p^2 + q^2 \le (\mathrm{rate}_a)^2``, smooth quadratic inequalities in place
    of explicit cones, so the model goes to a nonlinear interior-point solver
    (`Ipopt` by default) rather than a conic one. The flows, balance, LNC cuts,
    and the angle and storage constraints are all unchanged.

## See also

- Neighbours: [AC Rectangular (ACR)](@ref) (the flows are identical, only the
  status of the blocks changes) · [QC (Quadratic Convex)](@ref) (opens the
  blocks from the inside and reuses these LNC cuts) · [AC Polar](@ref).
- Original relaxation: radial conic programming, Jabr [Jabr2006](@cite).
- Code: constructor `socwr()`; smooth variant `socwrnl()`; source
  `src/modeling/branch/socwr.jl`.
