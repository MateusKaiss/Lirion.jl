# LPAC

The [DC](@ref) model delivered the cheapest possible solve, but at the price of an
amputation: no voltage magnitude, no reactive power, no losses. For transmission operation
under normal conditions this may suffice; for any study involving reactive support, voltage
limits, or a more heavily loaded network, the model simply does not see the phenomenon. The
LPAC model (*Linear Programming AC*) meets that gap with an intermediate question: can
voltage and reactive power be given back to the model **without** giving back the
non-convexity? The answer is yes, and the mechanism is a change of variable. Instead of
working with the magnitude ``v_i``, LPAC works with its deviation from nominal.

!!! info "Approximation, not relaxation"
    Like the [DC](@ref) model, LPAC is an **approximation**, not a relaxation: its optimum
    is neither a floor nor a ceiling on the true AC cost. But it errs far less, and over a
    much wider range of operating regimes.

## The voltage-deviation variable

For each bus, define

```math
\phi_i := v_i - 1,
```

the deviation of the magnitude from its nominal of ``1`` pu.

Why does this substitution help? In normal operation ``v_i`` lies between ``0.95`` and
``1.05``, so ``\phi_i`` is a small number, on the order of ``\pm 0.05``. And a small number
has a useful property: its square is negligible (``0.05^2 = 0.0025``, twenty times smaller).
It is this asymmetry that lets the flow expressions be linearized without significant error.
The nonlinearity of the magnitude products concentrates in the ``\phi^2`` terms, and they
vanish. Compare with the choice made by the [DC](@ref) model, which was more radical: there
``v_i = 1``, i.e. ``\phi_i = 0``, fixed. LPAC leaves ``\phi_i`` free as a variable and treats
only its products as negligible. One extra variable per bus, and reactive power returns to
the model. Reactive power rides on the voltage magnitudes, the very decoupling the [DC](@ref)
model invoked when it dropped ``q``. Giving the magnitude back as a free variable is what lets
reactive power act again.

## Derivation: three linearizations over the W-space flows

The structure that organizes the whole model is this. Recall the flows of the
[SOCWR (Jabr conic relaxation)](@ref) and [QC (Quadratic Convex)](@ref) models, which are
linear in the ``W`` blocks:

```math
p_{ij} = \frac{g + g^{fr}}{tm}\,w_{ii}
       + \frac{-g\,t_r + b\,t_i}{tm}\,w_r
       + \frac{-b\,t_r - g\,t_i}{tm}\,w_i .
```

LPAC reuses these formulas wholesale. It replaces each block by an expression linear in
``\phi`` and ``\theta``. Three substitutions, one per block.

This is the mirror image of the W-space relaxations. [SOCWR (Jabr conic relaxation)](@ref)
and [QC (Quadratic Convex)](@ref) kept the blocks as free variables and fenced them with a
cone or with envelopes, enlarging the feasible set to earn a lower bound. LPAC does the
opposite. It substitutes each block with an explicit expression in ``\phi``, ``\theta`` and
``cs``, so the block is no longer a free variable but a linear formula. Nothing is enlarged,
which is why LPAC is an approximation with no bound rather than a relaxation.

### Substitution 1: the self block ``w_{ii} = v_i^2``

```math
w_{ii} = v_i^2 = (1 + \phi_i)^2
= 1 + 2\phi_i + \underbrace{\phi_i^2}_{\text{negligible}}
\;\cong\; 1 + 2\phi_i .
```

### Substitution 2: the real cross block ``w_r = v_i v_j \cos\theta_{ij}``

First, the product of magnitudes:

```math
v_i v_j = (1 + \phi_i)(1 + \phi_j)
= 1 + \phi_i + \phi_j + \underbrace{\phi_i\phi_j}_{\text{negligible}}
\;\cong\; 1 + \phi_i + \phi_j .
```

Now introduce a new variable for the cosine, ``cs \cong \cos\theta_{ij}``:

```math
w_r \cong (1 + \phi_i + \phi_j)\,cs = cs + (\phi_i + \phi_j)\,cs .
```

Since ``cs`` is approximately ``1`` in the regime of interest (small angles), write
``cs = 1 + (cs - 1)`` with ``(cs - 1)`` small. The product ``(\phi_i + \phi_j)(cs - 1)`` is
small times small (second order) and vanishes:

```math
w_r \cong cs + \phi_i + \phi_j .
```

### Substitution 3: the imaginary cross block ``w_i = v_i v_j \sin\theta_{ij}``

```math
w_i \cong \underbrace{(1 + \phi_i + \phi_j)}_{\cong\, 1}\cdot
\underbrace{\sin\theta_{ij}}_{\cong\,\theta_{ij}}
\;\cong\; \theta_i - \theta_j .
```

Two approximations at once: the magnitude product drops to ``1`` (it is not worth keeping
the ``\phi`` here, because they would multiply ``\theta_{ij}``, itself already small, second
order again), and the sine becomes the angle, the same ``\sin x \cong x`` used by the
[DC](@ref) model.

### The resulting flows

Substituting the three into the W-space formulas yields the LPAC flows, linear in ``\phi``,
``\theta`` and ``cs``:

```math
p_{ij} = \frac{g + g^{fr}}{tm}(1 + 2\phi_i)
       + \frac{-g\,t_r + b\,t_i}{tm}(cs + \phi_i + \phi_j)
       + \frac{-b\,t_r - g\,t_i}{tm}\,\theta_{ij}. \tag{1}
```

```math
q_{ij} = -\frac{b + b^{fr}}{tm}(1 + 2\phi_i)
       - \frac{-b\,t_r - g\,t_i}{tm}(cs + \phi_i + \phi_j)
       + \frac{-g\,t_r + b\,t_i}{tm}\,\theta_{ij}. \tag{2}
```

```math
p_{ji} = (g + g^{to})(1 + 2\phi_j)
       + \frac{-g\,t_r - b\,t_i}{tm}(cs + \phi_i + \phi_j)
       - \frac{-b\,t_r + g\,t_i}{tm}\,\theta_{ij}. \tag{3}
```

```math
q_{ji} = -(b + b^{to})(1 + 2\phi_j)
       - \frac{-b\,t_r + g\,t_i}{tm}(cs + \phi_i + \phi_j)
       - \frac{-g\,t_r - b\,t_i}{tm}\,\theta_{ij}, \tag{4}
```

with ``\theta_{ij} = \theta_i - \theta_j``. Note what survived in contrast with the
[DC](@ref) model: the conductances ``g`` are present (so there are losses:
``p_{ij} \neq -p_{ji}``), the shunts are present, the complex tap is present, and ``q``
exists. The model keeps its AC structure; only the variables have been linearized. The
from-side coefficients all carry the ``1/tm`` of the complex tap, while the two to-side self
terms ``(g+g^{to})(1+2\phi_j)`` and ``-(b+b^{to})(1+2\phi_j)`` do not, exactly as in the
[SOCWR (Jabr conic relaxation)](@ref) flows they inherit (the tap sits on the *from* end).

## The cosine variable

The variable ``cs`` must behave like ``\cos\theta_{ij}``, but the equality
``cs = \cos\theta_{ij}`` is non-convex. It is the one nonlinearity the three substitutions
could not linearize away, the single trigonometric coupling left in the flows. How that
coupling is fenced decides whether the finished model is an LP or a QCQP, as the fidelity
note below explains. LPAC fences it two ways.

**Box bounds.** Exactly the three cases of the [QC (Quadratic Convex)](@ref) model, decided
by the position of the angle interval (the cosine is not monotone if the interval crosses
zero):

- If ``\theta^{\min} \ge 0``: ``\quad \cos\theta^{\max} \le cs \le \cos\theta^{\min}``.
- If ``\theta^{\max} \le 0``: ``\quad \cos\theta^{\min} \le cs \le \cos\theta^{\max}``.
- If it crosses zero: ``\quad \min(\cos\theta^{\min}, \cos\theta^{\max}) \le cs \le 1``.

**Quadratic cap.** And here the parabola of the [QC (Quadratic Convex)](@ref) model
reappears:

```math
cs \le 1 - \frac{1 - \cos m}{m^2}\,(\theta_i - \theta_j)^2,
\qquad m = \max(|\theta^{\min}|, |\theta^{\max}|).
```

## Constraints

**Voltage limits.** Directly on the deviation:

```math
V_i^{\min} - 1 \le \phi_i \le V_i^{\max} - 1. \tag{5}
```

**Thermal limit.** The cone returns, because reactive power is back (in the DC model the
limit was a box in ``p``):

```math
p_{ij}^2 + q_{ij}^2 \le \overline{S}_\ell^2,
\qquad p_{ji}^2 + q_{ji}^2 \le \overline{S}_\ell^2. \tag{6}
```

**Angular limits.**

```math
\theta^{\min} \le \theta_i - \theta_j \le \theta^{\max}. \tag{7}
```

## LPAC Formulation

The LPAC constructor (`lpac()`, model type `LPAC`) assembles the standard skeleton (buses,
branches, generators, storage, active/reactive balance, objective) and builds the branch
block from the linearized W-space flows derived above. With the cosine held by the parabola
and each thermal limit by a second-order cone, the problem is a **QCQP**, which Lirion solves
with `Ipopt` by default. The flow equations and the voltage, thermal and angle constraints
are Eqs. (1)–(7). This section gathers the complete problem around them (sets, parameters,
variables, objective, balance, and the remaining blocks) and continues the numbering
from (8).

### Sets

- ``\mathcal{N}``: set of buses.
- ``\mathcal{N}^{nv} \subseteq \mathcal{N}``: non-virtual buses, on which the power balance is
  enforced (`non_virtual_bus_ids`).
- ``\mathcal{L}``: set of active branches (in-service lines and transformers).
- ``\mathcal{G}``: set of existing generators.
- ``\mathcal{G}^c \subseteq \mathcal{G}``: subset of curtailable generators (renewables).
- ``\mathcal{S}``: set of storage units.
- ``\mathcal{T} = \{1,\dots,T\}``: set of time periods in the horizon.

### Indices and incidence functions

- ``b, i, j \in \mathcal{N}``: bus indices.
- ``\ell \in \mathcal{L}``: branch index.
- ``g \in \mathcal{G}``: generator index.
- ``s \in \mathcal{S}``: storage index.
- ``t \in \mathcal{T}``: time index.
- ``f(\ell),\ t(\ell)``: *from* and *to* buses of branch ``\ell`` (code: `br.from[l]`,
  `br.to[l]`).
- ``B(s)``: connection bus of storage ``s``.

### Parameters

**Branches.**

- ``g_\ell,\ b_\ell``: series conductance and susceptance of branch ``\ell``.
- ``g^{fr}_\ell,\ b^{fr}_\ell``: shunt conductance and susceptance, *from* side.
- ``g^{to}_\ell,\ b^{to}_\ell``: shunt conductance and susceptance, *to* side.
- ``t_{r,\ell},\ t_{i,\ell}``: real and imaginary parts of the complex tap
  ``T_\ell = \tau_\ell e^{j\Phi_\ell}`` (code: `br.tap_real`, `br.tap_imag`).
- ``tm_\ell = |T_\ell|^2 = \tau_\ell^2``: squared tap magnitude; the code multiplies the
  coefficients below by ``1/tm_\ell`` via `br.inv_tap[l]^2`.
- ``\overline{S}_\ell``: apparent-flow limit of the branch (`rate_a`).
- ``\theta_\ell^{\min},\ \theta_\ell^{\max}``: angular aperture limits (code: `br.ang_min`,
  `br.ang_max`; `NaN` means absent).

**Buses.**

- ``V_b^{\min},\ V_b^{\max}``: voltage magnitude limits at bus ``b``.
- ``g^{sh}_b,\ b^{sh}_b``: shunt conductance and susceptance of the bus.
- ``P^d_{b,t},\ Q^d_{b,t}``: active and reactive demand at bus ``b``.

**Generators.**

- ``P_g^{\min},\ P_g^{\max}``: static active-power limits.
- ``\overline{P}^{ts}_{g,t},\ \underline{P}^{ts}_{g,t}``: active-power time envelope.
- ``Q_g^{\min},\ Q_g^{\max}``: static reactive-power limits.
- ``\overline{Q}^{ts}_{g,t},\ \underline{Q}^{ts}_{g,t}``: reactive-power time envelope.
- ``c_{g,t,k}``: coefficients of the cost polynomial of ``g``.
- ``c^{curt}_g``: unit curtailment cost of ``g``.

**Storage.** LPAC reuses the **DC-style storage block** (`add_dc_storage_constraints!`), so
the converter-loss parameters ``r_s,x_s,Q^{loss}_s`` of the [AC Polar](@ref) model are
**not** used here.

- ``E_s^{\min},\ E_s^{\max}``: stored-energy limits; ``E^0_s``: initial energy.
- ``C_s^{\max},\ D_s^{\max}``: maximum charge and discharge power.
- ``\eta^c_s,\ \eta^d_s``: charge and discharge efficiencies.
- ``P^{loss}_s``: fixed active loss; ``\overline{S}_s``: inverter thermal limit.
- ``Q_s^{\min},\ Q_s^{\max}``: reactive-power limits of the storage.

### Decision variables

**Voltage.**

```math
\theta_{b,t} \in \mathbb{R} \;\text{: voltage angle at bus } b \text{ at time } t. \tag{8}
```

```math
\phi_{b,t} \in [V_b^{\min}-1,\ V_b^{\max}-1] \;\text{: voltage-magnitude deviation}
\;(v_b = 1 + \phi_b),\quad \text{see Eq. (5)}.
```

**Cosine surrogate.** One ``cs_{\ell,t}`` per branch, box-bounded by the angle regime
(exactly the three [QC (Quadratic Convex)](@ref) cases; a fourth, default branch applies
when the angle limits are absent):

```math
cs_{\ell,t} \in
\begin{cases}
[\cos\theta_\ell^{\max},\ \cos\theta_\ell^{\min}], & \theta_\ell^{\min}\ge 0,\\[2pt]
[\cos\theta_\ell^{\min},\ \cos\theta_\ell^{\max}], & \theta_\ell^{\max}\le 0,\\[2pt]
[\min(\cos\theta_\ell^{\min},\cos\theta_\ell^{\max}),\ 1], & \theta_\ell^{\min}<0<\theta_\ell^{\max},\\[2pt]
[-1,\ 1], & \theta_\ell^{\min},\theta_\ell^{\max}\ \text{absent (NaN)}.
\end{cases}
\tag{9}
```

**Generation.**

```math
p^g_{g,t} \in [P_g^{\min}, P_g^{\max}] \;\text{: active power generated by } g. \tag{10}
```

```math
Q^g_{g,t} \in [Q_g^{\min}, Q_g^{\max}] \;\text{: reactive power generated by } g. \tag{11}
```

```math
p^c_{g,t} \in [0, \overline{P}^{ts}_{g,t}] \;\text{: curtailed renewable power, } g \in \mathcal{G}^c. \tag{12}
```

**Storage** (only when ``\mathcal{S}\neq\emptyset``): ``E_{s,t}\in[E_s^{\min},E_s^{\max}]``,
``c_{s,t}\in[0,C_s^{\max}]``, ``d_{s,t}\in[0,D_s^{\max}]``,
``P^s_{s,t}\in[-\overline{S}_s,\overline{S}_s]``, ``Q^s_{s,t}\in[Q_s^{\min},Q_s^{\max}]``, and
``Q^{sc}_{s,t}\in[Q_s^{\min},Q_s^{\max}]`` (the converter-compensation variable `qsc` is
created but is **inert** under the DC-style storage block, as noted below).

### Objective function

Identical to the [AC Polar](@ref) model: fuel cost on the power in physical units
(``S^{base}\,p^g``) plus renewable curtailment cost.

```math
F = \sum_{t\in\mathcal{T}} \sum_{g\in\mathcal{G}} \sum_{k=1}^{m_{g,t}}
c_{g,t,k}\left(S^{base}\,p^g_{g,t}\right)^{m_{g,t} - k}. \tag{13}
```

```math
F^{curt} = \sum_{g\in\mathcal{G}^c} \sum_{t\in\mathcal{T}}
c^{curt}_g\,S^{base}\,p^c_{g,t}, \qquad \min\; F + F^{curt}. \tag{14}
```

### Branch flows

The four linearized W-space flows ``P^f_{\ell,t}, Q^f_{\ell,t}, P^t_{\ell,t},
Q^t_{\ell,t}`` are Eqs. (1)–(4), linear in ``\phi``, ``cs`` and ``\theta``, with
``\theta_{ij} = \theta_{f(\ell)} - \theta_{t(\ell)}`` and all *from*-side coefficients scaled
by ``1/tm_\ell``. The implementation reifies these four flows as explicit variables bounded to
``[-\overline{S}_\ell, \overline{S}_\ell]``, each tied by equality to its expression, so the
thermal cone (6) is a second-order cone in two named variables per side. Those box bounds are
already implied by (6); they condition the solver and add no law to the model.

### Power balance

**Active power.** The bus shunt enters as the linearized ``g^{sh}_b(1 + 2\phi_b)``:

```math
\sum_{g\in\mathcal{G}(b)} p^g_{g,t}
- \sum_{s\in\mathcal{S}(b)} P^s_{s,t}
- P^d_{b,t}
- \sum_{\ell:\,f(\ell)=b} P^f_{\ell,t}
- \sum_{\ell:\,t(\ell)=b} P^t_{\ell,t}
- g^{sh}_b\,(1 + 2\phi_{b,t}) = 0,
\qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{15}
```

**Reactive power.** The bus shunt enters as ``+\,b^{sh}_b(1 + 2\phi_b)``:

```math
\sum_{g\in\mathcal{G}(b)} Q^g_{g,t}
- \sum_{s\in\mathcal{S}(b)} Q^s_{s,t}
- Q^d_{b,t}
- \sum_{\ell:\,f(\ell)=b} Q^f_{\ell,t}
- \sum_{\ell:\,t(\ell)=b} Q^t_{\ell,t}
+ b^{sh}_b\,(1 + 2\phi_{b,t}) = 0,
\qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{16}
```

### Cosine envelope

Besides the box (9), ``cs`` is capped by the [QC (Quadratic Convex)](@ref) parabola, added
only when ``m_\ell = \max(|\theta_\ell^{\min}|,|\theta_\ell^{\max}|) > 0``:

```math
cs_{\ell,t} \le 1 - \frac{1 - \cos m_\ell}{m_\ell^2}\,
\big(\theta_{f(\ell),t} - \theta_{t(\ell),t}\big)^2,
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{17}
```

### Branch thermal and angle limits

The second-order cone (6) on each side, and the angle-aperture limits (7), each bound imposed
only when the corresponding ``\theta_\ell^{\min}/\theta_\ell^{\max}`` is present (not `NaN`).

### Reference bus

```math
\theta_{ref,t} = 0, \qquad t\in\mathcal{T}. \tag{18}
```

### Generation constraints

Time envelopes (added only when strictly tighter than the static bounds) and the renewable
curtailment identity. The per-period upper active bound is skipped for renewables, which take
the curtailment identity instead; the lower active bound and the reactive envelope apply to
every generator:

```math
p^g_{g,t} \le \overline{P}^{ts}_{g,t}\ (g\notin\mathcal{G}^c), \qquad
\underline{P}^{ts}_{g,t} \le p^g_{g,t}\ (g\in\mathcal{G}), \qquad
\underline{Q}^{ts}_{g,t} \le Q^g_{g,t} \le \overline{Q}^{ts}_{g,t},
\qquad t\in\mathcal{T}. \tag{19}
```

```math
p^g_{g,t} + p^c_{g,t} = \overline{P}^{ts}_{g,t},
\qquad g\in\mathcal{G}^c,\ t\in\mathcal{T}. \tag{20}
```

### Storage constraints (DC-style)

When ``\mathcal{S}\neq\emptyset`` the shared DC storage block is added. It differs from the
[AC Polar](@ref) storage in that the active injection carries **only the fixed loss**
``P^{loss}_s`` (no ``r_s|S|^2/V^2`` converter term) and there is **no reactive-coupling
equation** (the `qsc` variable is created but unused):

```math
E_{s,1} = E^0_s + \eta^c_s\,c_{s,1} - \frac{d_{s,1}}{\eta^d_s}, \qquad
E_{s,t} = E_{s,t-1} + \eta^c_s\,c_{s,t} - \frac{d_{s,t}}{\eta^d_s}\ (t\ge 2). \tag{21}
```

```math
c_{s,t}\,d_{s,t} = 0, \qquad
(P^s_{s,t})^2 + (Q^s_{s,t})^2 \le \overline{S}_s^2, \qquad
P^s_{s,t} + d_{s,t} - c_{s,t} = P^{loss}_s. \tag{22}
```

### Complete model

```math
\min_{\substack{\theta,\,\phi,\,cs,\,p^g,\,Q^g,\,p^c,\\ E,\,c,\,d,\,P^s,\,Q^s,\,Q^{sc}}}
\quad F + F^{curt} \qquad (\text{Eq. } 14)
```

subject to:

| Block | Equations |
|---|---|
| Flows | (1)–(4) |
| Balance | (15), (16) |
| Cosine surrogate | (9), (17) |
| Thermal / angle | (6), (7) |
| Reference bus | (18) |
| Generation | (19), (20) |
| Storage | (21), (22) |
| Bounds | (5), (8), (10)–(12) |

## What it captures / what it loses

!!! info "Philosophy"
    **Captures:** voltage magnitude (through the deviation ``\phi``), reactive power ``q``,
    ohmic losses (``p_{ij} \neq -p_{ji}``, since the conductances ``g`` survive), line
    shunts, and the complex transformer tap. In form the model is still AC. **Loses:** the
    exact cosine (replaced by the ``cs`` box and parabolic cap) and all second-order terms
    (``\phi^2``, ``\phi_i\phi_j``, ``(\phi_i+\phi_j)(cs-1)``, and the
    ``(\phi_i+\phi_j)\theta_{ij}`` of the imaginary block). Being an approximation, its optimum
    is neither a floor nor a ceiling on the true AC cost. But it errs far less than the DC
    model, and over a far wider set of regimes.

!!! note "When to trust it"
    LPAC is at its best where its linearization holds and no bound is required: a transmission
    network near nominal voltage, where ``\phi`` stays small enough that ``\phi^2`` is
    genuinely negligible, and angle differences stay modest so the cosine surrogate and
    ``\sin\theta \cong \theta`` both stay accurate. There it recovers what the [DC](@ref)
    model cannot see, the voltage magnitude, reactive power and losses, while staying convex
    and cheap to solve. Prefer it over DC whenever reactive support, voltage limits, or losses
    are part of the question, or whenever the network is loaded past the near-nominal regime
    the DC model assumes.

    It stops being reliable when ``\phi`` grows large, under heavy loading or voltage stress,
    where ``\phi^2`` is no longer negligible, and when angle differences widen past the range
    where the cosine surrogate is tight. When the answer must be a certified bound rather than
    an estimate, use a relaxation instead ([SOCWR (Jabr conic relaxation)](@ref),
    [QC (Quadratic Convex)](@ref), [SDPWRM](@ref)). When it must be exact, the [AC Polar](@ref)
    model.

## Lirion vs. the original

!!! warning "Fidelity note"
    The original LPAC article of Coffrin and Van Hentenryck (2014) fences the cosine with a
    **piecewise-linear approximation, a bundle of tangent lines**,
    ``cs \le \cos\hat\theta - \sin\hat\theta\,(\theta_{ij} - \hat\theta)``, one for each chosen
    angle point ``\hat\theta``. This keeps the model strictly **LP**, which is the reason for
    the name. Lirion instead implements the cosine with the **parabola**, more compact (one
    constraint instead of a bundle) and tighter near the centre, at the cost of leaving LP:
    with the parabola the model is a **QCQP**. The same divergence appears in the thermal
    limit. Lirion keeps the second-order **cone**. The original article writes the
    apparent-power limit as the same quadratic circle, an optional added constraint, and only
    remarks that it can be replaced by piecewise-linear cuts (a **polygon**) to stay inside
    LP. It does not build that polygon itself.

## See also

- Neighbours: [DC](@ref) (the amputation LPAC answers) · [NFA](@ref) (the other direction
  from the DC model) · exact reference [AC Polar](@ref).
- The W-space flows reused here: [SOCWR (Jabr conic relaxation)](@ref) ·
  [QC (Quadratic Convex)](@ref) (box cases and cosine parabola).
- Foundations: [Foundations](@ref).
- Reference implementation: constructor `lpac()`; source `src/modeling/branch/lpac.jl`.
