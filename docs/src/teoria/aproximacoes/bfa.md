# BFA (LinDistFlow)

The **branch-flow approximation** (BFA) is the last rung of the branch-flow family. It stands to
the branch-flow second-order-cone model (SOCBF) much as the [DC](@ref) model stands to the AC
model: a linearization of the model one level up. The step is narrow. One quantity leaves, the
squared branch current, and with it the only nonlinearity in the branch physics. In the
distribution literature the model is known as **LinDistFlow** [BaranWu1989](@cite), and it is
probably the most widely used linear network model in distribution optimization, for Volt/Var
control, distributed-generation siting, and feeder reconfiguration.
See [Kersting2012](@cite) for the distribution-system modelling context.

!!! info "Approximation, not relaxation"
    Like the [DC](@ref) and [LPAC](@ref) models, the BFA is an **approximation, not a
    relaxation**: its optimum is neither a floor nor a ceiling of the true cost.

## Derivation

### The branch-flow starting point

The branch-flow view describes each branch ``\ell`` by the flows ``p_{ij}, q_{ij}`` that enter it
at its *from* bus ``i`` and the flows ``p_{ji}, q_{ji}`` that enter at its *to* bus ``j``,
together with the squared magnitude of the current it carries, ``|I_\ell|^2``. Each bus holds a
single voltage state, the squared magnitude ``w_b = |V_b|^2``. Two exact relations tie a branch
to its two end buses: the power it loses, and the drop in squared voltage from one end to the
other. On a branch with no shunt admittance, where the algebra stays readable, they are

```math
p_{ij} + p_{ji} = r_\ell\,|I_\ell|^2,
\qquad
\frac{w_i}{tm_\ell} - w_j
  = 2\,(r_\ell\,p_{ij} + x_\ell\,q_{ij}) - (r_\ell^2 + x_\ell^2)\,|I_\ell|^2,
```

with ``tm_\ell`` the squared tap magnitude. The current magnitude has its own defining relation,
``|I_\ell|^2 = (p_{ij}^2 + q_{ij}^2)\big/(w_i/tm_\ell)``, which the SOCBF model carries as a
rotated second-order cone. That cone is the one nonlinearity among the branch relations, and it
is what the approximation targets. It is a different object from the apparent-power thermal limit
the formulation keeps, which stays a cone further down.

### One term fewer

With the branch-flow relations in hand, the approximation is a single one:

```math
|I_\ell|^2 \cong 0,
```

the squared branch current neglected wherever it appears. This does not mean the branch carries
no current. The flows stay; the only things dropped are the small terms built from
``|I_\ell|^2``, and with them the defining relation for ``|I_\ell|^2`` itself. The rest of the
derivation is a matter of seeing what survives in each equation.

!!! note "What justifies the discard"
    Compare the two terms of the voltage-drop equation:

    ```math
    \underbrace{2\,(r\,p_{ij} + x\,q_{ij})}_{\text{principal term}}
    \quad\text{versus}\quad
    \underbrace{(r^2 + x^2)\,|I_\ell|^2}_{\text{discarded term}}.
    ```

    The first is linear in the impedance, the second quadratic (``r^2+x^2``). In per unit,
    branch impedances are small numbers, a typical feeder having ``r`` in the range ``0.01`` to
    ``0.1``, so the square is an order of magnitude smaller. It is the same logic as the
    [LPAC](@ref) model discarding ``\phi^2``: keep the first order, throw away the second.

### Quantifying the error

The neglected ``|I_\ell|^2`` sits in the branch losses (the resistive ``r_\ell\,|I_\ell|^2`` and
the reactive ``x_\ell\,|I_\ell|^2``) and in the quadratic ``(r_\ell^2+x_\ell^2)\,|I_\ell|^2`` of
the voltage drop. With ``|I_\ell|^2 \cong (p^2+q^2)/w`` (``w`` the squared voltage magnitude) and
the voltage near ``1``, these sit at order ``r\cdot p`` in per unit relative to the terms they
correct, for a line with ``r\cong x`` and ``p\cong q``. The error therefore grows with the
resistance and with the loading:

- on a line with ``r=0.05`` and a flow of ``0.3`` pu, the neglected part is ``\cong 1.5\%`` of
  the drop;
- at a flow of ``1.0`` pu it is ``\cong 5\%``.

In the load-dominated regime these feeders normally run, the bias has a definite sign. Dropping
the branch losses lets each line carry only its net load, without the extra current the losses
really draw, so its flows run short of the true ones and the computed drop comes out low: **the
BFA underestimates the voltage drop**, and its voltages sit a little high. The quadratic term in
the drop, which is subtracted there, would on its own lengthen the computed drop and so pulls the
other way, but it offsets only half of the loss effect and cannot reverse the sign. Under heavy
distributed generation, where power flows back up the feeder, the sign can turn over with it.
This load-dominated bias is the familiar behaviour of LinDistFlow.

### The surviving branch relations

Besides the current-magnitude definition the starting point gave, the SOCBF ties each branch to
its end buses through the power losses and the voltage drop. Setting ``|I_\ell|^2 \cong 0``
removes the current-magnitude definition outright (there is nothing left to define) and trims the
other two. Restoring the line charging the starting point set aside, here is what the losses and
the drop keep.

**Losses.** The active- and reactive-power sums across the branch keep only their shunt terms,
the series losses going with ``|I_\ell|^2``:

```math
p_{ij} + p_{ji} = r\,|I_\ell|^2 + \frac{g^{fr}}{tm}\,w_i + g^{to}\,w_j
\;\longrightarrow\;
p_{ij} + p_{ji} = \frac{g^{fr}}{tm}\,w_i + g^{to}\,w_j,
```

```math
q_{ij} + q_{ji} = x\,|I_\ell|^2 - \frac{b^{fr}}{tm}\,w_i - b^{to}\,w_j
\;\longrightarrow\;
q_{ij} + q_{ji} = -\frac{b^{fr}}{tm}\,w_i - b^{to}\,w_j.
```

Here ``w_i`` and ``w_j`` are the squared voltage magnitudes ``|V_i|^2`` and ``|V_j|^2`` at the
*from* and *to* buses (the LinDistFlow voltage state, one value per bus), and ``tm = |T|^2`` is
the squared tap magnitude, as on the [AC Polar](@ref) and [LPAC](@ref) pages. The series losses
``r\,|I_\ell|^2`` and ``x\,|I_\ell|^2`` become zero, but the shunt losses stay. This contrasts
with the [DC](@ref) model. There ``p_{ij} = -p_{ji}`` exactly, the model was rigorously lossless,
whereas here ``p_{ij} + p_{ji} \neq 0`` in general, because the shunt still drains power, a drain
proportional to ``|V|^2``. In practice that drain is small: the shunt conductance ``g^{fr},
g^{to}`` is usually near zero (line charging is almost purely reactive), so the kept loss is tiny
and the BFA runs close to lossless, which is part of why its voltages come out high.

**Voltage drop.** Setting ``|I_\ell|^2 \cong 0`` removes the quadratic
``(r^2+x^2)\,|I_\ell|^2`` from the Ohm's-law relation. With the line charging back in, the exact
relation also carries a coefficient ``1 + 2(r\,g^{fr} - x\,b^{fr})`` on ``w_i/tm``; the BFA rounds
it to ``1``. This is a further, smaller simplification, built from the charging shunt rather than
the current, and it is identically zero on a branch without charging, which is why the shunt-free
derivation above saw only the single discard. What is left is a linear tie between the two bus
voltages:

```math
\frac{w_i}{tm} - w_j = 2\,(r\,p_{ij} + x\,q_{ij}),
```

the classical LinDistFlow drop: the squared magnitude falls by twice the real part of the
impedance times the flow.

## BFA Formulation

The BFA constructor (`bfa()`, model type `BFA`) assembles the standard skeleton (buses,
branches, generators, storage, active/reactive balance, objective) in **W-space**: the only
voltage variable is the squared magnitude ``w_b = |V_b|^2`` per bus. This is the per-bus voltage
state the W-space relaxations also use, but the BFA ties buses together through the linear
voltage drop (8) alone, not through the off-diagonal bus-pair products and Jabr cone of
[SOCWR (Jabr conic relaxation)](@ref). The branch relations above are linear, but the thermal
limit is kept as a second-order cone (see the fidelity note below), so the problem is an
**SOCP**, which Lirion solves with `Clarabel` by default. The cone raises the solver class, not
the model's standing: it is an exact operating limit, not a relaxed physical law, so the BFA
stays an approximation and its optimum is still no bound.

### Sets

- ``\mathcal{N}``: set of buses.
- ``\mathcal{N}^{nv} \subseteq \mathcal{N}``: non-virtual buses, on which the power balance is
  enforced (`non_virtual_bus_ids`).
- ``\mathcal{L}``: set of active branches.
- ``\mathcal{G}``: generators; ``\mathcal{G}^c \subseteq \mathcal{G}``: curtailable (renewable)
  generators.
- ``\mathcal{S}``: storage units.
- ``\mathcal{T} = \{1,\dots,T\}``: time periods.

### Indices and incidence functions

- ``b, i, j \in \mathcal{N}``; ``\ell \in \mathcal{L}``; ``g \in \mathcal{G}``;
  ``s \in \mathcal{S}``; ``t \in \mathcal{T}``.
- ``f(\ell),\ t(\ell)``: *from* and *to* buses of branch ``\ell`` (code: `br.from[l]`,
  `br.to[l]`).
- ``B(s)``: connection bus of storage ``s``.

### Parameters

- ``r_\ell,\ x_\ell``: series resistance and reactance of branch ``\ell``.
- ``g^{fr}_\ell,\ b^{fr}_\ell``; ``g^{to}_\ell,\ b^{to}_\ell``: shunt conductance and
  susceptance of the *from* and *to* ends.
- ``tm_\ell = |T_\ell|^2 = \tau_\ell^2``: squared tap magnitude; the code writes the tap
  denominator as ``\texttt{inv\_tap}^2 = 1/tm_\ell``.
- ``\overline{S}_\ell``: apparent-power (thermal) limit of the branch (`rate_a`).
- ``\theta^{\min}_\ell,\ \theta^{\max}_\ell``: angular-aperture limits (`ang_min`, `ang_max`;
  `NaN` means absent).
- ``V_b^{\min},\ V_b^{\max}``: voltage-magnitude limits; ``g^{sh}_b,\ b^{sh}_b``: bus shunt
  conductance and susceptance; ``P^d_{b,t},\ Q^d_{b,t}``: active and reactive demand.
- Generation and **DC-style storage** parameters exactly as in the [LPAC](@ref) model.

### Decision variables

```math
w_{b,t} \in [(V_b^{\min})^2,\ (V_b^{\max})^2]
  \;\text{: squared voltage magnitude } |V_b|^2\ (\texttt{v\_magnitude\_sq}). \tag{1}
```

```math
P^f_{\ell,t},\ Q^f_{\ell,t},\ P^t_{\ell,t},\ Q^t_{\ell,t} \in [-\overline{S}_\ell,\ \overline{S}_\ell]
  \;\text{: from- and to-end active and reactive flows.} \tag{2}
```

```math
p^g_{g,t} \in [P_g^{\min}, P_g^{\max}],\quad
Q^g_{g,t} \in [Q_g^{\min}, Q_g^{\max}],\quad
p^c_{g,t} \in [0, \overline{P}^{ts}_{g,t}]\,(g\in\mathcal{G}^c)
\;\text{: generation and curtailment.} \tag{3}
```

Storage, when ``\mathcal{S}\neq\emptyset``, adds the DC-style block
``E_{s,t}, c_{s,t}, d_{s,t}, P^s_{s,t}, Q^s_{s,t}, Q^{sc}_{s,t}`` of the [LPAC](@ref) model.
There is no squared-current variable: the ``|I_\ell|^2`` of SOCBF is exactly what the
approximation removed.

### Objective function

Identical to the [AC Polar](@ref) model: fuel cost on power in physical units plus renewable
curtailment cost.

```math
F = \sum_{t\in\mathcal{T}}\sum_{g\in\mathcal{G}}\sum_{k=1}^{m_{g,t}}
  c_{g,t,k}\,(S^{base}\,p^g_{g,t})^{m_{g,t}-k},\qquad
F^{curt} = \sum_{g\in\mathcal{G}^c}\sum_{t\in\mathcal{T}} c^{curt}_g\,S^{base}\,p^c_{g,t}. \tag{4}
```

```math
\min\; F + F^{curt}. \tag{5}
```

### Branch relations

The losses carry only the shunt terms, and the voltage drop is linear in the flows (the two
relations derived above), per branch ``\ell`` and time ``t``:

```math
P^f_{\ell,t} + P^t_{\ell,t} = \frac{g^{fr}_\ell}{tm_\ell}\,w_{f(\ell),t} + g^{to}_\ell\,w_{t(\ell),t}, \tag{6}
```

```math
Q^f_{\ell,t} + Q^t_{\ell,t} = -\frac{b^{fr}_\ell}{tm_\ell}\,w_{f(\ell),t} - b^{to}_\ell\,w_{t(\ell),t}, \tag{7}
```

```math
\frac{w_{f(\ell),t}}{tm_\ell} - w_{t(\ell),t} = 2\,(r_\ell\,P^f_{\ell,t} + x_\ell\,Q^f_{\ell,t}). \tag{8}
```

### Branch thermal limits

A second-order cone on each end (not a box), the one nonlinearity the model keeps. It is an
operating limit, not a branch relation, and a different cone from the current-magnitude one the
approximation dropped:

```math
(P^f_{\ell,t})^2 + (Q^f_{\ell,t})^2 \le \overline{S}_\ell^2,\qquad
(P^t_{\ell,t})^2 + (Q^t_{\ell,t})^2 \le \overline{S}_\ell^2. \tag{9}
```

### Angle limits

The angular aperture is enforced in branch-flow form, as a pair of linear cuts on the flows and
the from-bus voltage (reused from SOCBF, `_add_socbf_angle_constraints!`), posted per branch
only where ``\theta^{\min}_\ell`` or ``\theta^{\max}_\ell`` is set:

```math
\tan(\theta^{\min}_\ell)\,\mathrm{Re}_{\ell,t} \le \mathrm{Im}_{\ell,t} \le
\tan(\theta^{\max}_\ell)\,\mathrm{Re}_{\ell,t}, \tag{10}
```

where ``\mathrm{Re}_{\ell,t}`` and ``\mathrm{Im}_{\ell,t}`` are the real and imaginary parts of a
branch-flow angle surrogate, each linear in ``w_{f(\ell),t}``, ``P^f_{\ell,t}`` and
``Q^f_{\ell,t}``. Their ratio ``\mathrm{Im}/\mathrm{Re}`` stands in for the tangent of the angle
across the branch, so holding it between ``\tan\theta^{\min}`` and ``\tan\theta^{\max}`` bounds
the aperture. With no angle variable in W-space, these two linear combinations are what the model
carries in place of ``\theta_{f(\ell)} - \theta_{t(\ell)}``.

### Power balance

At each non-virtual bus the shunt rides on the squared magnitude ``w_{b,t}`` (not a constant, as
in the [DC](@ref) and [NFA](@ref) models, because here the voltage is a variable):

```math
\sum_{g\in\mathcal{G}(b)} p^g_{g,t}
  - \sum_{s\in\mathcal{S}(b)} P^s_{s,t}
  - P^d_{b,t}
  - \sum_{\ell:\,f(\ell)=b} P^f_{\ell,t}
  - \sum_{\ell:\,t(\ell)=b} P^t_{\ell,t}
  - g^{sh}_b\,w_{b,t} = 0,
\qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{11}
```

```math
\sum_{g\in\mathcal{G}(b)} Q^g_{g,t}
  - \sum_{s\in\mathcal{S}(b)} Q^s_{s,t}
  - Q^d_{b,t}
  - \sum_{\ell:\,f(\ell)=b} Q^f_{\ell,t}
  - \sum_{\ell:\,t(\ell)=b} Q^t_{\ell,t}
  + b^{sh}_b\,w_{b,t} = 0,
\qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{12}
```

There is no reference-voltage equation. The magnitudes are held by their bounds (1) and tied to
one another by the voltage drop (8), and with no angle variable anywhere in the model there is no
angle reference to set.

### Generation and storage

The active and reactive generation envelopes (applied as the tighter of the static and
per-period bounds) and the renewable curtailment identity are those of the [LPAC](@ref) model:

```math
p^g_{g,t} + p^c_{g,t} = \overline{P}^{ts}_{g,t}, \qquad g\in\mathcal{G}^c,\ t\in\mathcal{T}. \tag{13}
```

Storage, when present, is the shared **DC-style block** (`add_dc_storage_constraints!`): energy
balance, a binary-free charge/discharge exclusivity, and the converter limit, exactly as in the
[DC](@ref) and [LPAC](@ref) models. The reactive storage variables are created but stay
decoupled from the DC-style block.

### Complete model

```math
\min_{\substack{w,\,P^f,\,Q^f,\,P^t,\,Q^t,\,p^g,\,Q^g,\,p^c,\\
  E,\,c,\,d,\,P^s,\,Q^s,\,Q^{sc}}}
\quad F + F^{curt} \qquad (\text{Eq. } 5)
```

subject to:

| Block | Equations |
|---|---|
| Branch losses | (6), (7) |
| Voltage drop | (8) |
| Thermal (cone) | (9) |
| Angle cuts | (10) |
| Power balance | (11), (12) |
| Generation | (13) |
| Storage (DC-style) | see [DC](@ref) |
| Bounds | (1)–(3) |

## What it captures / what it loses

!!! info "Philosophy"
    **Captures:** the voltage-magnitude profile through the linear drop (8), reactive power,
    the shunt losses, and the apparent-power thermal limit as an exact cone. On a radial feeder
    near nominal voltage the result is close to the full AC solution.

    **Loses:** the resistive (Joule) loss ``r\,|I_\ell|^2`` and the quadratic voltage-drop term
    ``(r^2+x^2)\,|I_\ell|^2``, both dropped with ``|I_\ell|^2``. With the losses gone the branch
    flows carry only the net load, so the model **underestimates the voltage drop** and runs the
    voltages a little high; the error grows with resistance and loading. Being an approximation,
    its optimum is neither a floor nor a ceiling on the true AC cost.

!!! note "When to trust it"
    The BFA works best on radial or lightly meshed distribution feeders near nominal voltage,
    where the squared current is small. That is the regime LinDistFlow was built for. It stops
    being reliable as resistance and loading rise, where the dropped ``|I_\ell|^2`` terms grow
    and the voltage drop is increasingly underestimated. When a certified bound is needed rather
    than an estimate, use the SOCBF relaxation it linearizes; for the exact answer, the
    [AC Polar](@ref) model.

## Lirion vs. the original

!!! warning "Fidelity note"
    The textbook LinDistFlow of Baran and Wu [BaranWu1989](@cite) is a **linear** model, written
    for radial feeders and dropping the branch charging. It carries no apparent-power term of its
    own, so linear-programming implementations either leave the thermal limit out or bound the
    apparent power with a linear box or polygon. Lirion keeps the two features that make its
    branch-flow family uniform, so its BFA is slightly richer than the classical model and is
    **not** a pure LP:

    - the thermal limit stays the exact apparent-power **cone** (9) rather than a box or polygon,
      the same apparent-power limit PowerModels' `BFAPowerModel` enforces (PowerModels writes it
      as a convex quadratic, Lirion as the equivalent cone), which makes the model an **SOCP**
      solved with `Clarabel`;
    - the shunt terms ``g^{fr}, b^{fr}, g^{to}, b^{to}`` and the complex tap are kept, so the
      losses (6), (7) and the drop (8) carry the tap and shunt data a bare radial model omits.

    The ``\pm\overline{S}_\ell`` bounds on the flow variables (2) are left in place but are
    redundant once the cone (9) is posted.

## See also

- Neighbours: [DC](@ref) · [LPAC](@ref) · [NFA](@ref)
- Parent model: SOCBF, the branch-flow second-order-cone relaxation the BFA linearizes. It has
  no page of its own; its branch relations are the ones derived above.
- W-space voltage state: [SOCWR (Jabr conic relaxation)](@ref) and [SDPWRM](@ref), which also
  carry the voltages as squared magnitudes.
- Foundations: [Foundations](@ref).
- Reference: LinDistFlow [BaranWu1989](@cite); distribution-system context [Kersting2012](@cite).
- Code: constructor `bfa()`; source `src/modeling/branch/bfa.jl`.
