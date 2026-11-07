# DC

All models up to this point carry the full structure of the AC power flow: voltage
magnitudes, angles, reactive power, cones, envelopes. The DC model goes in the opposite
direction and asks: *how much physics can be thrown away and still leave a useful model?* By
discarding almost everything, what remains is a **linear** problem, the cheapest of all to
solve, and one that still reproduces the active-power flow of transmission networks
reasonably well. It is formulated in the bus angles ``\theta_i`` and the active generation
``p^g_g``.

!!! info "Approximation, not relaxation"
    The DC model is an **approximation**, not a relaxation. Relaxations
    ([SOCWR (Jabr conic relaxation)](@ref), [QC (Quadratic Convex)](@ref), [SDPWRM](@ref))
    *enlarge* the feasible region, so their optimum is a guaranteed **lower bound** (a floor)
    on the true cost. An approximation instead *replaces* the equations by similar ones: the
    new feasible region neither contains nor is contained in the true one, and the optimum is
    **neither a floor nor a ceiling**. There is no guarantee at all, only the expectation that
    the error is small.

## Derivation

The starting point is the exact active-power flow, derived on the [AC Polar](@ref) page (the
*from*-side active flow, written there with the inverse tap ``k_\ell = 1/\tau``):

```math
p_{ij} = \frac{g + g^{fr}}{tm}\,v_i^2
       - \frac{1}{\tau}\,v_i v_j\,(g\cos\delta + b\sin\delta).
```

Here ``tm = |T|^2 = \tau^2`` is the squared tap magnitude, ``\tau`` the tap ratio, and
``\delta = \theta_i - \theta_j - \theta^{shift}`` the effective angular aperture, all as on
the AC Polar page.

We apply four approximations, one at a time, and watch each piece fall away.

### Approximation 0: no transformer, no line charging

The DC model ignores taps and phase shifters (``\tau = 1``, ``\theta^{shift} = 0``), so
``tm = \tau^2 = 1`` and ``\delta = \theta_i - \theta_j = \theta_{ij}``. It also discards the
line charging (``g^{fr} = 0``). What remains is

```math
p_{ij} = g\,v_i^2 - v_i v_j\,(g\cos\theta_{ij} + b\sin\theta_{ij}).
```

!!! note "Justification"
    Most transmission lines carry no tap and no phase shifter, so on them ``\tau = 1`` and
    ``\theta^{shift} = 0`` hold exactly rather than as an approximation. The line charging
    enters the active flow only through its shunt conductance ``g^{fr}``, which is essentially
    zero since charging is a reactive effect, so dropping it costs almost nothing. The one
    part that can bite is a real tap or phase shifter where the network has one, and that is
    the first thing [LPAC](@ref) gives back.

### Approximation 1: all voltages at 1 pu

Set ``v_i = v_j = 1``:

```math
p_{ij} = g - (g\cos\theta_{ij} + b\sin\theta_{ij}).
```

!!! note "Physical justification"
    In well-operated transmission networks the magnitudes stay in a narrow band around
    nominal, typically between ``0.95`` and ``1.05`` pu. This is no accident: keeping the
    voltage near nominal is the explicit operational objective, sustained by regulators,
    capacitors and reactive dispatch. Treating ``v \approx 1`` therefore errs little.

### Approximation 2: small angles, ``\cos\theta_{ij} \approx 1``

Substituting ``\cos\theta_{ij} \approx 1``:

```math
p_{ij} = g - g\cdot 1 - b\sin\theta_{ij} = (g - g) - b\sin\theta_{ij} = -\,b\sin\theta_{ij}.
```

!!! note "Justification"
    Angle differences between neighbouring buses are small in normal operation, rarely beyond
    ``10^\circ`` to ``30^\circ``. And the cosine is flat near zero: its slope there is
    ``-\sin 0 = 0``, so near the peak the value barely moves. At ``10^\circ``,
    ``\cos = 0.985``; at ``20^\circ``, ``0.940``. Replacing it by ``1`` costs ``1.5\%`` in one
    case and ``6\%`` in the other.

The conductance cancelled **by itself**. It was not necessary to assume ``g \approx 0``: under
``v = 1`` and ``\cos \approx 1`` the self term (``g\,v_i^2``) and the real part of the cross
term (``v_i v_j\,g\cos``) become identical and annihilate. So the DC model is **born
lossless**: what leaves one bus is exactly what arrives at ``j``. This is not a modelling
choice, it is the raw result. The same four approximations on the *to*-side flow give
``p_{ji} = +\,b\sin\theta_{ij}``, the exact negative of ``p_{ij}``, so nothing is lost
crossing the line and the reverse flow is the identity ``P^t_\ell = -P^f_\ell`` rather than a
new equation.

### Approximation 3: ``\sin\theta_{ij} \approx \theta_{ij}``

```math
p_{ij} = -\,b\,(\theta_i - \theta_j).
```

This is the whole model: a single expression, linear in the angles.

!!! note "Justification"
    This is the geometric fact that ``\sin x \approx x`` near zero, the two curves hugging
    because the sine has an inflection at zero (zero curvature, so it is almost straight
    there). At ``10^\circ``, ``\sin = 0.1736`` against ``\theta = 0.1745`` rad, a ``0.5\%``
    error. At ``30^\circ``, ``0.5000`` against ``0.5236``, a ``4.5\%`` error.

!!! note "Reading the formula"
    The susceptance is ``b = -x/(r^2 + x^2)``, negative for the inductive lines a grid is
    built from, so ``-b = x/(r^2 + x^2) > 0``. Positive power flows from the higher angle to
    the lower one, downhill in angle. When the resistance ``r`` is small against the reactance
    ``x``, ``-b \approx 1/x`` and the flow is close to the angle difference divided by the
    line reactance. That is the mental picture most people keep of the DC model, and it comes
    straight out of ``p_{ij} = -b\,(\theta_i - \theta_j)``.

### Approximation 4: no reactive power

Approximation 1 froze the voltage magnitudes, and reactive power is what rides on them. Once
they are fixed, ``q`` has nothing left to move against, so the model drops it whole. The
justification is the classical decoupling. Near the nominal operating point active power is
governed chiefly by the angles and reactive power by the magnitudes, and since the magnitudes
were frozen at ``1.0``, reactive power loses its variables and leaves the model. Everything
that depended on it goes with it: the voltage limits, capacitors, reactive support, and the
reactive part of the thermal limit.

---

## DC Formulation

The model keeps only bus angles and active generation. There are no magnitudes, no reactive
power, and no ``W`` blocks. The flows are linear, the balance is linear, and the limits are
boxes, so the feasible set is a polyhedron. With a linear fuel cost the problem is a linear
program, and with a quadratic fuel cost a convex quadratic program. Either way it is the
easiest model in the package to solve, which is the whole reason the DC model exists. Lirion
solves it with `OSQP` by default.

### Sets

- ``\mathcal{N}``: buses; ``\mathcal{N}^{nv} \subseteq \mathcal{N}``: non-virtual buses, on
  which the power balance is enforced (`non_virtual_bus_ids`).
- ``\mathcal{L}``: active branches; ``f(\ell), t(\ell)``: the *from* and *to* buses of ``\ell``.
- ``\mathcal{G}``: generators; ``\mathcal{G}(b)``: those at bus ``b``.

The parameters ``b_\ell`` (series susceptance), ``P^d_b`` (active demand), ``g^{sh}_b`` (bus
shunt conductance), ``\overline{S}_\ell`` (apparent-power limit, `rate_a`) and the angular
limits ``\theta^{\min}_\ell, \theta^{\max}_\ell`` are as on the [AC Polar](@ref) page.

### Decision variables

```math
\theta_i \in \mathbb{R} \;\text{: voltage angle at bus } i. \tag{1}
```

```math
p^g_g \in [P_g^{\min}, P_g^{\max}] \;\text{: active power generated by } g. \tag{2}
```

### Branch flows

```math
P^f_\ell = -\,b_\ell\,(\theta_{f(\ell)} - \theta_{t(\ell)}),
\qquad P^t_\ell = -\,P^f_\ell,
\qquad \ell\in\mathcal{L}. \tag{3}
```

The reverse flow is the **exact negative** of the forward one. Lossless, as Approximation 2
derived.

### Power balance

One linear equality per bus. Everything that comes in must go out:

```math
\sum_{g\in\mathcal{G}(b)} p^g_g - P^d_b - g^{sh}_b
= \sum_{\ell:\,f(\ell)=b} P^f_\ell + \sum_{\ell:\,t(\ell)=b} P^t_\ell,
\qquad b\in\mathcal{N}^{nv}. \tag{4}
```

The constant ``g^{sh}_b`` is the bus shunt conductance. It drew ``g^{sh}_b\,v_b^2`` in the AC
model; with ``v_b = 1`` that collapses to a constant, folded in beside the load, the same way
the [NFA](@ref) model carries it. The balance is imposed on the non-virtual buses.

### Reference bus

```math
\theta_{ref} = 0. \tag{5}
```

The flows depend only on angle differences, so the whole set of angles can shift by a common
amount with nothing changing. Pinning one bus to ``\theta_{ref} = 0`` removes that free offset
and makes the solution unique.

### Branch thermal limits

Without reactive power the disc ``p^2 + q^2 \le \overline{S}_\ell^2`` degenerates into an
interval, a box in the active flow:

```math
-\,\overline{S}_\ell \le P^f_\ell \le \overline{S}_\ell,
\qquad \ell\in\mathcal{L}. \tag{6}
```

### Angular limits

```math
\theta^{\min}_\ell \le \theta_{f(\ell)} - \theta_{t(\ell)} \le \theta^{\max}_\ell,
\qquad \ell\in\mathcal{L}. \tag{7}
```

### Complete model

```math
\min_{\theta,\,p^g}\quad F \qquad (\text{generation cost})
```

subject to:

| Block | Equations |
|---|---|
| Flows | (3) |
| Balance | (4) |
| Reference | (5) |
| Thermal | (6) |
| Angle | (7) |
| Bounds | (1)–(2) |

The objective ``F`` is the generation cost every Lirion formulation minimizes: fuel cost on
the active generation ``p^g`` in physical units, plus renewable curtailment cost where the
system has curtailable generators. The DC builder inherits the same renewable-curtailment and
DC-style storage blocks as the richer models. Storage adds constraints but no cost term, so
the objective stays fuel plus curtailment. The sections above keep the core DC physics in view
and leave those shared blocks implicit.

## What it captures / what it loses

!!! info "Philosophy"
    **Captures:** the active-power flow pattern on a transmission network as a function of the
    bus angles (``P^f_\ell = -b_\ell\,\Delta\theta``), the nodal active balance, thermal box
    limits and angular apertures. **Loses:** voltage magnitudes (frozen at ``1`` pu), reactive
    power entirely, all losses (``P^t_\ell = -P^f_\ell`` exactly), transformer taps and phase
    shifts, the branch charging shunt ``g^{fr}``, and the reactive part of the thermal limit.
    The bus shunt conductance ``g^{sh}_b`` is kept, as a constant load. Being an approximation,
    its optimum is neither a floor nor a ceiling on the true AC cost.

!!! note "When to trust it"
    The DC model is at its best exactly where its four approximations were argued: a
    transmission network in normal operation, voltages near nominal, angle differences small.
    There it is fast and close enough to run thousands of times in planning and markets. It
    stops being reliable when the things it dropped start to matter, under heavy loading with
    wide angles, when losses or voltage support are the question, or on distribution feeders
    where resistance is not small against reactance. For those cases [LPAC](@ref) gives the
    physics back a layer at a time, and the AC models give all of it.

## See also

- Neighbours: [LPAC](@ref) (restores voltage, reactive power and losses) ·
  [NFA](@ref) (a true relaxation *of the DC model*).
- Exact starting point: [AC Polar](@ref) (the flow where the four approximations begin).
- Foundations: [Foundations](@ref).
