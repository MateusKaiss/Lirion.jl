# NFA

The **network flow approximation** (NFA) is the most aggressive linear model in the family:
it keeps only the bookkeeping of power and discards the physics of how power actually flows.
Starting from the [DC](@ref) model, it removes the last piece of electrical structure the DC
still carried (the coupling between branch flow and nodal angles) and treats the grid as a
transport network: a set of pipes with a maximum throughput, through which the dispatcher is
free to route power as it chooses. The only quantities left are the branch flows and the
generator outputs, and the only laws left are conservation of power at the buses and the
capacity of each line.

!!! info "Relaxation or approximation?"
    The NFA carries a **double classification**, and the word changes with the pair being
    compared.
    - **Relative to the AC model it is an approximation**. Its optimum is neither a floor
      nor a ceiling of the true cost; it simply swaps the network equations for looser
      bookkeeping and hopes the error is small.
    - **Relative to the [DC](@ref) model it is a genuine relaxation**. Every DC-feasible
      dispatch is also NFA-feasible, so ``\mathcal{F}_{DC}\subseteq\mathcal{F}_{NFA}`` and
      ``\text{cost}_{NFA}\le\text{cost}_{DC}``.

    See [Relation to the DC model](@ref) below. The [DC](@ref) page's "Approximation, not
    relaxation" note defines the distinction in general.

## Derivation

### Removing Kirchhoff's voltage law

The DC model preserves more than it appears to. Its flow equation

```math
p_\ell = -b_\ell\,(\theta_i - \theta_j)
```

says that the flow on every line is **dictated by the network, not chosen**: once the bus
angles are fixed, all flows are determined. Because the angles telescope around every loop,
this is exactly Kirchhoff's voltage law (KVL) at work. Power divides between parallel paths in
inverse proportion to their reactances, whether the operator wants it or not. Electricity
cannot choose where to go.

The NFA throws this away. The flow on each branch becomes a free decision variable, bounded
only by the line's capacity:

```math
\underbrace{p_\ell = -b_\ell\,(\theta_i-\theta_j)}_{\text{DC}}
\quad\longrightarrow\quad
\underbrace{p_\ell\ \text{free},\quad |p_\ell|\le \overline{S}_\ell}_{\text{NFA}}.
```

With the angles gone, everything that depended on them leaves the model too: the angle
variables ``\theta``, the angular-aperture limits, and the **reference bus**. In the
[DC](@ref) and [LPAC](@ref) models the fixed gauge ``\theta_{ref}=0`` was mandatory to remove
the global translation of the angles; with no angles there is nothing to pin, and the
constraint simply disappears.

### What remains: two accounting laws

Once KVL is removed, only two laws survive, and both are bookkeeping rather than
electromagnetism.

**1. Conservation at the nodes (KCL).** Whatever enters a bus must leave it:

```math
\sum_{g\in\mathcal{G}(b)} p_g - P^d_b - g^{sh}_b
= \sum_{\ell\in\mathcal{L}(b)} p_\ell,
\qquad b\in\mathcal{N}^{nv}. \tag{1}
```

The term ``g^{sh}_b`` is the bus shunt, exactly as in the [DC](@ref) model. With no voltage
variable, the power it drains becomes a constant, folded in like an additional load.

**2. Capacity on the arcs.** Each line has a maximum throughput:

```math
-\overline{S}_\ell \le p_\ell \le \overline{S}_\ell,
\qquad \ell\in\mathcal{L}. \tag{2}
```

A last identity closes the accounting without adding a law. With no losses, what enters one
end of a line leaves the other, so the reverse flow is the negation of the forward one:

```math
p_{ji} = -p_{ij},
\qquad \ell\in\mathcal{L}. \tag{3}
```

And that is all. This is a **network-flow problem**, one of the best-structured objects in
combinatorial optimization. Its constraint matrix is totally unimodular (every square submatrix
has determinant ``0`` or ``\pm 1``), the property that makes this family easy. In the classical
case, with integer supplies, demands and capacities and a linear cost, the vertices of the
feasible polyhedron are integer and the linear program coincides with the underlying
combinatorial problem. Whatever the cost, dedicated network-flow algorithms solve it far faster
than a general-purpose solver, at very large scale. Stripped of all its physics, the OPF becomes
a problem of logistics.

## NFA Formulation

The NFA constructor (`nfa()`, model type `NFA`) assembles the standard skeleton (buses,
branches, generators, storage, active-power balance, objective) and gives each branch a single
free flow variable, bounded by the line capacity. The bus routine
(`add_bus_variables!(…, NFA)`) creates nothing: no angle, no voltage magnitude, no reactive
power. The flows are linear, the balance is linear, and the limits are boxes, so with a linear
fuel cost the problem is a linear program and with a quadratic fuel cost a convex quadratic
program, which Lirion solves with `OSQP` by default. The balance and capacity laws, and the
flow-reversal identity, were derived above as Eqs. (1)–(3). This section assembles the complete problem around
them, continuing the numbering from (4).

### Sets

- ``\mathcal{N}``: set of buses.
- ``\mathcal{N}^{nv} \subseteq \mathcal{N}``: non-virtual buses, on which the power balance is
  enforced (`non_virtual_bus_ids`).
- ``\mathcal{L}``: set of active branches.
- ``\mathcal{G}``: set of existing generators.
- ``\mathcal{G}^c \subseteq \mathcal{G}``: subset of curtailable generators (renewables).
- ``\mathcal{S}``: set of storage units.
- ``\mathcal{T} = \{1,\dots,T\}``: set of time periods in the horizon.

### Indices and incidence functions

- ``b \in \mathcal{N}``: bus index.
- ``\ell \in \mathcal{L}``: branch index, with *from*/*to* buses ``f(\ell),\ t(\ell)``
  (code: `br.from[l]`, `br.to[l]`).
- ``g \in \mathcal{G}``: generator index; ``s \in \mathcal{S}``: storage index.
- ``t \in \mathcal{T}``: time index.

### Parameters

Because the NFA keeps no electrical structure, the branch series/shunt/tap data play no role;
the only branch datum used is the thermal rating.

- ``\overline{S}_\ell``: thermal (apparent-power) capacity of branch ``\ell`` (code:
  `br.rate_a`).
- ``g^{sh}_b``: bus shunt conductance, folded in as a **constant** load (no voltage to
  multiply it); ``P^d_{b,t}``: active demand at bus ``b``.
- ``P_g^{\min},\ P_g^{\max}``: static active-power limits of ``g``;
  ``\overline{P}^{ts}_{g,t},\ \underline{P}^{ts}_{g,t}``: active-power time envelope.
- ``c_{g,t,k}``: cost-polynomial coefficients; ``c^{curt}_g``: unit curtailment cost.
- **Storage** (DC-style block, as in [DC](@ref)): ``E_s^{\min},E_s^{\max},E^0_s``;
  ``C_s^{\max},D_s^{\max}``; ``\eta^c_s,\eta^d_s``; ``P^{loss}_s``; ``\overline{S}_s``;
  ``Q_s^{\min},Q_s^{\max}``.

### Decision variables

```math
p^f_{\ell,t} \in [-\overline{S}_\ell,\ \overline{S}_\ell] \;\text{: active-power flow on branch }
\ell \text{ (the derivation's } p_\ell \text{; free, bounded by capacity).} \tag{4}
```

```math
p^g_{g,t} \in [P_g^{\min}, P_g^{\max}] \;\text{: active power generated by } g. \tag{5}
```

```math
p^c_{g,t} \in [0, \overline{P}^{ts}_{g,t}] \;\text{: curtailed renewable power, }
g \in \mathcal{G}^c. \tag{6}
```

There are **no** ``\theta``, ``v`` or ``Q^g`` variables. If storage units are present, the
shared block additionally creates ``E_{s,t}, c_{s,t}, d_{s,t}, P^s_{s,t}`` and, even though
the NFA has no reactive balance, ``Q^s_{s,t}, Q^{sc}_{s,t}`` (see the storage note).

### Objective function

Identical to the whole package (cf. [AC Polar](@ref)): the NFA changes only the feasible set,
not the cost.

```math
F = \sum_{t\in\mathcal{T}} \sum_{g\in\mathcal{G}} \sum_{k=1}^{m_{g,t}}
c_{g,t,k}\left(S^{base}\,p^g_{g,t}\right)^{m_{g,t} - k}. \tag{7}
```

```math
F^{curt} = \sum_{g\in\mathcal{G}^c} \sum_{t\in\mathcal{T}}
c^{curt}_g\,S^{base}\,p^c_{g,t}, \qquad \min\; F + F^{curt}. \tag{8}
```

### Branch flows

A branch carries a single free flow; the reverse flow is its negation (lossless), Eq. (3):
``p^t_{\ell,t} = -\,p^f_{\ell,t}``.

### Active power balance

Conservation at every bus, Eq. (1); with the bus shunt as a constant and, when storage is
present, the injection ``P^s_{s,t}`` joining the balance exactly as in [AC Polar](@ref),
Eq. (20):

```math
\sum_{g\in\mathcal{G}(b)} p^g_{g,t}
- \sum_{s\in\mathcal{S}(b)} P^s_{s,t}
- P^d_{b,t}
- \sum_{\ell:\,f(\ell)=b} p^f_{\ell,t}
- \sum_{\ell:\,t(\ell)=b} p^t_{\ell,t}
- g^{sh}_b = 0,
\qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{9}
```

There is **no reactive balance**. The base model has no ``q``.

### Branch capacity

The thermal limit is a **box**, Eq. (2). With no reactive power the disk
``p^2 + q^2 \le \overline{S}_\ell^2`` degenerates into the interval
``-\overline{S}_\ell \le p^f_{\ell,t} \le \overline{S}_\ell`` (imposed directly as bounds on
the flow variable).

### No voltage / angle / KVL

The NFA omits everything that depends on nodal angles or voltages:

- **no voltage magnitude** and **no ``\phi``** (the bus routine creates nothing);
- **no angle variables ``\theta``** and hence **no angular-aperture limits**;
- **no reference bus** (with no angles there is nothing to pin, so the ``\theta_{ref}=0``
  gauge of the [DC](@ref)/[LPAC](@ref) models simply disappears);
- **no Kirchhoff voltage law** (the flow ``p^f_\ell`` is a free decision, not
  ``-b_\ell(\theta_i-\theta_j)``, the relaxation of the [DC](@ref) model);
- **no cones / no reactive power**.

### Generation constraints

Time envelopes, added only when a period limit is strictly tighter than the static bound (the
upper active bound is skipped for renewables, which take the curtailment identity instead),
and the renewable curtailment identity:

```math
p^g_{g,t} \le \overline{P}^{ts}_{g,t}\ (g\in\mathcal{G}\setminus\mathcal{G}^c), \qquad
\underline{P}^{ts}_{g,t} \le p^g_{g,t}\ (g\in\mathcal{G}),
\qquad t\in\mathcal{T}. \tag{10}
```

```math
p^g_{g,t} + p^c_{g,t} = \overline{P}^{ts}_{g,t},
\qquad g\in\mathcal{G}^c,\ t\in\mathcal{T}. \tag{11}
```

### Storage constraints (DC-style)

When ``\mathcal{S}\neq\emptyset`` the shared DC storage block (`add_dc_storage_constraints!`)
is added, exactly as in the [DC](@ref) model:

```math
E_{s,1} = E^0_s + \eta^c_s\,c_{s,1} - \frac{d_{s,1}}{\eta^d_s}, \qquad
E_{s,t} = E_{s,t-1} + \eta^c_s\,c_{s,t} - \frac{d_{s,t}}{\eta^d_s}\ (t\ge 2). \tag{12}
```

```math
c_{s,t}\,d_{s,t} = 0, \qquad
(P^s_{s,t})^2 + (Q^s_{s,t})^2 \le \overline{S}_s^2, \qquad
P^s_{s,t} + d_{s,t} - c_{s,t} = P^{loss}_s. \tag{13}
```

!!! note "Reactive storage variables are network-inert"
    The storage block is shared across every model, so a run with storage still creates the
    reactive injection ``Q^s_{s,t}`` and the auxiliary ``Q^{sc}_{s,t}``, bounded by
    ``[Q_s^{\min}, Q_s^{\max}]``, and still imposes the thermal cone
    ``(P^s_{s,t})^2 + (Q^s_{s,t})^2 \le \overline{S}_s^2``. The NFA has no reactive balance, so
    ``Q^s_{s,t}`` enters none of its network equations and the DC-style block leaves
    ``Q^{sc}_{s,t}`` unused. A network with storage is therefore not strictly active-only, but
    the reactive variables stay decoupled from everything the NFA solves.

### Complete model

```math
\min_{\,p^f,\,p^g,\,p^c,\ (E,\,c,\,d,\,P^s,\,Q^s,\,Q^{sc})}
\quad F + F^{curt} \qquad (\text{Eq. } 8)
```

subject to:

| Block | Equations |
|---|---|
| Active-power balance | (1) / (9) |
| Branch capacity (box) | (2) |
| Flow reversal (lossless) | (3) |
| Generation | (10), (11) |
| Storage (DC-style) | (12), (13) |
| Bounds | (4)–(6) |

The core (no-storage) problem is this **network-flow problem**, Eqs. (1)–(8), (10)–(11), with
a totally unimodular constraint matrix.

## What it captures / what it loses

!!! info "Philosophy"
    **Captures:** nodal power balance and line capacities, the bookkeeping of a transport
    network.

    **Loses:** every piece of electrical physics. No angle limits (there are no angles), no
    angular reference (nothing to pin), no voltage limits (there is no voltage), and no cones
    (there is no reactive power).

!!! note "When to use it"
    The NFA answers questions that sit above the physics: whether the network can move a given
    transfer at all, which corridors saturate first, where a new line would relieve the most
    congestion. It is the cheapest model in the package and scales to networks where
    even the DC solve is a bottleneck, which makes it the workhorse for fast feasibility
    screening and for transmission-expansion studies that weigh many candidate topologies at
    once. And since every DC-feasible dispatch is also NFA-feasible, its cost is a valid lower
    bound on the DC cost, useful inside a larger search.

    What it is never good for is telling the generators what to do. The NFA routes power the
    way a shipping network routes freight, free to choose a split that KVL would forbid, so a
    plan it certifies can be one the real grid cannot hold (see the worked example). Use it to
    screen and to bound. For an actual dispatch, step up to the [DC](@ref) model, which at least
    honours KVL, and to [LPAC](@ref) or the [AC Polar](@ref) model for voltage and reactive
    detail.

## Relation to the DC model

This is the structural property that sets the NFA apart from the two approximations before it.

Take any DC-feasible solution. Its flows satisfy the balance (the DC imposes it) and the
capacities (the DC imposes them), so it is automatically NFA-feasible. The converse is false
whenever the network contains a loop: the NFA accepts flow patterns that no set of angles
could ever produce. Therefore

```math
\mathcal{F}_{DC}\subseteq\mathcal{F}_{NFA}
\quad\Longrightarrow\quad
\text{cost}_{NFA}\le\text{cost}_{DC}.
```

So the NFA is a **relaxation of the DC model**. The [DC](@ref) and [LPAC](@ref) models are
approximations of the AC model (they err in either direction, with no guarantee), whereas the
NFA is an approximation relative to the AC model but a relaxation relative to the DC model. The
word changes with the pair being compared.

The guarantee reaches only as far as the DC model, and no further. The containment
``\mathcal{F}_{DC}\subseteq\mathcal{F}_{NFA}`` puts ``\text{cost}_{NFA}`` below the DC cost, but
the DC cost is itself neither a floor nor a ceiling of the true AC cost. A floor under a
quantity that is already off in an unknown direction says nothing certain about the AC optimum,
which is why the same model that is a rigorous relaxation of the DC is only an approximation of
the AC.

## Worked example: two parallel paths

Two buses joined by two parallel paths, a short one (reactance ``x=0.1``) and a long one
(``x=0.4``), each rated at ``100``. We want to move ``100\ \text{MW}`` across.

In the **DC** model the split is not a decision. KVL forces the flow to divide in inverse
proportion to the reactances, so the short path takes ``80`` MW and the long one ``20``. The
operator has no say in it.

| Path | Reactance ``x`` | DC flow (forced by KVL) |
|---|---|---|
| Short | ``0.1`` | ``80`` MW |
| Long  | ``0.4`` | ``20`` MW |

In the **NFA** that split is gone. The two flows are free variables, tied only by conservation
(they must sum to ``100``) and by capacity. Any routing that fits is allowed, and because the
network is lossless every routing costs the same, so the optimizer is indifferent among them.

The gap between the two models shows the moment a line gets tight. Suppose the short path is
rated at only ``60``. KVL still demands ``80`` MW through it, so the **DC** transfer becomes
**infeasible**, even though the two paths together can carry ``160``. The **NFA** just moves
the excess onto the long path, ``60`` MW short and ``40`` MW long, and the ``100`` MW goes
through. Any split with the short path at or below ``60`` works, ``60/40`` being the one that
loads the capped line to its limit.

!!! warning "Not an operating plan"
    The NFA has declared feasible a transfer the real network cannot hold. No setting of the
    angles produces a ``60/40`` split on these two lines; KVL allows only ``80/20``. This is
    the optimistic error a relaxation always risks, and the reason an NFA solution is a
    screening or bounding result and never a dispatch to act on.

## See also

- Neighbours: [DC](@ref) · [LPAC](@ref)
- Foundations: [Foundations](@ref) (notation)
- Reference model: [AC Polar](@ref) (the exact flow equations the whole family approximates)
- Code: constructor `nfa()`; source `src/modeling/branch/nfa.jl`.
