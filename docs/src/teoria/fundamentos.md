# Foundations

The theory pages each stand on their own, but they all share one network model and one
set of conventions. This page collects them, so the formulation pages do not repeat the
groundwork.

## Sets and indices

| Symbol | Meaning |
|---|---|
| ``\mathcal{N}`` | buses; ``\mathcal{N}^{nv}\subseteq\mathcal{N}`` the non-virtual buses where the power balance is enforced |
| ``\mathcal{L}`` | active branches; ``f(\ell),\ t(\ell)`` the *from* and *to* buses of branch ``\ell`` |
| ``\mathcal{G}`` | generators; ``\mathcal{G}(b)`` those at bus ``b``; ``\mathcal{G}^c`` the curtailable (renewable) ones |
| ``\mathcal{S}`` | storage units; ``\mathcal{S}(b)`` those at bus ``b`` |
| ``\mathcal{T}=\{1,\dots,T\}`` | time periods |

A **virtual bus** carries no load, generation, or shunt. It exists only to glue zones
together in a distributed solve, so the nodal balance skips it, which is what
``\mathcal{N}^{nv}`` selects.

## The branch model

Every branch is a ``\pi``-section behind an ideal complex transformer. From the
resistance ``r_\ell`` and reactance ``x_\ell`` the series admittance is

```math
y_\ell = g_\ell + j\,b_\ell,\qquad
g_\ell = \frac{r_\ell}{r_\ell^2 + x_\ell^2},\qquad
b_\ell = \frac{-x_\ell}{r_\ell^2 + x_\ell^2}.
```

The line charging sits at the two ends as shunt admittances ``y^{fr}_\ell =
g^{fr}_\ell + j\,b^{fr}_\ell`` and ``y^{to}_\ell = g^{to}_\ell + j\,b^{to}_\ell``. A
plain line splits its total charging in half, ``b^{fr}_\ell = b^{to}_\ell =
b^{sh}_\ell/2`` with ``g^{fr}_\ell = g^{to}_\ell = 0``. The transformer is a complex
tap on the *from* side,

```math
T_\ell = \tau_\ell\, e^{j\theta^{shift}_\ell},\qquad
tm_\ell = |T_\ell|^2 = \tau_\ell^2,
```

with tap ratio ``\tau_\ell`` and phase shift ``\theta^{shift}_\ell``. A branch with no
transformer has ``\tau_\ell = 1`` and ``\theta^{shift}_\ell = 0``.

## The branch admittance block

Put together, a branch ties its two terminal currents to its two terminal voltages
through a ``2\times 2`` block,

```math
\begin{bmatrix} I_f \\ I_t \end{bmatrix}
=
\begin{bmatrix} Y_{ff} & Y_{ft} \\ Y_{tf} & Y_{tt} \end{bmatrix}
\begin{bmatrix} V_f \\ V_t \end{bmatrix},
\qquad
\begin{aligned}
Y_{ff} &= \frac{y_\ell + y^{fr}_\ell}{tm_\ell}, &\quad
Y_{ft} &= -\frac{y_\ell}{\overline{T_\ell}}, \\[2pt]
Y_{tf} &= -\frac{y_\ell}{T_\ell}, &\quad
Y_{tt} &= y_\ell + y^{to}_\ell.
\end{aligned}
```

This block is the single object every formulation starts from. The AC pages expand it
into real and reactive flows, the W-space relaxations lift the voltage products
``V_i\,\overline{V_j}``, and the branch-flow models rewrite it around the squared
current magnitude. On a plain line (``T_\ell = 1``) it collapses to the textbook
``\pi``-model, ``Y_{ff} = y_\ell + y^{fr}_\ell``, ``Y_{ft} = Y_{tf} = -y_\ell``,
``Y_{tt} = y_\ell + y^{to}_\ell``.

## Conventions

- **Orientation.** Flows are measured at each end: ``P^f_\ell, Q^f_\ell`` at the *from*
  bus and ``P^t_\ell, Q^t_\ell`` at the *to* bus. Losses make ``P^f_\ell \neq
  -P^t_\ell`` in general.
- **Nodal balance.** At every non-virtual bus, generation minus storage minus load
  minus the flows leaving the bus minus the bus shunt ``y^{sh}_b|V_b|^2`` is zero.
- **Voltage state.** Polar models carry ``V_b`` and ``\theta_b``; rectangular models
  carry ``v^{re}_b, v^{im}_b``; the W-space relaxations carry the squared magnitude
  ``w_b = |V_b|^2`` and the bus-pair products. One bus is pinned as the angle
  reference, ``\theta_{ref} = 0``.
- **Per unit.** Everything is in per unit on the system base, so impedances, powers,
  and voltages are dimensionless and the relations above carry no scaling constants.
