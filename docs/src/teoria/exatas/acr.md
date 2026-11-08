# AC Rectangular (ACR)

The AC rectangular model (ACR) is **not a new model**. It is the [AC Polar](@ref)
model seen through a different coordinate system. The physics stays the same, from the
complex power injection ``S = V\,I^{*}`` to the branch admittance matrix; only the way we
write the voltage changes. In the polar model the nodal voltage is a
magnitude and an angle, ``V_i = v_i\,e^{j\theta_i}``; here it is its Cartesian
components,

```math
V_i = v^{\mathrm{re}}_i + j\,v^{\mathrm{im}}_i .
```

It is the **same complex point**, rewritten. Consequently the ACR is still the exact,
nonlinear (nonconvex) representation of power flow, and it minimizes the same generation
cost subject to the same network equations. The only change is that the decision variables are now the real
and imaginary parts of the voltage, ``v^{\mathrm{re}}_{b,t}`` and ``v^{\mathrm{im}}_{b,t}``,
and every trigonometric term is replaced by a bilinear expression in those two variables.

## Derivation

### Rectangular voltage

Euler's identity relates the two coordinate systems directly:

```math
V_i = v_i\,e^{j\theta_i} = v_i\,(\cos\theta_i + j\,\sin\theta_i) = v^{\mathrm{re}}_i + j\,v^{\mathrm{im}}_i ,
```

so that the Cartesian components are the projections of the phasor onto the real and
imaginary axes,

```math
v^{\mathrm{re}}_i = v_i\cos\theta_i, \qquad v^{\mathrm{im}}_i = v_i\sin\theta_i .
```

This is the same complex point already drawn for the [AC Polar](@ref) model. The
rectangular decomposition is exactly the projection of the voltage phasor onto the real
and imaginary axes shown there. From now on the decision variables are
``v^{\mathrm{re}}_i`` and ``v^{\mathrm{im}}_i``.

### The four fundamental combinations

All branch flows of the polar model can be written in terms of four scalar combinations of
the end-bus voltages, the **W blocks**, each of which equals a polar quantity.

**Squared magnitudes (combinations 1 and 2).** By Pythagoras, ``v_i = \sqrt{(v^{\mathrm{re}}_i)^2 + (v^{\mathrm{im}}_i)^2}``,
hence

```math
W_{ii} = v_i^2 = (v^{\mathrm{re}}_i)^2 + (v^{\mathrm{im}}_i)^2, \qquad
W_{jj} = v_j^2 = (v^{\mathrm{re}}_j)^2 + (v^{\mathrm{im}}_j)^2 .
```

**Cross terms (combinations 3 and 4).** In the polar model the product of the two end
voltages is ``V_i\,V_j^{*} = v_i v_j\,e^{j(\theta_i-\theta_j)} = v_i v_j(\cos\theta_{ij} + j\sin\theta_{ij})``,
with ``\theta_{ij} = \theta_i - \theta_j``. Writing the same product in rectangular
coordinates,

```math
V_i\,V_j^{*} = (v^{\mathrm{re}}_i + j\,v^{\mathrm{im}}_i)(v^{\mathrm{re}}_j - j\,v^{\mathrm{im}}_j)
= \underbrace{(v^{\mathrm{re}}_i v^{\mathrm{re}}_j + v^{\mathrm{im}}_i v^{\mathrm{im}}_j)}_{\mathrm{Re}(V_i V_j^{*})}
+ j\,\underbrace{(v^{\mathrm{im}}_i v^{\mathrm{re}}_j - v^{\mathrm{re}}_i v^{\mathrm{im}}_j)}_{\mathrm{Im}(V_i V_j^{*})} .
```

Matching real and imaginary parts term by term gives the remaining two blocks:

```math
W_r = v_i v_j\cos\theta_{ij} = v^{\mathrm{re}}_i v^{\mathrm{re}}_j + v^{\mathrm{im}}_i v^{\mathrm{im}}_j, \qquad
W_i = v_i v_j\sin\theta_{ij} = v^{\mathrm{im}}_i v^{\mathrm{re}}_j - v^{\mathrm{re}}_i v^{\mathrm{im}}_j .
```

The four blocks ``W_{ii},W_{jj},W_r,W_i`` are the heart of the rectangular rewrite. Each
is an exact algebraic image of a polar quantity, but written without any trigonometry.

### Rewriting the flows in W

Start from the *from*-side flows of the polar model,

```math
p_{ij} = \frac{g + g^{fr}}{tm}\,v_i^2 - \frac{1}{\tau}\,v_i v_j\,(g\cos\delta + b\sin\delta),
```
```math
q_{ij} = -\frac{b + b^{fr}}{tm}\,v_i^2 - \frac{1}{\tau}\,v_i v_j\,(-b\cos\delta + g\sin\delta),
```

where ``\delta = \theta_{ij} - \theta_{shift}`` is the effective angular aperture. The goal
is to express everything in ``W_{ii},W_{jj},W_r,W_i`` alone. Expanding ``\delta`` with the
angle-difference identities,

```math
\cos\delta = \cos\theta_{ij}\cos\theta_{shift} + \sin\theta_{ij}\sin\theta_{shift}, \qquad
\sin\delta = \sin\theta_{ij}\cos\theta_{shift} - \cos\theta_{ij}\sin\theta_{shift},
```

the active-power trigonometric factor becomes

```math
g\cos\delta + b\sin\delta
= (g\cos\theta_{shift} - b\sin\theta_{shift})\cos\theta_{ij}
+ (g\sin\theta_{shift} + b\cos\theta_{shift})\sin\theta_{ij} .
```

Now trade the phase shift for the **components of the complex tap**
``T = \tau\,e^{j\theta_{shift}} = \tau(\cos\theta_{shift} + j\sin\theta_{shift})``, namely

```math
t_r = \mathrm{Re}(T) = \tau\cos\theta_{shift}, \qquad t_i = \mathrm{Im}(T) = \tau\sin\theta_{shift}.
```

Substituting ``\cos\theta_{shift} = t_r/\tau`` and ``\sin\theta_{shift} = t_i/\tau`` and
recalling ``tm = |T|^2 = \tau^2``:

```math
p_{ij} = \frac{g + g^{fr}}{tm}\,\underbrace{v_i^2}_{W_{ii}}
- \frac{1}{\tau^2}\Big[(g\,t_r - b\,t_i)\,\underbrace{v_i v_j\cos\theta_{ij}}_{W_r}
+ (g\,t_i + b\,t_r)\,\underbrace{v_i v_j\sin\theta_{ij}}_{W_i}\Big] ,
```

which collects into a purely linear expression in the blocks:

```math
p_{ij} = \frac{g + g^{fr}}{tm}\,W_{ii}
+ \frac{-g\,t_r + b\,t_i}{tm}\,W_r
+ \frac{-b\,t_r - g\,t_i}{tm}\,W_i .
```

The same expansion applied to the reactive *from* flow and to both *to*-side flows yields:

```math
q_{ij} = -\frac{b + b^{fr}}{tm}\,W_{ii}
- \frac{-b\,t_r - g\,t_i}{tm}\,W_r
+ \frac{-g\,t_r + b\,t_i}{tm}\,W_i ,
```
```math
p_{ji} = (g + g^{to})\,W_{jj}
+ \frac{-g\,t_r - b\,t_i}{tm}\,W_r
+ \frac{-g\,t_i + b\,t_r}{tm}\,W_i ,
```
```math
q_{ji} = -(b + b^{to})\,W_{jj}
- \frac{g\,t_i - b\,t_r}{tm}\,W_r
- \frac{-g\,t_r - b\,t_i}{tm}\,W_i .
```

Every flow is now **linear in the W blocks**. This is the key structural fact that the convex
relaxations build on later. In the exact model the blocks are still quadratic functions of
``v^{\mathrm{re}},v^{\mathrm{im}}``, so the model remains nonconvex; what changed is that the
nonconvexity is now confined to the four block definitions.

### Main constraints

**Thermal limit.** Identical in form to the polar model. The apparent-flow magnitude is
squared to keep the constraint smooth:

```math
p_{ij}^2 + q_{ij}^2 \le rate_a^2, \qquad p_{ji}^2 + q_{ji}^2 \le rate_a^2 .
```

**Angle limit.** The rectangular model has no angle variable, so the angular aperture is
imposed through its tangent. Since

```math
\tan\theta_{ij} = \frac{\sin\theta_{ij}}{\cos\theta_{ij}}
= \frac{v_i v_j\sin\theta_{ij}}{v_i v_j\cos\theta_{ij}} = \frac{W_i}{W_r},
```

the bounds ``\theta^{\min} \le \theta_i - \theta_j \le \theta^{\max}`` translate into two
linear inequalities in the blocks:

```math
\theta_i - \theta_j \le \theta^{\max} \;\Longrightarrow\; W_i \le W_r\tan\theta^{\max},
```
```math
\theta^{\min} \le \theta_i - \theta_j \;\Longrightarrow\; W_r\tan\theta^{\min} \le W_i .
```

---

## AC Rectangular Formulation

The formal model reuses the polar sets, parameters, generation and storage blocks
verbatim; only the two voltage variables, the branch flows, the shunt term of the balance,
the voltage-magnitude bounds, the angle limit and the reference-bus condition are rewritten
in rectangular coordinates. In the formal notation the derivation's short symbols map one
to one: ``k_\ell = 1/\mathrm{tap}`` (so ``k_\ell^2 = 1/tm``), ``T^{\mathrm{re}}_\ell = t_r``,
``T^{\mathrm{im}}_\ell = t_i``, and the blocks become the per-branch quantities
``|V_f|^2 = W_{ii}``, ``|V_t|^2 = W_{jj}``, ``C_{\ell,t} = W_r``, ``S_{\ell,t} = W_i``.

### Sets

- ``\mathcal{N}``: set of buses.
- ``\mathcal{N}^{nv} \subseteq \mathcal{N}``: non-virtual buses, on which the nodal power balance is enforced.
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
- ``f(\ell),\ t(\ell)``: *from* and *to* buses of branch ``\ell``.
- ``B(s)``: connection bus of storage ``s``.

### Parameters

**Branches.**

- ``g_\ell,\ b_\ell``: series conductance and susceptance of branch ``\ell``.
- ``g^{fr}_\ell,\ b^{fr}_\ell``: shunt conductance and susceptance, *from* side.
- ``g^{to}_\ell,\ b^{to}_\ell``: shunt conductance and susceptance, *to* side.
- ``k_\ell``: inverse tap ratio (``1/\mathrm{tap}``).
- ``T^{\mathrm{re}}_\ell,\ T^{\mathrm{im}}_\ell``: components of the complex tap,
  ``\mathrm{tap}\cos\Phi_\ell`` and ``\mathrm{tap}\sin\Phi_\ell``, with ``\Phi_\ell`` the transformer phase shift.
- ``\overline{S}_\ell``: apparent-flow limit of the branch.
- ``\delta_\ell^{\min},\ \delta_\ell^{\max}``: angular aperture limits of the branch.

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

**Storage.**

- ``E_s^{\min},\ E_s^{\max}``: stored-energy limits.
- ``E^0_s``: initial energy of unit ``s``.
- ``C_s^{\max},\ D_s^{\max}``: maximum charge and discharge power.
- ``\eta^c_s,\ \eta^d_s``: charge and discharge efficiencies.
- ``r_s,\ x_s``: converter resistance and reactance.
- ``P^{loss}_s,\ Q^{loss}_s``: fixed active and reactive losses.
- ``\overline{S}_s``: inverter thermal limit.
- ``Q_s^{\min},\ Q_s^{\max}``: reactive-power limits of the storage.

### Decision variables

The only difference from the polar formulation is in the first two variables. Instead of
``(\theta_{b,t}, V_{b,t})``, the voltage is described by its real part ``v^{\mathrm{re}}_{b,t}``
and imaginary part ``v^{\mathrm{im}}_{b,t}``. The remaining generation and storage variables
are unchanged.

**Voltage.**

```math
v^{\mathrm{re}}_{b,t} \in [-V_b^{\max}, V_b^{\max}] \;\text{: real part of the voltage at bus } b \text{ at time } t. \tag{1}
```

```math
v^{\mathrm{im}}_{b,t} \in [-V_b^{\max}, V_b^{\max}] \;\text{: imaginary part of the voltage at bus } b. \tag{2}
```

**Generation.**

```math
p^g_{g,t} \in [P_g^{\min}, P_g^{\max}] \;\text{: active power generated by } g. \tag{3}
```

```math
Q^g_{g,t} \in [Q_g^{\min}, Q_g^{\max}] \;\text{: reactive power generated by } g. \tag{4}
```

```math
p^c_{g,t} \in [0, \overline{P}^{ts}_{g,t}] \;\text{: curtailed renewable power, } g \in \mathcal{G}^c. \tag{5}
```

**Storage.**

```math
E_{s,t} \in [E_s^{\min}, E_s^{\max}] \;\text{: energy stored in } s \text{ at the end of } t. \tag{6}
```

```math
c_{s,t} \in [0, C_s^{\max}] \;\text{: charging power of } s. \tag{7}
```

```math
d_{s,t} \in [0, D_s^{\max}] \;\text{: discharging power of } s. \tag{8}
```

```math
P^s_{s,t} \in [-\overline{S}_s, \overline{S}_s] \;\text{: net active power exchanged by } s \text{, positive when charging (drawn from the grid).} \tag{9}
```

```math
Q^s_{s,t} \in [Q_s^{\min}, Q_s^{\max}] \;\text{: net reactive power exchanged by } s \text{, positive when absorbed from the grid.} \tag{10}
```

```math
Q^{sc}_{s,t} \in [Q_s^{\min}, Q_s^{\max}] \;\text{: internal reactive compensation of the converter of } s. \tag{11}
```

### Objective function

Identical to the polar model. The fuel-cost polynomial is evaluated on the power in physical
units, plus the curtailment cost:

```math
F = \sum_{t\in\mathcal{T}} \sum_{g\in\mathcal{G}} \sum_{k=1}^{m_{g,t}}
c_{g,t,k}\left(S^{base}\,p^g_{g,t}\right)^{m_{g,t} - k}, \tag{12}
```

```math
F^{curt} = \sum_{g\in\mathcal{G}^c} \sum_{t\in\mathcal{T}} c^{curt}_g\,S^{base}\,p^c_{g,t}. \tag{13}
```

```math
\min\; F + F^{curt}. \tag{14}
```

### Branch flows

Here lies the central difference from the polar formulation. Instead of ``\delta_{\ell,t}``
and trigonometric functions, four auxiliary expressions bilinear in
``v^{\mathrm{re}},v^{\mathrm{im}}`` are defined. Writing
``v^{\mathrm{re}}_f \equiv v^{\mathrm{re}}_{f(\ell),t}``,
``v^{\mathrm{re}}_t \equiv v^{\mathrm{re}}_{t(\ell),t}`` (and likewise for
``v^{\mathrm{im}}``):

**Squared magnitudes.**

```math
|V_f|^2_{\ell,t} = (v^{\mathrm{re}}_f)^2 + (v^{\mathrm{im}}_f)^2, \qquad
|V_t|^2_{\ell,t} = (v^{\mathrm{re}}_t)^2 + (v^{\mathrm{im}}_t)^2 . \tag{15}
```

**Cross terms.**

```math
C_{\ell,t} = v^{\mathrm{re}}_f v^{\mathrm{re}}_t + v^{\mathrm{im}}_f v^{\mathrm{im}}_t, \qquad
S_{\ell,t} = v^{\mathrm{im}}_f v^{\mathrm{re}}_t - v^{\mathrm{re}}_f v^{\mathrm{im}}_t . \tag{16}
```

The identification with the polar form is exact. If ``V_f = v_f e^{j\theta_f}`` and
``V_t = v_t e^{j\theta_t}``, then ``C_{\ell,t} = v_f v_t\cos(\theta_f - \theta_t)`` and
``S_{\ell,t} = v_f v_t\sin(\theta_f - \theta_t)``. It is the same model, written without
trigonometry. With ``k_\ell = 1/\mathrm{tap}`` and the complex tap
``(T^{\mathrm{re}}_\ell, T^{\mathrm{im}}_\ell)``, the four flows are:

```math
P^f_{\ell,t} = (g_\ell + g^{fr}_\ell)\,k_\ell^2\,|V_f|^2_{\ell,t}
+ (-g_\ell T^{\mathrm{re}}_\ell + b_\ell T^{\mathrm{im}}_\ell)\,k_\ell^2\,C_{\ell,t}
+ (-b_\ell T^{\mathrm{re}}_\ell - g_\ell T^{\mathrm{im}}_\ell)\,k_\ell^2\,S_{\ell,t}, \tag{17}
```
```math
Q^f_{\ell,t} = -(b_\ell + b^{fr}_\ell)\,k_\ell^2\,|V_f|^2_{\ell,t}
- (-b_\ell T^{\mathrm{re}}_\ell - g_\ell T^{\mathrm{im}}_\ell)\,k_\ell^2\,C_{\ell,t}
+ (-g_\ell T^{\mathrm{re}}_\ell + b_\ell T^{\mathrm{im}}_\ell)\,k_\ell^2\,S_{\ell,t}, \tag{18}
```
```math
P^t_{\ell,t} = (g_\ell + g^{to}_\ell)\,|V_t|^2_{\ell,t}
+ (-g_\ell T^{\mathrm{re}}_\ell - b_\ell T^{\mathrm{im}}_\ell)\,k_\ell^2\,C_{\ell,t}
- (-b_\ell T^{\mathrm{re}}_\ell + g_\ell T^{\mathrm{im}}_\ell)\,k_\ell^2\,S_{\ell,t}, \tag{19}
```
```math
Q^t_{\ell,t} = -(b_\ell + b^{to}_\ell)\,|V_t|^2_{\ell,t}
- (-b_\ell T^{\mathrm{re}}_\ell + g_\ell T^{\mathrm{im}}_\ell)\,k_\ell^2\,C_{\ell,t}
- (-g_\ell T^{\mathrm{re}}_\ell - b_\ell T^{\mathrm{im}}_\ell)\,k_\ell^2\,S_{\ell,t} . \tag{20}
```

### Power balance

Structurally identical to the polar model; only the shunt term changes, now using
``(v^{\mathrm{re}}_{b,t})^2 + (v^{\mathrm{im}}_{b,t})^2`` in place of ``V_{b,t}^2``.

**Active power.**

```math
\sum_{g\in\mathcal{G}(b)} p^g_{g,t}
- \sum_{s\in\mathcal{S}(b)} P^s_{s,t}
- P^d_{b,t}
- \sum_{\ell:\,f(\ell)=b} P^f_{\ell,t}
- \sum_{\ell:\,t(\ell)=b} P^t_{\ell,t}
- g^{sh}_b\big((v^{\mathrm{re}}_{b,t})^2 + (v^{\mathrm{im}}_{b,t})^2\big) = 0,
\qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{21}
```

**Reactive power.**

```math
\sum_{g\in\mathcal{G}(b)} Q^g_{g,t}
- \sum_{s\in\mathcal{S}(b)} Q^s_{s,t}
- Q^d_{b,t}
- \sum_{\ell:\,f(\ell)=b} Q^f_{\ell,t}
- \sum_{\ell:\,t(\ell)=b} Q^t_{\ell,t}
+ b^{sh}_b\big((v^{\mathrm{re}}_{b,t})^2 + (v^{\mathrm{im}}_{b,t})^2\big) = 0,
\qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{22}
```

### Network constraints

**Voltage magnitude limits.** The magnitude bound is imposed on the sum of squares:

```math
(V_b^{\min})^2 \le (v^{\mathrm{re}}_{b,t})^2 + (v^{\mathrm{im}}_{b,t})^2 \le (V_b^{\max})^2,
\qquad b\in\mathcal{N},\ t\in\mathcal{T}. \tag{23}
```

Geometrically, the feasible voltage of each bus lives in a ring (annulus) of inner radius
``V_b^{\min}`` and outer radius ``V_b^{\max}``.

**Branch thermal limits.**

```math
(P^f_{\ell,t})^2 + (Q^f_{\ell,t})^2 \le \overline{S}_\ell^2, \qquad
(P^t_{\ell,t})^2 + (Q^t_{\ell,t})^2 \le \overline{S}_\ell^2,
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{24}
```

**Angular aperture limits.** Imposed through the tangent, as a pair of linear inequalities
in the cross terms:

```math
\tan(\delta_\ell^{\min})\,C_{\ell,t} \le S_{\ell,t} \le \tan(\delta_\ell^{\max})\,C_{\ell,t},
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{25}
```

These rows are added only for branches that declare an aperture limit; where a side is left
open the inequality is dropped. In the implementation the tangent is evaluated with its
argument clipped just inside ``\pm\pi/2``, which keeps the coefficient finite without
changing any admissible aperture.

### Reference bus

The polar condition ``\theta_{\mathrm{ref}} = 0`` becomes the vanishing of the imaginary
part of the reference voltage:

```math
v^{\mathrm{im}}_{\mathrm{ref},t} = 0, \qquad t\in\mathcal{T}. \tag{26}
```

### Generation constraints

Time envelope (added only when the period limit is below the static limit) and renewable
curtailment, identical to the polar model:

```math
p^g_{g,t} \le \overline{P}^{ts}_{g,t},\quad g\in\mathcal{G}\setminus\mathcal{G}^c; \qquad
p^g_{g,t} \ge \underline{P}^{ts}_{g,t},\quad g\in\mathcal{G}, \tag{27}
```
```math
\underline{Q}^{ts}_{g,t} \le Q^g_{g,t} \le \overline{Q}^{ts}_{g,t},\quad g\in\mathcal{G}; \qquad
p^g_{g,t} + p^c_{g,t} = \overline{P}^{ts}_{g,t},\quad g\in\mathcal{G}^c, \tag{28}
```

for all ``t\in\mathcal{T}``.

### Storage constraints

Identical to the polar formulation, except that the connection-bus voltage appearing in the
converter losses is written as ``(v^{\mathrm{re}}_{B(s),t})^2 + (v^{\mathrm{im}}_{B(s),t})^2``
instead of ``V_{B(s),t}^2``.

**Energy balance.**

```math
E_{s,1} = E^0_s + \eta^c_s\,c_{s,1} - \frac{d_{s,1}}{\eta^d_s},\quad s\in\mathcal{S}; \qquad
E_{s,t} = E_{s,t-1} + \eta^c_s\,c_{s,t} - \frac{d_{s,t}}{\eta^d_s},\quad t = 2,\dots,T. \tag{29}
```

**Non-simultaneity of charging and discharging.**

```math
c_{s,t}\,d_{s,t} = 0, \qquad s\in\mathcal{S},\ t\in\mathcal{T}. \tag{30}
```

**Inverter thermal limit.**

```math
(P^s_{s,t})^2 + (Q^s_{s,t})^2 \le \overline{S}_s^2, \qquad s\in\mathcal{S},\ t\in\mathcal{T}. \tag{31}
```

**Converter losses (AC coupling).**

```math
P^s_{s,t} + d_{s,t} - c_{s,t}
= P^{loss}_s + r_s\,\frac{(P^s_{s,t})^2 + (Q^s_{s,t})^2}{(v^{\mathrm{re}}_{B(s),t})^2 + (v^{\mathrm{im}}_{B(s),t})^2},
\qquad s\in\mathcal{S},\ t\in\mathcal{T}, \tag{32}
```
```math
Q^s_{s,t} = Q^{sc}_{s,t} + Q^{loss}_s
+ x_s\,\frac{(P^s_{s,t})^2 + (Q^s_{s,t})^2}{(v^{\mathrm{re}}_{B(s),t})^2 + (v^{\mathrm{im}}_{B(s),t})^2},
\qquad s\in\mathcal{S},\ t\in\mathcal{T}. \tag{33}
```

### Complete model

```math
\min_{\substack{v^{\mathrm{re}},\,v^{\mathrm{im}},\,p^g,\,Q^g,\,p^c,\\ E,\,c,\,d,\,P^s,\,Q^s,\,Q^{sc}}}
\quad F + F^{curt} \qquad (\text{Eq. } 14)
```

subject to:

| Block | Equations |
|---|---|
| Flows | (15)–(20) |
| Balance | (21), (22) |
| Network | (23)–(26) |
| Generation | (27), (28) |
| Storage | (29)–(33) |
| Bounds | (1)–(11) |

