# IVR

The IVR model (current–voltage, rectangular) is the third exact formulation. In the
previous models, [AC Polar](@ref) and [AC Rectangular (ACR)](@ref), the branch current
was computed from the admittance matrix and **immediately eliminated** by substituting it
into ``S = V\,I^{*}`` to obtain powers. The current was an intermediate step that
disappeared, and the only decision variables were voltages. The IVR model does the
opposite. Instead of eliminating the current, it **promotes it to a decision variable**.
Each branch gains four new variables, the real and imaginary parts of the current at each
port,

- ``c^{fr}_r,\ c^{fr}_i``: current at the *from* port;
- ``c^{to}_r,\ c^{to}_i``: current at the *to* port,

while the voltages stay rectangular, ``v^{\mathrm{re}}, v^{\mathrm{im}}``. Because current
and voltage are now both variables, the relation between them, Ohm's law, must now be
written as a constraint. It is no longer a substitution but an equation of the model. The
consequence is structural. Ohm's law ``I = Y\,V`` is **linear** in ``V``, so keeping the
current as a variable turns it into a set of **linear equalities**. The model stays exact
and nonconvex, but the nonconvexity is pushed out of the network equations into three
places, as shown at the end.

## Derivation

### Ohm's law in rectangular coordinates

Start from the branch admittance matrix already derived for the [AC Polar](@ref) model:

```math
\begin{bmatrix} I_i \\ I_j \end{bmatrix}
=
\begin{bmatrix} Y_{ff} & Y_{ft} \\ Y_{tf} & Y_{tt} \end{bmatrix}
\begin{bmatrix} V_i \\ V_j \end{bmatrix},
\qquad
\begin{aligned}
I_i &= Y_{ff}\,V_i + Y_{ft}\,V_j, \\
I_j &= Y_{tf}\,V_i + Y_{tt}\,V_j.
\end{aligned}
```

Each block is complex (``Y = Y^{\mathrm{re}} + j\,Y^{\mathrm{im}}``) and so is each voltage
(``V = v^{\mathrm{re}} + j\,v^{\mathrm{im}}``). Expanding the *from*-port current ``I_i`` into real and imaginary
parts,

```math
c^{fr}_r + j\,c^{fr}_i
= (Y^{\mathrm{re}}_{ff} + j\,Y^{\mathrm{im}}_{ff})(v^{\mathrm{re}}_i + j\,v^{\mathrm{im}}_i)
+ (Y^{\mathrm{re}}_{ft} + j\,Y^{\mathrm{im}}_{ft})(v^{\mathrm{re}}_j + j\,v^{\mathrm{im}}_j),
```

and separating the two parts gives the *from*-port current:

```math
c^{fr}_r = Y^{\mathrm{re}}_{ff}\,v^{\mathrm{re}}_i - Y^{\mathrm{im}}_{ff}\,v^{\mathrm{im}}_i
+ Y^{\mathrm{re}}_{ft}\,v^{\mathrm{re}}_j - Y^{\mathrm{im}}_{ft}\,v^{\mathrm{im}}_j,
```
```math
c^{fr}_i = Y^{\mathrm{im}}_{ff}\,v^{\mathrm{re}}_i + Y^{\mathrm{re}}_{ff}\,v^{\mathrm{im}}_i
+ Y^{\mathrm{im}}_{ft}\,v^{\mathrm{re}}_j + Y^{\mathrm{re}}_{ft}\,v^{\mathrm{im}}_j .
```

The same expansion of ``I_j`` gives the *to*-port current:

```math
c^{to}_r = Y^{\mathrm{re}}_{tf}\,v^{\mathrm{re}}_i - Y^{\mathrm{im}}_{tf}\,v^{\mathrm{im}}_i
+ Y^{\mathrm{re}}_{tt}\,v^{\mathrm{re}}_j - Y^{\mathrm{im}}_{tt}\,v^{\mathrm{im}}_j,
```
```math
c^{to}_i = Y^{\mathrm{im}}_{tf}\,v^{\mathrm{re}}_i + Y^{\mathrm{re}}_{tf}\,v^{\mathrm{im}}_i
+ Y^{\mathrm{im}}_{tt}\,v^{\mathrm{re}}_j + Y^{\mathrm{re}}_{tt}\,v^{\mathrm{im}}_j .
```

These are **four linear equalities per branch**, all of Ohm's law, with no nonlinearity.

### The Y coefficients

The blocks are exactly the ``Y`` blocks of the AC model. The tap sits on the *from* side, so
``Y_{ff}`` carries the factor ``1/tm`` while ``Y_{tt}`` is unscaled. For ``Y_{ft}`` and
``Y_{tf}`` the numerator and denominator are rationalized by the conjugate of the tap, which
turns the ``e^{\pm j\theta_{shift}}`` of the AC blocks into the real pair ``(t_r, t_i)``.
Split into real and imaginary parts,

```math
Y^{\mathrm{re}}_{ff} = \frac{g + g^{fr}}{tm}, \qquad\qquad
Y^{\mathrm{im}}_{ff} = \frac{b + b^{fr}}{tm},
```
```math
Y^{\mathrm{re}}_{ft} = -\frac{g\,t_r - b\,t_i}{tm}, \qquad\qquad
Y^{\mathrm{im}}_{ft} = -\frac{g\,t_i + b\,t_r}{tm},
```
```math
Y^{\mathrm{re}}_{tf} = -\frac{g\,t_r + b\,t_i}{tm}, \qquad\qquad
Y^{\mathrm{im}}_{tf} = \frac{g\,t_i - b\,t_r}{tm},
```
```math
Y^{\mathrm{re}}_{tt} = g + g^{to}, \qquad\qquad
Y^{\mathrm{im}}_{tt} = b + b^{to},
```

where ``tm = |T|^2 = \tau^2``.

### Current balance (Kirchhoff's current law)

The IVR does not write a power balance. It writes a **current balance**. Since the variable
is current, the Kirchhoff law imposed is the nodal (current) law:

```math
\sum I_{\mathrm{gen}} = \sum I_{\mathrm{branch}} + I_{\mathrm{load}} + I_{\mathrm{shunt}} .
```

The branch currents are the Ohm-law variables above, and the generated currents are chosen
by the solver. The load and shunt terms must be derived.

### The load term

The load is given in power, but the balance asks for current. Convert it through
``S = V\,I^{*}``:

```math
I_{\mathrm{load}} = \left(\frac{S}{V}\right)^{*} = \frac{S^{*}}{V^{*}}
= \frac{p_d - j\,q_d}{v^{\mathrm{re}} - j\,v^{\mathrm{im}}} .
```

Rationalizing by the conjugate of the denominator,

```math
I_{\mathrm{load}} = \frac{(p_d - j\,q_d)(v^{\mathrm{re}} + j\,v^{\mathrm{im}})}{(v^{\mathrm{re}})^2 + (v^{\mathrm{im}})^2},
```

and separating real and imaginary parts:

```math
I_{\mathrm{load}} = \underbrace{\frac{p_d\,v^{\mathrm{re}} + q_d\,v^{\mathrm{im}}}{(v^{\mathrm{re}})^2 + (v^{\mathrm{im}})^2}}_{\text{real part}}
+ j\,\underbrace{\frac{p_d\,v^{\mathrm{im}} - q_d\,v^{\mathrm{re}}}{(v^{\mathrm{re}})^2 + (v^{\mathrm{im}})^2}}_{\text{imaginary part}} .
```

The division by ``(v^{\mathrm{re}})^2 + (v^{\mathrm{im}})^2`` is the first place where the
constant-power load reintroduces a nonconvexity.

### The shunt term

The bus shunt is an admittance ``Y_s = g_s + j\,b_s`` connected from the bus to ground. The
current it draws is pure Ohm's law (voltage times admittance):

```math
I_{\mathrm{shunt}} = Y_s\,V = (g_s + j\,b_s)(v^{\mathrm{re}} + j\,v^{\mathrm{im}})
= \underbrace{(g_s\,v^{\mathrm{re}} - b_s\,v^{\mathrm{im}})}_{\text{real part}}
+ j\,\underbrace{(g_s\,v^{\mathrm{im}} + b_s\,v^{\mathrm{re}})}_{\text{imaginary part}} .
```

### The generator bridge: power ↔ current

The balance is in current, so the generator injects a current (``c_{rg}, c_{ig}``). But its
limits are natural in power (``P_g^{\min} \le p_g \le P_g^{\max}``), and the cost in the
objective is also a function of ``p_g``. A bridge linking the generator's current to its
power comes from ``S_g = V\,I_g^{*}`` applied to the generator:

```math
S_g = V\,I_g^{*} = (v^{\mathrm{re}} + j\,v^{\mathrm{im}})(c_{rg} - j\,c_{ig})
= \underbrace{(v^{\mathrm{re}}\,c_{rg} + v^{\mathrm{im}}\,c_{ig})}_{p_g}
+ j\,\underbrace{(v^{\mathrm{im}}\,c_{rg} - v^{\mathrm{re}}\,c_{ig})}_{q_g} .
```

Hence

```math
p_g = v^{\mathrm{re}}\,c_{rg} + v^{\mathrm{im}}\,c_{ig}, \qquad
q_g = v^{\mathrm{im}}\,c_{rg} - v^{\mathrm{re}}\,c_{ig} .
```

### Main constraints

**Thermal limit.** It is imposed on apparent power, ``|S| \le rate_a``. The IVR has no
power variable, only voltage and current, so ``S`` is converted into them. From
``S = V\,I^{*}``, taking the modulus,

```math
|S| = |V\,I^{*}| = |V|\,|I^{*}| = |V|\,|I|, \qquad |S|^2 = |V|^2\,|I|^2,
```

which gives, per port, a **quartic** inequality:

```math
\big((v^{\mathrm{re}}_f)^2 + (v^{\mathrm{im}}_f)^2\big)\big((c^{fr}_r)^2 + (c^{fr}_i)^2\big) \le rate_a^2,
```
```math
\big((v^{\mathrm{re}}_t)^2 + (v^{\mathrm{im}}_t)^2\big)\big((c^{to}_r)^2 + (c^{to}_i)^2\big) \le rate_a^2 .
```

**Angle limit.** Identical to the [AC Rectangular (ACR)](@ref) model. It acts through the
tangent on the voltage cross terms
``W_r = v^{\mathrm{re}}_i v^{\mathrm{re}}_j + v^{\mathrm{im}}_i v^{\mathrm{im}}_j`` and
``W_i = v^{\mathrm{im}}_i v^{\mathrm{re}}_j - v^{\mathrm{re}}_i v^{\mathrm{im}}_j``, written
``C_{\ell,t}`` and ``S_{\ell,t}`` in the formulation:

```math
W_r\tan(\theta^{\min}) \le W_i \le \tan(\theta^{\max})\,W_r .
```

### Where the nonconvexity migrates

In the IVR the nonconvexity does not disappear. It **migrates**. The difficulty
of ``S = V\,I^{*}`` moves to three places where power meets current:

1. **Thermal limit** → becomes quartic (degree-4).
2. **Load in the balance** → becomes a nonconvex ratio.
3. **Generator bridge** (power ↔ current) → becomes bilinear.

Ohm's law itself is now linear. The price is that these three couplings carry the
nonconvexity that the voltage-only models hid inside the flow expressions.

---

## IVR Formulation

The IVR reuses the sets, parameters, objective, generation and storage blocks of the
[AC Rectangular (ACR)](@ref) model verbatim. What changes is the network layer: each
branch carries four **current** variables tied to the voltages by the linear Ohm's law
(Eqs. 21–24), each generator carries two current variables tied to its power by the
bilinear bridge (Eqs. 27–28), and the nodal law is a **current** balance (Eqs. 25–26)
instead of a power balance.

### Sets

- ``\mathcal{N}``: set of buses.
- ``\mathcal{N}^{nv} \subseteq \mathcal{N}``: non-virtual buses, on which the nodal current balance is enforced.
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
- ``\mathcal{G}(b)``: generators connected to bus ``b``.
- ``B(g)``: connection bus of generator ``g``.
- ``B(s)``: connection bus of storage ``s``.

### Parameters

**Branches.**

- ``g_\ell,\ b_\ell``: series conductance and susceptance of branch ``\ell``.
- ``g^{fr}_\ell,\ b^{fr}_\ell``: shunt conductance and susceptance, *from* side.
- ``g^{to}_\ell,\ b^{to}_\ell``: shunt conductance and susceptance, *to* side.
- ``k_\ell``: inverse tap ratio (``1/\mathrm{tap}``), so ``k_\ell^2 = 1/tm = 1/|T_\ell|^2``.
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

**Voltage** (rectangular, as in the ACR model).

```math
v^{\mathrm{re}}_{b,t} \in [-V_b^{\max}, V_b^{\max}] \;\text{: real part of the voltage at bus } b. \tag{1}
```

```math
v^{\mathrm{im}}_{b,t} \in [-V_b^{\max}, V_b^{\max}] \;\text{: imaginary part of the voltage at bus } b. \tag{2}
```

**Branch currents** (the new variables that define the IVR), in code
``c^{fr}_r,c^{fr}_i,c^{to}_r,c^{to}_i`` = `cr_fr`, `ci_fr`, `cr_to`, `ci_to`. They are
*free*, with no bound on the current itself. The apparent-power limit is carried by
the thermal constraint (Eqs. 30–31).

```math
c^{fr}_{r,\ell,t} \in \mathbb{R} \;\text{: real part of the current at the from port.} \tag{3}
```

```math
c^{fr}_{i,\ell,t} \in \mathbb{R} \;\text{: imaginary part of the current at the from port.} \tag{4}
```

```math
c^{to}_{r,\ell,t} \in \mathbb{R} \;\text{: real part of the current at the to port.} \tag{5}
```

```math
c^{to}_{i,\ell,t} \in \mathbb{R} \;\text{: imaginary part of the current at the to port.} \tag{6}
```

**Generation.**

```math
p^g_{g,t} \in [P_g^{\min}, P_g^{\max}] \;\text{: active power generated by } g. \tag{7}
```

```math
Q^g_{g,t} \in [Q_g^{\min}, Q_g^{\max}] \;\text{: reactive power generated by } g. \tag{8}
```

```math
p^c_{g,t} \in [0, \overline{P}^{ts}_{g,t}] \;\text{: curtailed renewable power, } g \in \mathcal{G}^c. \tag{9}
```

**Generator currents** (free, one pair per existing generator), in code
``c^{\mathrm{g}}_r, c^{\mathrm{g}}_i`` = `crg`, `cig`.

```math
c^{\mathrm{g}}_{r,g,t} \in \mathbb{R} \;\text{: real part of the generator current.} \tag{10}
```

```math
c^{\mathrm{g}}_{i,g,t} \in \mathbb{R} \;\text{: imaginary part of the generator current.} \tag{11}
```

**Storage** (identical to the ACR / AC storage variables).

```math
E_{s,t} \in [E_s^{\min}, E_s^{\max}] \;\text{: energy stored in } s \text{ at the end of } t. \tag{12}
```

```math
c_{s,t} \in [0, C_s^{\max}] \;\text{: charging power of } s. \tag{13}
```

```math
d_{s,t} \in [0, D_s^{\max}] \;\text{: discharging power of } s. \tag{14}
```

```math
P^s_{s,t} \in [-\overline{S}_s, \overline{S}_s] \;\text{: net active power exchanged by } s \text{, positive when charging.} \tag{15}
```

```math
Q^s_{s,t} \in [Q_s^{\min}, Q_s^{\max}] \;\text{: net reactive power exchanged by } s \text{, positive when absorbed.} \tag{16}
```

```math
Q^{sc}_{s,t} \in [Q_s^{\min}, Q_s^{\max}] \;\text{: internal reactive compensation of the converter of } s. \tag{17}
```

### Objective function

Identical to the ACR and AC models. The fuel-cost polynomial is evaluated on the power in
physical units (``S^{base}\,p^g``), plus the curtailment cost:

```math
F = \sum_{t\in\mathcal{T}} \sum_{g\in\mathcal{G}} \sum_{k=1}^{m_{g,t}}
c_{g,t,k}\left(S^{base}\,p^g_{g,t}\right)^{m_{g,t} - k}. \tag{18}
```

```math
F^{curt} = \sum_{g\in\mathcal{G}^c} \sum_{t\in\mathcal{T}}
c^{curt}_g\,S^{base}\,p^c_{g,t}. \tag{19}
```

Hence:

```math
\min\; F + F^{curt}. \tag{20}
```

### Branch currents (Ohm's law)

Each ``Y`` block is the AC admittance block, split into real and imaginary parts. With
``k_\ell = 1/\mathrm{tap}`` and the complex tap ``(T^{\mathrm{re}}_\ell, T^{\mathrm{im}}_\ell)``:

```math
Y^{\mathrm{re}}_{ff} = (g_\ell + g^{fr}_\ell)\,k_\ell^2, \qquad
Y^{\mathrm{im}}_{ff} = (b_\ell + b^{fr}_\ell)\,k_\ell^2,
```

```math
Y^{\mathrm{re}}_{ft} = -(g_\ell T^{\mathrm{re}}_\ell - b_\ell T^{\mathrm{im}}_\ell)\,k_\ell^2, \qquad
Y^{\mathrm{im}}_{ft} = -(g_\ell T^{\mathrm{im}}_\ell + b_\ell T^{\mathrm{re}}_\ell)\,k_\ell^2,
```

```math
Y^{\mathrm{re}}_{tf} = -(g_\ell T^{\mathrm{re}}_\ell + b_\ell T^{\mathrm{im}}_\ell)\,k_\ell^2, \qquad
Y^{\mathrm{im}}_{tf} = (g_\ell T^{\mathrm{im}}_\ell - b_\ell T^{\mathrm{re}}_\ell)\,k_\ell^2,
```

```math
Y^{\mathrm{re}}_{tt} = g_\ell + g^{to}_\ell, \qquad
Y^{\mathrm{im}}_{tt} = b_\ell + b^{to}_\ell.
```

Writing ``v^{\mathrm{re}}_f \equiv v^{\mathrm{re}}_{f(\ell),t}``,
``v^{\mathrm{re}}_t \equiv v^{\mathrm{re}}_{t(\ell),t}`` (and likewise for ``v^{\mathrm{im}}``),
Ohm's law is the set of **four linear equalities per branch**:

```math
c^{fr}_{r,\ell,t} = Y^{\mathrm{re}}_{ff}\,v^{\mathrm{re}}_f - Y^{\mathrm{im}}_{ff}\,v^{\mathrm{im}}_f
+ Y^{\mathrm{re}}_{ft}\,v^{\mathrm{re}}_t - Y^{\mathrm{im}}_{ft}\,v^{\mathrm{im}}_t,
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{21}
```

```math
c^{fr}_{i,\ell,t} = Y^{\mathrm{im}}_{ff}\,v^{\mathrm{re}}_f + Y^{\mathrm{re}}_{ff}\,v^{\mathrm{im}}_f
+ Y^{\mathrm{im}}_{ft}\,v^{\mathrm{re}}_t + Y^{\mathrm{re}}_{ft}\,v^{\mathrm{im}}_t,
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{22}
```

```math
c^{to}_{r,\ell,t} = Y^{\mathrm{re}}_{tf}\,v^{\mathrm{re}}_f - Y^{\mathrm{im}}_{tf}\,v^{\mathrm{im}}_f
+ Y^{\mathrm{re}}_{tt}\,v^{\mathrm{re}}_t - Y^{\mathrm{im}}_{tt}\,v^{\mathrm{im}}_t,
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{23}
```

```math
c^{to}_{i,\ell,t} = Y^{\mathrm{im}}_{tf}\,v^{\mathrm{re}}_f + Y^{\mathrm{re}}_{tf}\,v^{\mathrm{im}}_f
+ Y^{\mathrm{im}}_{tt}\,v^{\mathrm{re}}_t + Y^{\mathrm{re}}_{tt}\,v^{\mathrm{im}}_t,
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{24}
```

### Current balance (Kirchhoff's current law)

In place of a power balance, the IVR writes a **current** balance at each non-virtual bus
``b\in\mathcal{N}^{nv}``, split into its real and imaginary parts. The generator currents inject into the bus. The
branch currents (summed over the branches for which ``b`` is the *from* or the *to* end),
the constant-power load current and the shunt current flow out:

```math
\sum_{g\in\mathcal{G}(b)} c^{\mathrm{g}}_{r,g,t}
= \sum_{\ell:\,f(\ell)=b} c^{fr}_{r,\ell,t}
+ \sum_{\ell:\,t(\ell)=b} c^{to}_{r,\ell,t}
+ \frac{P^d_{b,t}\,v^{\mathrm{re}}_{b,t} + Q^d_{b,t}\,v^{\mathrm{im}}_{b,t}}
{(v^{\mathrm{re}}_{b,t})^2 + (v^{\mathrm{im}}_{b,t})^2}
+ g^{sh}_b\,v^{\mathrm{re}}_{b,t} - b^{sh}_b\,v^{\mathrm{im}}_{b,t},
\qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{25}
```

```math
\sum_{g\in\mathcal{G}(b)} c^{\mathrm{g}}_{i,g,t}
= \sum_{\ell:\,f(\ell)=b} c^{fr}_{i,\ell,t}
+ \sum_{\ell:\,t(\ell)=b} c^{to}_{i,\ell,t}
+ \frac{P^d_{b,t}\,v^{\mathrm{im}}_{b,t} - Q^d_{b,t}\,v^{\mathrm{re}}_{b,t}}
{(v^{\mathrm{re}}_{b,t})^2 + (v^{\mathrm{im}}_{b,t})^2}
+ g^{sh}_b\,v^{\mathrm{im}}_{b,t} + b^{sh}_b\,v^{\mathrm{re}}_{b,t},
\qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{26}
```

### Generator power–current coupling

The bridge ``S_g = V\,I_g^{*}`` links each generator's current to its power at the
connection bus ``B(g)``. It is bilinear in ``(v, c^{\mathrm{g}})``:

```math
p^g_{g,t} = v^{\mathrm{re}}_{B(g),t}\,c^{\mathrm{g}}_{r,g,t}
+ v^{\mathrm{im}}_{B(g),t}\,c^{\mathrm{g}}_{i,g,t},
\qquad g\in\mathcal{G},\ t\in\mathcal{T}. \tag{27}
```

```math
Q^g_{g,t} = v^{\mathrm{im}}_{B(g),t}\,c^{\mathrm{g}}_{r,g,t}
- v^{\mathrm{re}}_{B(g),t}\,c^{\mathrm{g}}_{i,g,t},
\qquad g\in\mathcal{G},\ t\in\mathcal{T}. \tag{28}
```

### Network constraints

**Voltage magnitude limits.** Imposed on the sum of squares, so each bus voltage lives in a
ring of inner radius ``V_b^{\min}`` and outer radius ``V_b^{\max}``:

```math
(V_b^{\min})^2 \le (v^{\mathrm{re}}_{b,t})^2 + (v^{\mathrm{im}}_{b,t})^2 \le (V_b^{\max})^2,
\qquad b\in\mathcal{N},\ t\in\mathcal{T}. \tag{29}
```

**Branch thermal limits.** With no power variable, ``|S|^2 = |V|^2\,|I|^2`` makes the
apparent-power limit a **quartic** inequality per port:

```math
\big((v^{\mathrm{re}}_f)^2 + (v^{\mathrm{im}}_f)^2\big)
\big((c^{fr}_{r,\ell,t})^2 + (c^{fr}_{i,\ell,t})^2\big) \le \overline{S}_\ell^2,
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{30}
```

```math
\big((v^{\mathrm{re}}_t)^2 + (v^{\mathrm{im}}_t)^2\big)
\big((c^{to}_{r,\ell,t})^2 + (c^{to}_{i,\ell,t})^2\big) \le \overline{S}_\ell^2,
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{31}
```

**Angular aperture limits.** Imposed through the tangent on the voltage cross terms, as in the
ACR model. With

```math
C_{\ell,t} = v^{\mathrm{re}}_f v^{\mathrm{re}}_t + v^{\mathrm{im}}_f v^{\mathrm{im}}_t, \qquad
S_{\ell,t} = v^{\mathrm{im}}_f v^{\mathrm{re}}_t - v^{\mathrm{re}}_f v^{\mathrm{im}}_t,
```

the aperture bounds are the pair of inequalities

```math
\tan(\delta_\ell^{\min})\,C_{\ell,t} \le S_{\ell,t} \le \tan(\delta_\ell^{\max})\,C_{\ell,t},
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{32}
```

Each side is added only when its limit is set: the lower bound when ``\delta_\ell^{\min}`` is
defined, the upper bound when ``\delta_\ell^{\max}`` is defined.

### Reference bus

The polar condition ``\theta_{\mathrm{ref}} = 0`` becomes the vanishing of the imaginary
part of the reference voltage:

```math
v^{\mathrm{im}}_{\mathrm{ref},t} = 0, \qquad t\in\mathcal{T}. \tag{33}
```

### Generation constraints

Active/reactive time envelope (each row added only when the period limit is strictly
tighter than the static bound) and renewable curtailment, identical to the ACR model:

```math
p^g_{g,t} \le \overline{P}^{ts}_{g,t},\quad g\in\mathcal{G}\setminus\mathcal{G}^c; \tag{34}
```

```math
p^g_{g,t} \ge \underline{P}^{ts}_{g,t},\quad g\in\mathcal{G}; \tag{35}
```

```math
\underline{Q}^{ts}_{g,t} \le Q^g_{g,t} \le \overline{Q}^{ts}_{g,t},\quad g\in\mathcal{G}; \tag{36}
```

```math
p^g_{g,t} + p^c_{g,t} = \overline{P}^{ts}_{g,t},\quad g\in\mathcal{G}^c, \tag{37}
```

for all ``t\in\mathcal{T}``.

### Storage constraints

The IVR dispatches the **same storage block as the ACR model** (converter losses use the
rectangular connection-bus voltage ``(v^{\mathrm{re}}_{B(s),t})^2 + (v^{\mathrm{im}}_{B(s),t})^2``).

!!! note "Storage coupling in IVR"
    The constraints below bind each unit's own ``P^s`` and ``Q^s``. In the current
    implementation these are not injected into the current balance (25)–(26), so an IVR
    storage unit is modeled locally and does not exchange at its bus. No case in this
    work uses storage, and wiring the storage current into the balance is left to a
    future issue.

**Energy balance.**

```math
E_{s,1} = E^0_s + \eta^c_s\,c_{s,1} - \frac{d_{s,1}}{\eta^d_s},\quad s\in\mathcal{S}; \tag{38}
```

```math
E_{s,t} = E_{s,t-1} + \eta^c_s\,c_{s,t} - \frac{d_{s,t}}{\eta^d_s},\quad s\in\mathcal{S},\ t = 2,\dots,T. \tag{39}
```

**Non-simultaneity of charging and discharging.**

```math
c_{s,t}\,d_{s,t} = 0, \qquad s\in\mathcal{S},\ t\in\mathcal{T}. \tag{40}
```

**Inverter thermal limit.**

```math
(P^s_{s,t})^2 + (Q^s_{s,t})^2 \le \overline{S}_s^2, \qquad s\in\mathcal{S},\ t\in\mathcal{T}. \tag{41}
```

**Converter losses (AC coupling).**

```math
P^s_{s,t} + d_{s,t} - c_{s,t}
= P^{loss}_s + r_s\,\frac{(P^s_{s,t})^2 + (Q^s_{s,t})^2}
{(v^{\mathrm{re}}_{B(s),t})^2 + (v^{\mathrm{im}}_{B(s),t})^2},
\qquad s\in\mathcal{S},\ t\in\mathcal{T}. \tag{42}
```

```math
Q^s_{s,t} = Q^{sc}_{s,t} + Q^{loss}_s
+ x_s\,\frac{(P^s_{s,t})^2 + (Q^s_{s,t})^2}
{(v^{\mathrm{re}}_{B(s),t})^2 + (v^{\mathrm{im}}_{B(s),t})^2},
\qquad s\in\mathcal{S},\ t\in\mathcal{T}. \tag{43}
```

### Complete model

```math
\min_{\substack{v^{\mathrm{re}},\,v^{\mathrm{im}},\,c^{fr}_r,\,c^{fr}_i,\,c^{to}_r,\,c^{to}_i,\\
p^g,\,Q^g,\,p^c,\,c^{\mathrm{g}}_r,\,c^{\mathrm{g}}_i,\\
E,\,c,\,d,\,P^s,\,Q^s,\,Q^{sc}}}
\quad F + F^{curt} \qquad (\text{Eq. } 20)
```

subject to:

| Block | Equations |
|---|---|
| Ohm's law (branch currents) | (21)–(24) |
| Current balance | (25), (26) |
| Generator coupling | (27), (28) |
| Network | (29)–(32) |
| Reference bus | (33) |
| Generation | (34)–(37) |
| Storage | (38)–(43) |
| Bounds | (1)–(17) |

The power–current coupling ``S = V\,I^{*}`` carries its nonconvexity into three places, as
anticipated by the derivation: the quartic thermal limits (30)–(31), the constant-power load
ratio in the current balance (25)–(26), and the bilinear generator coupling (27)–(28). Ohm's
law (21)–(24) is linear. The voltage lower bound (29) and the angle limits (32) stay
nonconvex as well, exactly as in the ACR model. Promoting the current neither adds to nor
removes those.
