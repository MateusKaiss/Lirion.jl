# AC Polar

The AC model is the exact, nonlinear representation of power flow, formulated in
**polar coordinates**, with the nodal voltages described by a magnitude ``V_{b,t}`` and an
angle ``\theta_{b,t}``. It minimizes the generation cost subject to the network equations
written **without any simplification**.

## Derivation

Deriving the complete model requires three physical laws.

### 1. Complex power injection

The power flowing from a node into a branch is

```math
S = V \cdot I^{*},
```

where ``V`` is the phasor voltage at the node and ``I`` the current leaving it into the
branch.

### 2. Complex Ohm's law

The branch current is

```math
I = Y\,(V_i - V_j),
```

where ``Y`` is the series admittance connecting buses ``i`` and ``j``, ``V_i`` is the
voltage at bus ``i`` and ``V_j`` the voltage at bus ``j``.

### 3. π model

A real line is not merely a series impedance: it has distributed capacitance (the line
"charging"), modeled as shunt admittances ``y^{fr}``, ``y^{to}`` at the ends, and it may
include a transformer with a **complex tap** on the *from* side. A branch ``\ell``
connecting buses ``i`` and ``j`` and containing an ideal transformer is thus represented
by:

```@raw html
<figure style="text-align:center;margin:1.6rem 0">
<svg viewBox="0 0 525 225" role="img" aria-label="Pi-model of a branch with an ideal transformer on the from side" style="max-width:520px;width:100%;height:auto;color:inherit" xmlns="http://www.w3.org/2000/svg" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round" font-family="Georgia, 'Times New Roman', serif">
  <!-- bus i -->
  <line x1="55" y1="50" x2="55" y2="185"/>
  <line x1="55" y1="50" x2="41" y2="36"/>
  <text x="31" y="44" fill="currentColor" stroke="none" font-style="italic" font-size="20">i</text>
  <!-- bus i to transformer -->
  <line x1="55" y1="120" x2="90" y2="120"/>
  <!-- ideal transformer -->
  <circle cx="107" cy="120" r="17"/>
  <circle cx="131" cy="120" r="17"/>
  <!-- transformer to node A -->
  <line x1="148" y1="120" x2="233" y2="120"/>
  <circle cx="233" cy="120" r="3.5" fill="currentColor" stroke="none"/>
  <!-- shunt y^fr -->
  <line x1="233" y1="120" x2="233" y2="151"/>
  <rect x="221" y="151" width="24" height="33"/>
  <line x1="233" y1="184" x2="233" y2="197"/>
  <line x1="220" y1="197" x2="246" y2="197"/>
  <line x1="224" y1="201.5" x2="242" y2="201.5"/>
  <line x1="228" y1="206" x2="238" y2="206"/>
  <text x="253" y="173" fill="currentColor" stroke="none" font-style="italic" font-size="16">y<tspan font-size="11" dy="-6">fr</tspan></text>
  <!-- node A to series box -->
  <line x1="233" y1="120" x2="270" y2="120"/>
  <!-- series admittance -->
  <rect x="270" y="105" width="78" height="30"/>
  <text x="309" y="99" fill="currentColor" stroke="none" font-style="italic" font-size="18" text-anchor="middle">y</text>
  <!-- series box to node B -->
  <line x1="348" y1="120" x2="385" y2="120"/>
  <circle cx="385" cy="120" r="3.5" fill="currentColor" stroke="none"/>
  <!-- shunt y^to -->
  <line x1="385" y1="120" x2="385" y2="151"/>
  <rect x="373" y="151" width="24" height="33"/>
  <line x1="385" y1="184" x2="385" y2="197"/>
  <line x1="372" y1="197" x2="398" y2="197"/>
  <line x1="376" y1="201.5" x2="394" y2="201.5"/>
  <line x1="380" y1="206" x2="390" y2="206"/>
  <text x="405" y="173" fill="currentColor" stroke="none" font-style="italic" font-size="16">y<tspan font-size="11" dy="-6">to</tspan></text>
  <!-- node B to bus j -->
  <line x1="385" y1="120" x2="470" y2="120"/>
  <!-- bus j -->
  <line x1="470" y1="50" x2="470" y2="185"/>
  <line x1="470" y1="50" x2="484" y2="36"/>
  <text x="488" y="44" fill="currentColor" stroke="none" font-style="italic" font-size="20">j</text>
</svg>
<figcaption style="font-size:0.9em;opacity:0.72;margin-top:0.5rem">
&pi;-model of a branch \(\ell\) connecting buses \(i\) and \(j\), with an ideal transformer and shunt admittances \(y^{fr}\), \(y^{to}\).
</figcaption>
</figure>
```

The voltage in polar coordinates is the heart of the model:

```math
V_i = v_i\,e^{j\theta_i} = v_i\,(\cos\theta_i + j\,\sin\theta_i),
```

where ``v_i = |V_i|`` is the magnitude and ``\theta_i`` the angle. Geometrically, it is a
point in the complex plane whose projections are ``v_i\cos\theta_i`` (real axis) and
``v_i\sin\theta_i`` (imaginary axis):

```@raw html
<figure style="text-align:center;margin:1.6rem 0">
<svg viewBox="0 0 390 285" role="img" aria-label="Voltage phasor in the complex plane" style="max-width:380px;width:100%;height:auto;color:inherit" xmlns="http://www.w3.org/2000/svg" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round" font-family="Georgia, 'Times New Roman', serif">
  <!-- Im axis -->
  <line x1="80" y1="235" x2="80" y2="45"/>
  <line x1="80" y1="45" x2="74" y2="57"/>
  <line x1="80" y1="45" x2="86" y2="57"/>
  <text x="58" y="42" fill="currentColor" stroke="none" font-style="italic" font-size="18">Im</text>
  <!-- Re axis -->
  <line x1="80" y1="235" x2="360" y2="235"/>
  <line x1="360" y1="235" x2="348" y2="229"/>
  <line x1="360" y1="235" x2="348" y2="241"/>
  <text x="362" y="240" fill="currentColor" stroke="none" font-style="italic" font-size="18">Re</text>
  <!-- projections (dashed) -->
  <g stroke-dasharray="4 3" stroke-width="1.6">
    <line x1="235" y1="120" x2="235" y2="235"/>
    <line x1="235" y1="120" x2="80" y2="120"/>
  </g>
  <!-- angle arc theta_i -->
  <path d="M 125 235 A 45 45 0 0 0 116 208" fill="none" stroke-width="1.8"/>
  <text x="129" y="224" fill="currentColor" stroke="none" font-style="italic" font-size="15">&theta;<tspan font-size="11" dy="4">i</tspan></text>
  <!-- voltage phasor -->
  <line x1="80" y1="235" x2="235" y2="120" stroke-width="2.8"/>
  <circle cx="235" cy="120" r="4" fill="currentColor" stroke="none"/>
  <text x="138" y="156" fill="currentColor" stroke="none" font-style="italic" font-size="17">V<tspan font-size="11" dy="4">i</tspan></text>
  <!-- projection labels -->
  <text x="210" y="256" fill="currentColor" stroke="none" font-style="italic" font-size="15">v<tspan font-size="11" dy="-6">re</tspan><tspan font-size="11" dy="10">i</tspan></text>
  <text x="14" y="125" fill="currentColor" stroke="none" font-style="italic" font-size="15">v<tspan font-size="11" dy="-6">im</tspan><tspan font-size="11" dy="10">i</tspan></text>
</svg>
<figcaption style="font-size:0.9em;opacity:0.72;margin-top:0.5rem">
Voltage phasor \(V_i = v_i\,e^{j\theta_i}\) in the complex plane, with real and imaginary projections \(v^{\mathrm{re}}_i\), \(v^{\mathrm{im}}_i\).
</figcaption>
</figure>
```

### Series admittance from ``r`` and ``x``

The branch has resistance ``r`` and reactance ``x``, hence series impedance ``z = r + jx``.
The admittance is its inverse:

```math
y = \frac{1}{z} = \frac{1}{r + jx} = \frac{r - jx}{(r + jx)(r - jx)} = \frac{r - jx}{r^2 + x^2}.
```

Separating real and imaginary parts,

```math
g = \frac{r}{r^2 + x^2}, \qquad b = -\,\frac{x}{r^2 + x^2}, \qquad y = g + jb,
```

where ``g`` is the conductance, ``b`` the susceptance and ``y`` the series admittance.

### Complex transformer tap

The transformer sits on the *from* end, with a complex transformation ratio

```math
T = \tau\,e^{j\theta_{shift}}
  = \underbrace{\tau\cos\theta_{shift}}_{\tau_r}
  + j\,\underbrace{\tau\sin\theta_{shift}}_{\tau_i},
```

where ``\tau`` is the tap magnitude and ``\theta_{shift}`` the phase-shift angle. We further
define

```math
tm = |T|^2 = \tau^2.
```

### Branch admittance matrix

**Ideal transformation relation.** By power conservation across the ideal transformer
(``V_i\,I_i^{*} = V'\,I'^{*}``):

```math
V' = \frac{V_i}{T}, \qquad I_i = \frac{I'}{T^{*}},
```

where ``V'`` and ``I'`` are the voltage and current at the internal node.

**Kirchhoff at the internal node.** From the diagram, the current entering the internal
node is the one flowing through the local shunt plus the one crossing the series branch;
the same holds at node ``j``:

```math
I' = y^{fr}\,V' + y\,(V' - V_j), \qquad
I_j = y^{to}\,V_j + y\,(V_j - V').
```

Collecting ``I'`` and substituting ``V' = V_i/T`` into ``T^{*} I_i = I'``:

```math
I' = (y + y^{fr})\,V' - y\,V_j
\;\Longrightarrow\;
T^{*} I_i = (y + y^{fr})\,\frac{V_i}{T} - y\,V_j,
```

and isolating ``I_i`` (and, analogously, ``I_j``):

```math
I_i = \frac{(y + y^{fr})\,V_i}{T\,T^{*}} - \frac{y}{T^{*}}\,V_j, \qquad
I_j = -\frac{y}{T}\,V_i + (y + y^{to})\,V_j.
```

Identifying the matrix entries,

```math
Y_{ff} = \frac{y + y^{fr}}{tm}, \qquad
Y_{ft} = -\frac{y}{T^{*}}, \qquad
Y_{tf} = -\frac{y}{T}, \qquad
Y_{tt} = y + y^{to},
```

we obtain the branch admittance matrix:

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

### Power on the *from* side

The total power leaving bus ``i`` splits into a self term and a cross term:

```math
S_{ij} = V_i\,I_i^{*} = V_i\,(Y_{ff}\,V_i + Y_{ft}\,V_j)^{*}
= \underbrace{V_i\,Y_{ff}^{*}\,V_i^{*}}_{\text{self term}}
+ \underbrace{V_i\,Y_{ft}^{*}\,V_j^{*}}_{\text{cross term}}.
```

**Self term.** With ``V_i\,V_i^{*} = v_i^2`` and
``Y_{ff}^{*} = \dfrac{(g+g^{fr}) - j(b+b^{fr})}{tm}``:

```math
V_i\,Y_{ff}^{*}\,V_i^{*} = \frac{(g+g^{fr})\,v_i^2}{tm} - j\,\frac{(b+b^{fr})\,v_i^2}{tm}.
```

**Cross term.** With ``Y_{ft}^{*} = -\dfrac{y^{*}}{T} = -\dfrac{y^{*}}{\tau}\,e^{-j\theta_{shift}}``,
``V_i V_j^{*} = v_i v_j\,e^{j(\theta_i - \theta_j)}`` and introducing the **effective angular
aperture** ``\delta = \theta_i - \theta_j - \theta_{shift}``:

```math
V_i\,Y_{ft}^{*}\,V_j^{*} = -\frac{y^{*}}{\tau}\,v_i\,v_j\,e^{j\delta}.
```

Expanding with ``y^{*} = g - jb`` and ``e^{j\delta} = \cos\delta + j\sin\delta``:

```math
V_i\,Y_{ft}^{*}\,V_j^{*}
= -\frac{v_i\,v_j}{\tau}\Big[(g\cos\delta + b\sin\delta)
+ j\,(g\sin\delta - b\cos\delta)\Big].
```

Separating the real (active) and imaginary (reactive) parts:

```math
p_{ij} = \frac{(g+g^{fr})\,v_i^2}{tm} - \frac{v_i\,v_j}{\tau}\,(g\cos\delta + b\sin\delta),
```

```math
q_{ij} = -\frac{(b+b^{fr})\,v_i^2}{tm} - \frac{v_i\,v_j}{\tau}\,(g\sin\delta - b\cos\delta).
```

### Power on the *to* side

Analogously, the total power leaving bus ``j`` is

```math
S_{ji} = V_j\,I_j^{*} = V_j\,(Y_{tf}\,V_i + Y_{tt}\,V_j)^{*}
= \underbrace{V_j\,Y_{tt}^{*}\,V_j^{*}}_{\text{self term}}
+ \underbrace{V_j\,Y_{tf}^{*}\,V_i^{*}}_{\text{cross term}}.
```

**Self term.** With ``V_j\,V_j^{*} = v_j^2`` and ``Y_{tt}^{*} = (g+g^{to}) - j(b+b^{to})``:

```math
V_j\,Y_{tt}^{*}\,V_j^{*} = (g+g^{to})\,v_j^2 - j\,(b+b^{to})\,v_j^2.
```

**Cross term.** With ``Y_{tf}^{*} = -\dfrac{y^{*}}{T^{*}} = -\dfrac{y^{*}}{\tau}\,e^{j\theta_{shift}}``,
``V_j V_i^{*} = v_j v_i\,e^{j(\theta_j - \theta_i)}`` and ``-\delta = \theta_j - \theta_i + \theta_{shift}``:

```math
V_j\,Y_{tf}^{*}\,V_i^{*} = -\frac{y^{*}}{\tau}\,v_i\,v_j\,e^{-j\delta}.
```

Expanding with ``y^{*} = g - jb`` and ``e^{-j\delta} = \cos\delta - j\sin\delta``:

```math
V_j\,Y_{tf}^{*}\,V_i^{*}
= -\frac{v_i\,v_j}{\tau}\Big[(g\cos\delta - b\sin\delta)
+ j\,(-g\sin\delta - b\cos\delta)\Big],
```

whence

```math
p_{ji} = (g+g^{to})\,v_j^2 - \frac{v_j\,v_i}{\tau}\,(g\cos\delta - b\sin\delta),
```

```math
q_{ji} = -(b+b^{to})\,v_j^2 - \frac{v_i\,v_j}{\tau}\,(-g\sin\delta - b\cos\delta).
```

### Main constraints

**Thermal limit.** The complex power on the line is ``S_{ij} = P_{ij} + jQ_{ij}``, with
magnitude ``|S_{ij}| = \sqrt{P_{ij}^2 + Q_{ij}^2}``. The square root makes the constraint
nonsmooth (not differentiable at zero and numerically more expensive). Since both sides
are nonnegative, squaring preserves the inequality:

```math
p_{ij}^2 + q_{ij}^2 \le rate_a^2, \qquad p_{ji}^2 + q_{ji}^2 \le rate_a^2.
```

**Angle limit.**

```math
\theta^{\min} \le \theta_i - \theta_j \le \theta^{\max}.
```

---

## AC Polar Formulation

### Sets

- ``\mathcal{N}``: set of buses.
- ``\mathcal{N}^{nv} \subseteq \mathcal{N}``: non-virtual buses, those where the nodal balance is imposed. Virtual buses added when the network is split into zones are excluded.
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
- ``k_\ell``: inverse tap ratio (``1/\mathrm{tap}``), the derivation's ``1/\tau``, so ``k_\ell^2 = 1/tm``.
- ``\Phi_\ell``: transformer phase shift, the derivation's ``\theta_{shift}``.
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

**Voltage.**

```math
\theta_{b,t} \in \mathbb{R} \;\text{: voltage angle at bus } b \text{ at time } t. \tag{1}
```

```math
V_{b,t} \in [V_b^{\min}, V_b^{\max}] \;\text{: voltage magnitude at bus } b. \tag{2}
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

The fuel cost is evaluated on the power in physical units (``S^{base}\,p^g``), with the
polynomial written in decreasing order of degree: for a vector of ``m_{g,t}`` coefficients,
term ``k`` has exponent ``m_{g,t} - k``.

```math
F = \sum_{t\in\mathcal{T}} \sum_{g\in\mathcal{G}} \sum_{k=1}^{m_{g,t}}
c_{g,t,k}\left(S^{base}\,p^g_{g,t}\right)^{m_{g,t} - k}. \tag{12}
```

The curtailment cost, when curtailable generators are present:

```math
F^{curt} = \sum_{g\in\mathcal{G}^c} \sum_{t\in\mathcal{T}}
c^{curt}_g\,S^{base}\,p^c_{g,t}. \tag{13}
```

Hence:

```math
\min\; F + F^{curt}. \tag{14}
```

### Branch flows

Define the branch effective angular aperture, with the transformer phase shift already
subtracted:

```math
\delta_{\ell,t} = \theta_{f(\ell),t} - \theta_{t(\ell),t} - \Phi_\ell,
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{15}
```

The four flows (active and reactive, in each direction) are:

```math
P^f_{\ell,t} = (g_\ell + g^{fr}_\ell)\,(k_\ell\,V_{f(\ell),t})^2
- k_\ell\,V_{f(\ell),t}\,V_{t(\ell),t}\,
(g_\ell\cos\delta_{\ell,t} + b_\ell\sin\delta_{\ell,t}). \tag{16}
```

```math
P^t_{\ell,t} = (g_\ell + g^{to}_\ell)\,V_{t(\ell),t}^2
- k_\ell\,V_{t(\ell),t}\,V_{f(\ell),t}\,
(g_\ell\cos\delta_{\ell,t} - b_\ell\sin\delta_{\ell,t}). \tag{17}
```

```math
Q^f_{\ell,t} = -(b_\ell + b^{fr}_\ell)\,(k_\ell\,V_{f(\ell),t})^2
- k_\ell\,V_{f(\ell),t}\,V_{t(\ell),t}\,
(-b_\ell\cos\delta_{\ell,t} + g_\ell\sin\delta_{\ell,t}). \tag{18}
```

```math
Q^t_{\ell,t} = -(b_\ell + b^{to}_\ell)\,V_{t(\ell),t}^2
- k_\ell\,V_{t(\ell),t}\,V_{f(\ell),t}\,
(-b_\ell\cos\delta_{\ell,t} - g_\ell\sin\delta_{\ell,t}). \tag{19}
```

### Power balance

The physical idea is simple and inviolable: at each bus, everything that comes in must go
out.

**Active power.**

```math
\sum_{g\in\mathcal{G}(b)} p^g_{g,t}
- \sum_{s\in\mathcal{S}(b)} P^s_{s,t}
- P^d_{b,t}
- \sum_{\ell:\,f(\ell)=b} P^f_{\ell,t}
- \sum_{\ell:\,t(\ell)=b} P^t_{\ell,t}
- g^{sh}_b\,V_{b,t}^2 = 0,
\qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{20}
```

**Reactive power.**

```math
\sum_{g\in\mathcal{G}(b)} Q^g_{g,t}
- \sum_{s\in\mathcal{S}(b)} Q^s_{s,t}
- Q^d_{b,t}
- \sum_{\ell:\,f(\ell)=b} Q^f_{\ell,t}
- \sum_{\ell:\,t(\ell)=b} Q^t_{\ell,t}
+ b^{sh}_b\,V_{b,t}^2 = 0,
\qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{21}
```

### Branch thermal limits

```math
(P^f_{\ell,t})^2 + (Q^f_{\ell,t})^2 \le \overline{S}_\ell^2,
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{22}
```

```math
(P^t_{\ell,t})^2 + (Q^t_{\ell,t})^2 \le \overline{S}_\ell^2,
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{23}
```

### Angular aperture limits

Imposed on the **raw** difference ``\theta_f - \theta_t`` (without ``\Phi_\ell``):

```math
\delta_\ell^{\min} \le \theta_{f(\ell),t} - \theta_{t(\ell),t} \le \delta_\ell^{\max},
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{24}
```

### Reference bus

```math
\theta_{ref,t} = 0, \qquad t\in\mathcal{T}. \tag{25}
```

### Generation constraints

Time envelope; added only when the period limit ``\overline{P}^{ts}_{g,t}`` is below
``P_g^{\max}``.

```math
p^g_{g,t} \le \overline{P}^{ts}_{g,t},
\qquad g\in\mathcal{G}\setminus\mathcal{G}^c,\ t\in\mathcal{T}. \tag{26}
```

```math
p^g_{g,t} \ge \underline{P}^{ts}_{g,t},
\qquad g\in\mathcal{G},\ t\in\mathcal{T}. \tag{27}
```

```math
\underline{Q}^{ts}_{g,t} \le Q^g_{g,t} \le \overline{Q}^{ts}_{g,t},
\qquad g\in\mathcal{G},\ t\in\mathcal{T}. \tag{28}
```

**Renewable generation curtailment.**

```math
p^g_{g,t} + p^c_{g,t} = \overline{P}^{ts}_{g,t},
\qquad g\in\mathcal{G}^c,\ t\in\mathcal{T}. \tag{29}
```

### Storage constraints

**Energy balance.**

```math
E_{s,1} = E^0_s + \eta^c_s\,c_{s,1} - \frac{d_{s,1}}{\eta^d_s},
\qquad s\in\mathcal{S}. \tag{30}
```

```math
E_{s,t} = E_{s,t-1} + \eta^c_s\,c_{s,t} - \frac{d_{s,t}}{\eta^d_s},
\qquad s\in\mathcal{S},\ t = 2,\dots,T. \tag{31}
```

**Non-simultaneity of charging and discharging.**

```math
c_{s,t}\,d_{s,t} = 0, \qquad s\in\mathcal{S},\ t\in\mathcal{T}. \tag{32}
```

**Inverter thermal limit.**

```math
(P^s_{s,t})^2 + (Q^s_{s,t})^2 \le \overline{S}_s^2,
\qquad s\in\mathcal{S},\ t\in\mathcal{T}. \tag{33}
```

**Converter losses (AC coupling).** The ohmic loss is proportional to the square of the
current, ``|S|^2/V^2``, evaluated at the connection bus ``B(s)``:

```math
P^s_{s,t} + d_{s,t} - c_{s,t}
= P^{loss}_s + r_s\,\frac{(P^s_{s,t})^2 + (Q^s_{s,t})^2}{V_{B(s),t}^2},
\qquad s\in\mathcal{S},\ t\in\mathcal{T}. \tag{34}
```

```math
Q^s_{s,t} = Q^{sc}_{s,t} + Q^{loss}_s
+ x_s\,\frac{(P^s_{s,t})^2 + (Q^s_{s,t})^2}{V_{B(s),t}^2},
\qquad s\in\mathcal{S},\ t\in\mathcal{T}. \tag{35}
```

### Complete model

```math
\min_{\substack{\theta,\,V,\,p^g,\,Q^g,\,p^c,\\ E,\,c,\,d,\,P^s,\,Q^s,\,Q^{sc}}}
\quad F + F^{curt} \qquad (\text{Eq. } 14)
```

subject to:

| Block | Equations |
|---|---|
| Flows | (15)–(19) |
| Balance | (20), (21) |
| Network | (22)–(25) |
| Generation | (26)–(29) |
| Storage | (30)–(35) |
| Bounds | (1)–(11) |