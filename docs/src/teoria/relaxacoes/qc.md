# QC (Quadratic Convex)

*Quadratic-Convex relaxation, recursive McCormick (QCRM).*

The QC model is a convex relaxation of the AC power flow in **W-space**, built on the
same four flow expressions as [SOCWR (Jabr conic relaxation)](@ref) (which in turn reuses the rectangular
[AC Rectangular (ACR)](@ref) flows). The difference is where the convexity comes from.
SOCWR treats the cross-blocks ``W_r``, ``W_i`` as opaque variables tied *from the outside*
by a single rotated cone. The QC relaxation **attacks from the inside**. It opens the
internal structure of each block and relaxes every piece separately: the magnitude product,
the cosine, the sine, and the products that rebuild the blocks. In this variant QC does not
carry the Jabr voltage-product cone at all. It trades that one cone for the inside family, so
the two relaxations are different rather than nested. Its full name is **QCRM**,
Quadratic-Convex Recursive-McCormick, *recursive* because the blocks are nested products (a
product of a product), each wrapped in its own McCormick envelope [Coffrin2018](@cite).

!!! info "Relaxation or approximation?"
    QC is a **convex relaxation**: it enlarges the feasible set, so its optimum is a
    guaranteed **lower bound** on the true AC optimum. It is a *different* relaxation from
    SOCWR, not a strict tightening of it. QC drops the Jabr voltage-product cone and fences
    each piece of the blocks on its own, so the two feasible sets are not nested in general.
    See [Foundations](@ref).

## Derivation

### Intermediate variables

Recall what the real block ``W_r`` actually is:

```math
W_r = \underbrace{v_i\,v_j}_{\text{product of magnitudes}} \cdot
      \underbrace{\cos\theta_{ij}}_{\text{cosine of the angle}}.
```

SOCWR never looked inside this expression. The QC relaxation decomposes it: it creates a
variable for each piece and fences each piece with its own convex relaxation. To do so it
introduces the following intermediate variables (the identifiers are the ones used in the
code):

| Variable | Represents |
|---|---|
| `td`  | ``\theta_i - \theta_j`` |
| `vv`  | ``v_i \cdot v_j`` |
| `cs`  | envelope of ``\cos\theta_{ij}`` |
| `si`  | envelope of ``\sin\theta_{ij}`` |
| `ccm` | squared branch current |

The bus variables are magnitude, angle, and the lifted squares ``W_{ii}`` and ``W_{jj}``
(the from- and to-bus squared magnitudes). The branch variables are ``W_r`` and ``W_i``. The idea is to **rebuild the blocks as a chain of
products**:

```math
W_r = \mathrm{vv}\cdot\mathrm{cs}, \qquad
W_i = \mathrm{vv}\cdot\mathrm{si}, \qquad
\mathrm{vv} = v_i\cdot v_j.
```

### Branch flows

The QC flows are **exactly** those of SOCWR, which in turn were those of the ACR model.
With ``tm = |T|^2 = \tau^2`` the squared tap magnitude and ``t_r``, ``t_i`` the real and
imaginary parts of the complex tap ``T``:

```math
p_{ij} = \frac{(g + g^{fr})\,W_{ii}}{tm}
       + \frac{(-g\,t_r + b\,t_i)\,W_r}{tm}
       + \frac{(-b\,t_r - g\,t_i)\,W_i}{tm},
```

```math
q_{ij} = -\frac{(b + b^{fr})}{tm}\,W_{ii}
       - \frac{(-b\,t_r - g\,t_i)}{tm}\,W_r
       + \frac{(-g\,t_r + b\,t_i)}{tm}\,W_i,
```

```math
p_{ji} = (g + g^{to})\,W_{jj}
       + \frac{(-g\,t_r - b\,t_i)}{tm}\,W_r
       + \frac{(b\,t_r - g\,t_i)}{tm}\,W_i,
```

```math
q_{ji} = -(b + b^{to})\,W_{jj}
       - \frac{(-b\,t_r + g\,t_i)}{tm}\,W_r
       + \frac{(g\,t_r + b\,t_i)}{tm}\,W_i.
```

### Bounds on the intermediate variables

Before anything is relaxed, each new variable needs bounds, because McCormick requires
bounds on the factors.

**`td`.** Directly from the branch angle limits:

```math
\theta^{\min} \le \mathrm{td} \le \theta^{\max}.
```

**`vv`.** A product of two positive factors, so minimum with minimums, maximum with
maximums:

```math
v_i^{\min}\,v_j^{\min} \le \mathrm{vv} \le v_i^{\max}\,v_j^{\max}.
```

**`si`.** The sine is increasing on ``[-\pi/2, \pi/2]`` (its derivative, the cosine, is
positive there); an increasing function maps the endpoints of the interval onto the
endpoints of the range (a single case):

```math
\sin\theta^{\min} \le \mathrm{si} \le \sin\theta^{\max}.
```

**`cs`.** The cosine is **not** monotone on intervals that cross zero (it rises on the
left, falls on the right, peaks at zero), so its bounds depend on where the interval lies:

```math
\begin{aligned}
\theta^{\min} \ge 0 &:\quad \cos\theta^{\max} \le \mathrm{cs} \le \cos\theta^{\min},\\
\theta^{\max} \le 0 &:\quad \cos\theta^{\min} \le \mathrm{cs} \le \cos\theta^{\max},\\
\text{crosses zero} &:\quad \min(\cos\theta^{\min}, \cos\theta^{\max}) \le \mathrm{cs} \le 1.
\end{aligned}
```

**`wr`, `wi`.** The magnitudes ride along, ``W_r = v_i v_j \cos\theta_{ij}`` and
``W_i = v_i v_j \sin\theta_{ij}`` are products ``a\cdot s`` with ``a = v_i v_j \ge 0``. The
cosine is ``\ge 0`` on ``[-\pi/2, \pi/2]``; the sine changes sign between cases:

```math
\begin{aligned}
\theta^{\min} \ge 0 &:\quad
  v_i^{\min} v_j^{\min}\cos\theta^{\max} \le W_r \le v_i^{\max} v_j^{\max}\cos\theta^{\min},\\
 &\phantom{:}\quad
  v_i^{\min} v_j^{\min}\sin\theta^{\min} \le W_i \le v_i^{\max} v_j^{\max}\sin\theta^{\max};\\[4pt]
\theta^{\max} \le 0 &:\quad
  v_i^{\min} v_j^{\min}\cos\theta^{\min} \le W_r \le v_i^{\max} v_j^{\max}\cos\theta^{\max},\\
 &\phantom{:}\quad
  v_i^{\max} v_j^{\max}\sin\theta^{\min} \le W_i \le v_i^{\min} v_j^{\min}\sin\theta^{\max};\\[4pt]
\text{crosses zero} &:\quad
  v_i^{\min} v_j^{\min}\min(\cos\theta^{\min}, \cos\theta^{\max}) \le W_r \le v_i^{\max} v_j^{\max},\\
 &\phantom{:}\quad
  v_i^{\max} v_j^{\max}\sin\theta^{\min} \le W_i \le v_i^{\max} v_j^{\max}\sin\theta^{\max}.
\end{aligned}
```

``\mathrm{ccm}`` is deliberately left out of this section: its bounds only make sense once
it is defined (Constraint 8).

## Constraints

Each constraint is numbered continuously within the model.

### Constraint 1: secant of the magnitude

SOCWR and QC promoted ``W_{ii}`` to an independent variable: it is no longer computed from
``v_i`` but chosen by the optimizer. A free ``W_{ii}`` can *lie*. The optimizer may pick a
value that does not correspond to any true ``v_i^2``. We therefore need constraints that
pin ``W_{ii}`` near ``v_i^2``. The true relation is not convex as an equality:

```math
W_{ii} = v_i^2
```

is a parabola, and equality with a curve is non-convex. Instead of forcing ``W_{ii}`` onto
the parabola, we trap it between a floor and a ceiling.

**Floor.** ``W_{ii} \ge v_i^2``, easy and convex (the region above a parabola is the
interior of a bowl, which is convex).

**Ceiling.** ``W_{ii} \le v_i^2`` would mark the region *below* the parabola, which is
non-convex (the exterior of the bowl). So the parabola cannot be used as a ceiling. A
**line** can: below a line is a half-plane, always convex. The line must (i) stay above the
parabola on all of ``[v_i^{\min}, v_i^{\max}]`` (otherwise it would cut valid physical
points) and (ii) be as low as possible (the tighter the band, the better the relaxation).
The **secant** through the two endpoints of the parabola satisfies both. Its endpoints are
``P_1 = (v_i^{\min}, (v_i^{\min})^2)`` and ``P_2 = (v_i^{\max}, (v_i^{\max})^2)``, so the
slope is

```math
m = \frac{(v_i^{\max})^2 - (v_i^{\min})^2}{v_i^{\max} - v_i^{\min}} = v_i^{\max} + v_i^{\min},
```

and the secant line reads

```math
W_{ii} = (v_i^{\min} + v_i^{\max})\,v_i - v_i^{\min}\,v_i^{\max}.
```

The complete band:

```math
v_i^2 \le W_{ii} \le (v_i^{\min} + v_i^{\max})\,v_i - v_i^{\min}\,v_i^{\max}.
```

### Constraint 2: exact angle coupling

```math
\theta_i - \theta_j = \mathrm{td}.
```

An exact, linear equality. This is the piece SOCWR did not have: there the blocks lived
loose, with no link to the bus angles. Here the bus angle ``\theta`` ties the branch
difference ``\mathrm{td}``.

!!! note "This restores angle cycle consistency"
    Because ``\mathrm{td}`` comes from ``\theta_i - \theta_j`` and all envelopes rest on
    ``\mathrm{td}``, the angle differences telescope around any loop: ``\sum \mathrm{td} = 0``
    on a cycle, automatically, since each ``\theta_i`` is a single number. This returns to QC
    the angle cycle consistency that SOCWR lost. The blocks ``W_r, W_i`` still reach those
    angles only through the relaxed envelopes, so phasor consistency itself stays relaxed.

### Constraint 3: cosine envelope

We want ``\mathrm{cs} = \cos(\mathrm{td})``, but equality with a curve is non-convex. We
fence ``\mathrm{cs}`` between a convex floor and a convex ceiling. Since the cosine is
concave on ``[-\pi/2, \pi/2]``, concavity fixes the roles: the ceiling is a curve above,
the floor is a secant below. Define ``m = \max(|\theta^{\min}|, |\theta^{\max}|)``.

**Ceiling (the parabola).** Take an inverted parabola with its peak where the cosine peaks.
The cosine peaks at ``\mathrm{td} = 0`` with value ``\cos 0 = 1``, so the parabola must have
peak ``1`` at ``\mathrm{td} = 0``. From the general form ``C_s(\mathrm{td}) = a\,\mathrm{td}^2 + b\,\mathrm{td} + c``,
imposing ``C_s(0) = 1`` gives ``c = 1``, and a zero slope at the peak
(``C_s'(\mathrm{td}) = 2a\,\mathrm{td} + b``, ``C_s'(0) = b = 0``) leaves

```math
C_s = 1 - k\,\mathrm{td}^2, \qquad k > 0.
```

The coefficient ``k`` is fixed by requiring the parabola to touch the cosine at the ends
``\mathrm{td} = \pm m``: as low as it can go while still touching both tips (any larger
``k`` would cut the cosine; any smaller ``k`` would leave it higher and looser). At
``\mathrm{td} = m`` the cosine is ``\cos m``, so

```math
1 - k\,m^2 = \cos m \;\Longrightarrow\; k = \frac{1 - \cos m}{m^2},
```

and the ceiling is

```math
\mathrm{cs} \le 1 - \frac{1 - \cos m}{m^2}\,\mathrm{td}^2.
```

The expression for ``k`` is used only for ``m > 0``. When the aperture collapses to a point
(``m = 0``, no angle freedom on the branch) the envelope degenerates and ``\mathrm{cs} = 1``
is fixed.

The two curves coincide at three points, ``\mathrm{td} \in \{-m, 0, +m\}``. That they touch
at the tips does not by itself guarantee the parabola stays above in between. The proof is
below.

!!! note "Proof that the parabola envelops the cosine"
    Define the gap

    ```math
    f(\mathrm{td}) = \left(1 - \frac{1 - \cos m}{m^2}\,\mathrm{td}^2\right) - \cos(\mathrm{td}),
    ```

    and show ``f(\mathrm{td}) \ge 0`` for all ``\mathrm{td} \in [-m, m]``. By construction
    ``f`` vanishes at three points: ``f(0) = 0`` (both equal ``1``), ``f(m) = 0`` (both
    equal ``\cos m``) and ``f(-m) = 0`` (by symmetry, both functions are even). The question
    is the sign of ``f`` between these zeros. Differentiate twice:

    ```math
    f'(\mathrm{td}) = -2k\,\mathrm{td} + \sin(\mathrm{td}), \qquad
    f''(\mathrm{td}) = \cos(\mathrm{td}) - 2k.
    ```

    The term ``\cos(\mathrm{td})`` ranges between ``\cos m`` (at the tips) and ``1`` (at the
    centre), while ``2k = 2(1-\cos m)/m^2`` is a constant. We locate ``2k`` by proving, for
    ``m \in (0, \pi/2]``,

    ```math
    \cos m < \frac{2(1 - \cos m)}{m^2} < 1.
    ```

    *Part 1 (``2k < 1``).* This is ``1 - \cos m < m^2/2``. From the geometric fact
    ``\sin x < x`` for ``x > 0`` (the vertical segment is shorter than the arc on the unit
    circle), integrate from ``0`` to ``m``:

    ```math
    \int_0^m \sin x\,dx < \int_0^m x\,dx \;\Longrightarrow\; 1 - \cos m < \frac{m^2}{2}.
    ```

    *Part 2 (``2k > \cos m``).* This is ``2(1 - \cos m) > m^2 \cos m``. Define
    ``g(m) = 2(1 - \cos m) - m^2 \cos m`` with ``g(0) = 0`` and differentiate:

    ```math
    g'(m) = 2\sin m - \big(2m\cos m - m^2 \sin m\big)
          = \underbrace{2(\sin m - m\cos m)}_{A} + \underbrace{m^2 \sin m}_{B}.
    ```

    ``B > 0`` for ``0 < m \le \pi/2``. For ``A``, set ``h(m) = \sin m - m\cos m`` with
    ``h(0) = 0``; then ``h'(m) = \cos m - (\cos m - m\sin m) = m\sin m > 0`` on
    ``0 < m < \pi/2``, so ``h`` starts at zero and grows, giving ``A > 0``. With ``A > 0``
    and ``B > 0``, ``g' > 0``: ``g`` starts at zero and grows, so ``g(m) > 0``.

    *Conclusion.* With ``\cos m < 2k < 1`` located, the sign of
    ``f'' = \cos(\mathrm{td}) - 2k`` changes once inside the interval: near the centre
    (``\cos\mathrm{td}`` close to ``1 > 2k``) we have ``f'' > 0``, and ``f`` leaves the zero
    at ``\mathrm{td} = 0`` curving upward; near the tips (``\cos\mathrm{td}`` close to
    ``\cos m < 2k``) we have ``f'' < 0``, and ``f`` comes back down to zero at ``\pm m`` from
    above. It rises from the central zero, falls to the tip zeros, and never crosses:
    ``f \ge 0`` on all of ``[-m, m]``. The parabola envelops the cosine. ∎

**Floor (the secant).** The line through the two endpoints of the cosine on the interval,
``P_1 = (\theta^{\min}, \cos\theta^{\min})`` and ``P_2 = (\theta^{\max}, \cos\theta^{\max})``,
with slope ``m_{\text{sec}} = (\cos\theta^{\max} - \cos\theta^{\min})/(\theta^{\max} - \theta^{\min})``:

```math
\mathrm{cs} \ge \frac{\cos\theta^{\max} - \cos\theta^{\min}}{\theta^{\max} - \theta^{\min}}
  \,(\mathrm{td} - \theta^{\min}) + \cos\theta^{\min}.
```

The complete frame: secant ``\le \mathrm{cs} \le`` inverted parabola.

### Constraint 4: sine envelope

The sine is nearly linear near zero and, unlike the cosine, changes concavity in the middle
of the window: it is **concave** for ``\theta > 0``, **convex** for ``\theta < 0``, with an
**inflection** at ``\theta = 0``. Because the envelope technique (a curve on one side, a
secant on the other) needs a *fixed* concavity, the sine splits into three cases. The tool
is the tangent line: since ``\frac{d}{d\theta}\sin\theta = \cos\theta``, the tangent at
``(\alpha, \sin\alpha)`` is

```math
\mathrm{si} = \cos(\alpha)\,(\mathrm{td} - \alpha) + \sin\alpha.
```

A concave function lies **below** every one of its tangents; a convex one lies **above**.
With ``m = \max(|\theta^{\min}|, |\theta^{\max}|)``:

**Case 1 (interval all positive, ``\theta^{\min} \ge 0``).** The sine is concave
throughout; the envelope matches the cosine's. Floor = secant (concave ⇒ chord below the
curve); ceiling = tangent at ``\alpha = m/2`` (concave ⇒ tangent above the curve):

```math
\mathrm{si} \le \cos(m/2)\,(\mathrm{td} - m/2) + \sin(m/2).
```

The arc centre ``m/2`` minimizes the worst gap, but the tangent is valid at any point on the
concave side.

**Case 2 (interval all negative, ``\theta^{\max} \le 0``).** The sine is convex throughout,
so ceiling and floor swap. Ceiling = secant (convex ⇒ chord above); floor = tangent at
``\alpha = -m/2`` (convex ⇒ tangent below):

```math
\mathrm{si} \le \frac{\sin\theta^{\min} - \sin\theta^{\max}}{\theta^{\min} - \theta^{\max}}
  \,(\mathrm{td} - \theta^{\min}) + \sin\theta^{\min},
```

```math
\mathrm{si} \ge \cos(m/2)\,(\mathrm{td} + m/2) - \sin(m/2).
```

**Case 3 (interval crosses zero, ``\theta^{\min} < 0 < \theta^{\max}``).** The interval
contains the inflection, so no secant is guaranteed on either side; use **two tangents**,
one per side:

```math
\mathrm{si} \le \cos(m/2)\,(\mathrm{td} - m/2) + \sin(m/2) \quad \text{(tangent at } +m/2,\ \text{ceiling)},
```

```math
\mathrm{si} \ge \cos(m/2)\,(\mathrm{td} + m/2) - \sin(m/2) \quad \text{(tangent at } -m/2,\ \text{floor)}.
```

The tangent on the concave side binds from above, the one on the convex side from below,
and each stays valid across the whole interval.

!!! note "Proof of the tangent bound"
    Take the ceiling ``T(\mathrm{td}) = \cos(\alpha)\,(\mathrm{td} - \alpha) + \sin\alpha``
    with ``\alpha = m/2``.

    **Step 1 (gap function).**

    ```math
    \varphi(\mathrm{td}) = T(\mathrm{td}) - \sin(\mathrm{td}).
    ```

    We want ``\varphi \ge 0`` on all of ``[-m, m]``. By construction ``\varphi(\alpha) = 0``
    and ``\varphi'(\alpha) = 0``.

    **Step 2 (derivative of the gap).**

    ```math
    \varphi'(\mathrm{td}) = \cos(\alpha) - \cos(\mathrm{td}).
    ```

    This is a comparison of cosines; since the cosine is even and decreasing on
    ``[0, \pi/2]``, ``\cos(\mathrm{td}) = \cos(|\mathrm{td}|)`` and the comparison depends
    only on the distance to zero:

    ```math
    |\mathrm{td}| < \alpha \;\to\; \cos(\mathrm{td}) > \cos(\alpha) \;\to\; \varphi' < 0,
    \qquad
    |\mathrm{td}| > \alpha \;\to\; \varphi' > 0.
    ```

    **Step 3 (the gap map).** The sign forces the shape of ``\varphi`` on ``[-m, m]``:
    rising on ``[-m, -\alpha]``, falling on ``[-\alpha, \alpha]``, rising on ``[\alpha, m]``,
    an *up, down, up* profile. A function of this shape attains its minimum only at the end of
    the descent (``\mathrm{td} = \alpha``, where ``\varphi = 0``) or at the left tip
    ``\mathrm{td} = -m``; everything else is above these two:

    ```math
    \min_{[-m,\,m]} \varphi = \min\{\varphi(-m),\, 0\}.
    ```

    **Step 4 (the left tip is non-negative).** Only ``\varphi(-m) \ge 0`` remains. Write
    ``u = m/2`` and use ``\sin m = 2\sin u\cos u`` to factor out the positive cosine:

    ```math
    \varphi(-m) = \sin u + \sin 2u - 3u\cos u = \cos u\,\big(\tan u + 2\sin u - 3u\big).
    ```

    Let ``h(u) = \tan u + 2\sin u - 3u``. Then ``h(0) = 0`` and
    ``h'(u) = \sec^2 u + 2\cos u - 3`` with ``h'(0) = 0``; differentiating once more,

    ```math
    h''(u) = 2\sec^2 u\,\tan u - 2\sin u = 2\sin u\,\frac{1 - \cos^3 u}{\cos^3 u} > 0
    \quad\text{on } \left(0, \tfrac{\pi}{4}\right].
    ```

    So ``h'`` rises from ``0`` and stays positive, ``h`` rises from ``0`` and stays positive,
    and ``\varphi(-m) = \cos u\,h(u) \ge 0``. Hence ``\varphi \ge 0`` on all of ``[-m, m]``:
    the tangent at ``m/2`` stays above the sine across the window. The all-negative case is
    its mirror image, the floor tangent at ``-m/2`` bounding the sine from below. ∎

### Constraint 5: McCormick products (the full envelope)

For a product ``z = x\cdot y`` with ``x \in [x^{\min}, x^{\max}]`` and
``y \in [y^{\min}, y^{\max}]``, the complete envelope is four linear inequalities, each born
from a product of two factors of known sign:

```math
\begin{aligned}
\text{Floor 1 (both at minimum)}: &\quad z \ge x^{\min} y + y^{\min} x - x^{\min} y^{\min},\\
\text{Floor 2 (both at maximum)}: &\quad z \ge x^{\max} y + y^{\max} x - x^{\max} y^{\max},\\
\text{Ceiling 1 (opposite sides)}: &\quad z \le x^{\min} y + y^{\max} x - x^{\min} y^{\max},\\
\text{Ceiling 2 (opposite sides)}: &\quad z \le x^{\max} y + y^{\min} x - x^{\max} y^{\min}.
\end{aligned}
```

Together they are the convex hull of the graph of ``z = xy`` over the box
``[x^{\min},x^{\max}]\times[y^{\min},y^{\max}]``, the tightest convex envelope of a single
bilinear term. QC instantiates this envelope three times, the **recursive** part, because
the output of one product feeds the bounds of the next.

**Product 1, ``\mathrm{vv} = v_i\cdot v_j``** (``z = \mathrm{vv}``, ``x = v_i``, ``y = v_j``):

```math
\begin{aligned}
\mathrm{vv} &\ge v_i^{\min} v_j + v_j^{\min} v_i - v_i^{\min} v_j^{\min},\\
\mathrm{vv} &\ge v_i^{\max} v_j + v_j^{\max} v_i - v_i^{\max} v_j^{\max},\\
\mathrm{vv} &\le v_i^{\min} v_j + v_j^{\max} v_i - v_i^{\min} v_j^{\max},\\
\mathrm{vv} &\le v_i^{\max} v_j + v_j^{\min} v_i - v_i^{\max} v_j^{\min}.
\end{aligned}
```

**Product 2, ``W_r = \mathrm{vv}\cdot\mathrm{cs}``** (``z = W_r``, ``x = \mathrm{vv}``,
``y = \mathrm{cs}``), with ``\mathrm{vv} \in [\mathrm{vv}^{\min}, \mathrm{vv}^{\max}]``
(``\mathrm{vv}^{\min} = v_i^{\min} v_j^{\min}``, ``\mathrm{vv}^{\max} = v_i^{\max} v_j^{\max}``)
and ``\mathrm{cs} \in [\mathrm{cs}^{\min}, \mathrm{cs}^{\max}]`` (the cosine-envelope bounds,
the three cases above):

```math
\begin{aligned}
W_r &\ge \mathrm{vv}^{\min}\,\mathrm{cs} + \mathrm{cs}^{\min}\,\mathrm{vv} - \mathrm{vv}^{\min}\,\mathrm{cs}^{\min},\\
W_r &\ge \mathrm{vv}^{\max}\,\mathrm{cs} + \mathrm{cs}^{\max}\,\mathrm{vv} - \mathrm{vv}^{\max}\,\mathrm{cs}^{\max},\\
W_r &\le \mathrm{vv}^{\min}\,\mathrm{cs} + \mathrm{cs}^{\max}\,\mathrm{vv} - \mathrm{vv}^{\min}\,\mathrm{cs}^{\max},\\
W_r &\le \mathrm{vv}^{\max}\,\mathrm{cs} + \mathrm{cs}^{\min}\,\mathrm{vv} - \mathrm{vv}^{\max}\,\mathrm{cs}^{\min}.
\end{aligned}
```

The ``\mathrm{vv}`` that enters here is itself already tied by the four cuts of Product 1.

**Product 3, ``W_i = \mathrm{vv}\cdot\mathrm{si}``** (``z = W_i``, ``x = \mathrm{vv}``,
``y = \mathrm{si}``), with ``\mathrm{si} \in [\sin\theta^{\min}, \sin\theta^{\max}]``:

```math
\begin{aligned}
W_i &\ge \mathrm{vv}^{\min}\,\mathrm{si} + \mathrm{si}^{\min}\,\mathrm{vv} - \mathrm{vv}^{\min}\,\mathrm{si}^{\min},\\
W_i &\ge \mathrm{vv}^{\max}\,\mathrm{si} + \mathrm{si}^{\max}\,\mathrm{vv} - \mathrm{vv}^{\max}\,\mathrm{si}^{\max},\\
W_i &\le \mathrm{vv}^{\min}\,\mathrm{si} + \mathrm{si}^{\max}\,\mathrm{vv} - \mathrm{vv}^{\min}\,\mathrm{si}^{\max},\\
W_i &\le \mathrm{vv}^{\max}\,\mathrm{si} + \mathrm{si}^{\min}\,\mathrm{vv} - \mathrm{vv}^{\max}\,\mathrm{si}^{\min}.
\end{aligned}
```

This is the full *recursive-McCormick* chain: nested products, each with its envelope of
four, the ``\mathrm{cs}`` and ``\mathrm{si}`` envelopes feeding the bounds.

### Constraint 6: thermal limit

```math
p_{ij}^2 + q_{ij}^2 \le \mathrm{rate}_a^2, \qquad p_{ji}^2 + q_{ji}^2 \le \mathrm{rate}_a^2.
```

### Constraint 7: angle limits via tangent

```math
\tan(\theta^{\min})\,W_r \le W_i \le \tan(\theta^{\max})\,W_r.
```

The angle limit becomes a linear cut on the blocks. It is valid because ``W_r > 0`` in the
physical region (``v_i v_j > 0`` and ``\cos\theta_{ij} > 0`` for ``|\theta_{ij}| < \pi/2``).

### Constraint 8: the squared current (`ccm`)

**Motivation.** The optimizer might find blocks ``W_r, W_i, W_{ii}, W_{jj}`` and flows
``p, q`` that satisfy the thermal limit yet correspond to no coherent current. The physical
current–voltage–power relation is one more tie the relaxation can exploit; QC injects it
through a new variable.

**The branch model.** Recall the branch (see [AC Polar](@ref)): an ideal transformer on the
*from* end, an internal node at voltage ``V' = V_i/T``, a series element ``y`` and a local
shunt ``y^{fr} = g_{fr} + j\,b_{fr}``. The total current entering the internal node is

```math
I' = \underbrace{y\,(V' - V_j)}_{I_{\text{series}}} + \underbrace{y^{fr}\,V'}_{I_{\text{shunt}}}.
```

Define

```math
\mathrm{ccm} := |I'|^2
```

(the squared current at the internal node; since ``I_i = I'/T^*``, it equals
``\tau^2 |I_i|^2``, with ``\tau = |T|``).

**Step 1 (the bridge with the flows).** By ideal-transformer conservation,
``S_{ij} = V_i I_i^* = V' I'^*``; taking magnitudes, ``|S_{ij}|^2 = |V'|^2 |I'|^2``. With
``|S_{ij}|^2 = p_{ij}^2 + q_{ij}^2`` and ``|V'|^2 = |V_i|^2/\tau^2 = W_{ii}/tm``, this is
physically the equality ``p_{ij}^2 + q_{ij}^2 = (W_{ii}/tm)\,\mathrm{ccm}``. Imposing it
exactly would bring back non-convexity (a product of variables on one side), so we relax it:

```math
p_{ij}^2 + q_{ij}^2 \le \frac{W_{ii}}{tm}\,\mathrm{ccm},
```

a rotated second-order cone, the same geometry as the Jabr tie in [SOCWR (Jabr conic relaxation)](@ref).

**Step 2, expressing `ccm` in the blocks (series current).** On its own, Step 1 says
nothing (``\mathrm{ccm}`` could inflate freely); the force comes from the algebraic equation
tying it to the blocks. Start with the series current:

```math
|I_{\text{series}}|^2 = |y|^2\,|V' - V_j|^2
  = (g^2 + b^2)\Big(|V'|^2 + |V_j|^2 - 2\,\mathrm{Re}(V' V_j^*)\Big).
```

Open each piece in the blocks: ``|V'|^2 = W_{ii}/tm``, ``|V_j|^2 = W_{jj}``, and, using
``V_i V_j^* = W_r + j\,W_i`` and ``1/T = T^*/|T|^2 = (t_r - j\,t_i)/tm``,

```math
\mathrm{Re}(V' V_j^*) = \mathrm{Re}\!\left(\frac{V_i V_j^*}{T}\right) = \frac{t_r\,W_r + t_i\,W_i}{tm}.
```

Collecting,

```math
|I_{\text{series}}|^2 = (g^2 + b^2)\left(\frac{W_{ii}}{tm} + W_{jj} - \frac{2\,(t_r\,W_r + t_i\,W_i)}{tm}\right).
```

Linear in the blocks, and this is why the computation works: the squared series current, which
looked quadratic, is linear in ``W_{ii}, W_{jj}, W_r, W_i``.

**Step 3 (the shunt correction).** ``\mathrm{ccm}`` is ``|I'|^2``, not ``|I_{\text{series}}|^2``.
From ``S_{ij} = V'(I')^*`` we get ``I' = S_{ij}^*/(V')^*``, and
``I_{\text{series}} = I' - y^{fr} V'``. Expanding
``|A - B|^2 = |A|^2 - 2\,\mathrm{Re}(AB^*) + |B|^2`` with ``A = S_{ij}^*/(V')^*`` and
``B = y^{fr} V'``:

- ``|A|^2 = |S_{ij}|^2/|V'|^2 = \mathrm{ccm}``,
- ``|B|^2 = |y^{fr}|^2 |V'|^2 = (g_{fr}^2 + b_{fr}^2)\,W_{ii}/tm``,
- ``AB^* = S_{ij}^*\,(y^{fr})^* = (S_{ij}\,y^{fr})^*``, so ``\mathrm{Re}(AB^*) = \mathrm{Re}(S_{ij}\,y^{fr}) = g_{fr}\,p_{ij} - b_{fr}\,q_{ij}``.

Hence

```math
|I_{\text{series}}|^2 = \mathrm{ccm} - 2\,(g_{fr}\,p_{ij} - b_{fr}\,q_{ij}) + (g_{fr}^2 + b_{fr}^2)\,\frac{W_{ii}}{tm}.
```

Isolating ``\mathrm{ccm}`` and substituting Step 2 gives the closed form:

```math
\mathrm{ccm} = (g^2 + b^2)\left(\frac{W_{ii}}{tm} + W_{jj} - \frac{2\,(t_r\,W_r + t_i\,W_i)}{tm}\right)
  - (g_{fr}^2 + b_{fr}^2)\,\frac{W_{ii}}{tm} + 2\,(g_{fr}\,p_{ij} - b_{fr}\,q_{ij}).
```

**Ceiling of `ccm`.** From ``|S_{ij}| \le \mathrm{rate}_a`` and ``|S_{ij}| = |V'|\,|I'|``:

```math
|I'| = \frac{|S_{ij}|}{|V'|} \le \frac{\mathrm{rate}_a}{v_i^{\min}/\tau} = \frac{\mathrm{rate}_a\,\tau}{v_i^{\min}}
  \;\Longrightarrow\;
\mathrm{ccm} \le \left(\frac{\mathrm{rate}_a\,\tau}{v_i^{\min}}\right)^2.
```

### Constraint 9: Lifted Nonlinear Cuts (LNC)

The QC model in Lirion keeps the **same two LNC cuts as SOCWR**. Their defining quantities are

```math
\phi = \frac{\theta^{\max} + \theta^{\min}}{2}, \qquad
d = \frac{\theta^{\max} - \theta^{\min}}{2},
```

```math
s_f = v_f^{\min} + v_f^{\max}, \qquad s_t = v_t^{\min} + v_t^{\max}.
```

!!! note "The two LNC cuts come from the SOCWR page"
    QC posts the same pair of lifted nonlinear cuts as
    [SOCWR (Jabr conic relaxation)](@ref), on the shared ``W_r, W_i`` blocks. They are
    derived in full there (the McCormick bridge and the explicit Cut 1 / Cut 2), so this page
    cross-references them rather than repeating the construction.

## QC Formulation

The formulation below collects, in consolidated form and **without re-deriving**, exactly
the blocks the model builds. The derivation above justifies each piece. The model is assembled
as an **NLP** and solved with Ipopt.

### Sets

- ``\mathcal{N}``: set of buses.
- ``\mathcal{N}^{nv} \subseteq \mathcal{N}``: non-virtual buses, on which the power balance
  is enforced (`non_virtual_bus_ids`).
- ``\mathcal{L}``: set of active branches (in-service lines and transformers).
- ``\mathcal{P}``: set of ordered ``(f,t)`` **bus pairs** carrying at least one active
  branch. The cross-blocks ``W_r, W_i`` live on this set, one product per pair, shared by
  all parallel branches (`add_rt_variables_with_bounds!`).
- ``\mathcal{G}``: set of existing generators; ``\mathcal{G}^c \subseteq \mathcal{G}``:
  curtailable (renewable) generators.
- ``\mathcal{S}``: set of storage units.
- ``\mathcal{T} = \{1,\dots,T\}``: time periods.

### Indices and incidence functions

- ``b, i, j \in \mathcal{N}``; ``\ell \in \mathcal{L}``; ``p \in \mathcal{P}``;
  ``g \in \mathcal{G}``; ``s \in \mathcal{S}``; ``t \in \mathcal{T}``.
- ``f(\ell),\ t(\ell)``: *from* and *to* buses of branch ``\ell``.
- ``p(\ell) \in \mathcal{P}``: the bus pair of branch ``\ell`` (``\mathrm{pair\_of}`` in the
  code); ``W_r, W_i`` of ``\ell`` are the pair variables of ``p(\ell)``.
- ``B(s)``: connection bus of storage ``s``.

### Parameters

**Branches.**

- ``g_\ell,\ b_\ell``: series conductance and susceptance (`g`, `b`).
- ``g^{fr}_\ell,\ b^{fr}_\ell``; ``g^{to}_\ell,\ b^{to}_\ell``: shunt conductance/susceptance
  of the *from* and *to* ends (`g_fr,b_fr,g_to,b_to`).
- ``\tau_\ell``: tap magnitude (`tap_position`); ``k_\ell = 1/\tau_\ell`` (`inv_tap`);
  ``\Phi_\ell``: phase shift (`tap_angle`).
- ``t_r = \tau_\ell\cos\Phi_\ell`` (`tap_real`), ``t_i = \tau_\ell\sin\Phi_\ell``
  (`tap_imag`); ``tm_\ell = |T_\ell|^2 = \tau_\ell^2`` the squared tap magnitude. Every tap
  denominator on this page is ``1/tm_\ell``, which the code writes as ``\texttt{inv\_tap}^2``
  in the flows and as ``1/\texttt{tap\_position}^2`` in the current tie.
- ``\mathrm{rate}_{a,\ell}``: branch apparent-power (thermal) limit (`rate_a`).
- ``\delta_\ell^{\min},\ \delta_\ell^{\max}``: angular-aperture limits (`ang_min`,
  `ang_max`), written ``\theta^{\min},\theta^{\max}`` in the derivation above.

**Buses.**

- ``V_b^{\min},\ V_b^{\max}``: voltage-magnitude limits.
- ``g^{sh}_b,\ b^{sh}_b``: bus shunt conductance and susceptance.
- ``P^d_{b,t},\ Q^d_{b,t}``: active and reactive demand.

**Generators.** As in [AC Polar](@ref): ``P_g^{\min/\max}``, ``Q_g^{\min/\max}``, the
active/reactive time envelopes ``\overline{P}^{ts}_{g,t}, \underline{P}^{ts}_{g,t}``,
``\overline{Q}^{ts}_{g,t}, \underline{Q}^{ts}_{g,t}``, the cost coefficients ``c_{g,t,k}``
and the unit curtailment cost ``c^{curt}_g``.

**Storage.** As in [AC Polar](@ref): ``E_s^{\min/\max}``, ``E^0_s``, charge/discharge ramps
``C_s^{\max}, D_s^{\max}``, efficiencies ``\eta^c_s, \eta^d_s``, converter ``r_s, x_s``,
fixed losses ``P^{loss}_s, Q^{loss}_s``, inverter limit ``\overline{S}_s``
(`thermal_limit`), and ``Q_s^{\min/\max}``. The ``\mathrm{ccm}^s`` ceiling uses a single
**global** minimum voltage ``V^{\min}`` (`model.ext[:v_min]`, default ``0.9``), not the
per-bus bound.

### Decision variables

**Voltage / W-space.**

```math
\theta_{b,t} \in \mathbb{R} \;\text{: voltage angle at bus } b. \tag{1}
```

```math
V_{b,t} \in [V_b^{\min}, V_b^{\max}] \;\text{: voltage magnitude.} \tag{2}
```

```math
W_{b,t} \in [(V_b^{\min})^2, (V_b^{\max})^2]
  \;\text{: lifted squared magnitude } v_b^2 \ (\texttt{v\_magnitude\_sq}). \tag{3}
```

**Branch / bus-pair blocks.** The cross-blocks are per bus pair ``p``; the rest are per
branch ``\ell``. For ``\ell`` write ``W_r \equiv W^r_{p(\ell),t}``,
``W_i \equiv W^i_{p(\ell),t}``, ``W_{ii} \equiv W_{f(\ell),t}``,
``W_{jj} \equiv W_{t(\ell),t}``.

```math
W^r_{p,t} \in [W^{r,\min}_p, W^{r,\max}_p],\qquad
W^i_{p,t} \in [W^{i,\min}_p, W^{i,\max}_p]
  \;\text{: cross blocks } v_iv_j\cos\theta_{ij},\ v_iv_j\sin\theta_{ij}. \tag{4}
```

```math
\mathrm{td}_{\ell,t} \in [\delta_\ell^{\min}, \delta_\ell^{\max}]
  \;\text{: branch angle difference.} \tag{5}
```

```math
\mathrm{vv}_{\ell,t} \in [V_i^{\min}V_j^{\min}, V_i^{\max}V_j^{\max}]
  \;\text{: magnitude product } v_iv_j. \tag{6}
```

```math
\mathrm{cs}_{\ell,t} \in [\mathrm{cs}^{\min}_\ell, \mathrm{cs}^{\max}_\ell]
  \;\text{: cosine envelope,} \qquad
\mathrm{si}_{\ell,t} \in [\sin\delta_\ell^{\min}, \sin\delta_\ell^{\max}]
  \;\text{: sine envelope,} \tag{7}
```

with the cosine box set by the three-case rule
(``\delta^{\min}\!\ge 0``: ``[\cos\delta^{\max},\cos\delta^{\min}]``;
``\delta^{\max}\!\le 0``: ``[\cos\delta^{\min},\cos\delta^{\max}]``;
crossing zero: ``[\min(\cos\delta^{\min},\cos\delta^{\max}),\,1]``).

```math
\mathrm{ccm}_{\ell,t} \in \Big[0,\ \big(\mathrm{rate}_{a,\ell}\,\tau_\ell / V_i^{\min}\big)^2\Big]
  \;\text{: squared internal-node current.} \tag{8}
```

The ``W_r, W_i`` pair bounds are the PowerModels voltage-product bounds computed per bus
pair on the **tightest** aggregated aperture (``\max`` of ``\delta^{\min}``, ``\min`` of
``\delta^{\max}`` across parallel branches), the per-pair version of the ``W_r, W_i``
bounds derived above (`_rt_pair_product_bounds`). The ``\mathrm{ccm}`` ceiling is applied
only when ``\mathrm{rate}_{a,\ell}`` is finite and positive.

**Generation.**

```math
p^g_{g,t} \in [P_g^{\min}, P_g^{\max}],\qquad
Q^g_{g,t} \in [Q_g^{\min}, Q_g^{\max}],\qquad
p^c_{g,t} \in [0, \overline{P}^{ts}_{g,t}]\ (g\in\mathcal{G}^c). \tag{9}
```

**Storage** (identical to the [SOCWR (Jabr conic relaxation)](@ref) W-space storage model,
see the note under *Storage constraints*).

```math
E_{s,t}\in[E_s^{\min},E_s^{\max}],\
c_{s,t}\in[0,C_s^{\max}],\
d_{s,t}\in[0,D_s^{\max}],\
P^s_{s,t}\in[-\overline{S}_s,\overline{S}_s], \tag{10}
```

```math
Q^s_{s,t}\in[Q_s^{\min},Q_s^{\max}],\
Q^{sc}_{s,t}\in[Q_s^{\min},Q_s^{\max}],\
\mathrm{ccm}^s_{s,t}\in\big[0,(\overline{S}_s/V^{\min})^2\big],\
z_{s,t}\in\{0,1\}. \tag{11}
```

With storage present, the binary ``z_{s,t}`` lifts the instance to a mixed-integer NLP.
Without storage, QC is the pure NLP solved directly by Ipopt.

### Objective function

Identical to [AC Polar](@ref) (`src/modeling/objective/generation_cost.jl`): the fuel cost
on physical power ``S^{base}p^g`` plus the curtailment cost when curtailable units exist.

```math
F = \sum_{t\in\mathcal{T}}\sum_{g\in\mathcal{G}}\sum_{k=1}^{m_{g,t}}
  c_{g,t,k}\,(S^{base}p^g_{g,t})^{m_{g,t}-k}, \qquad
F^{curt} = \sum_{g\in\mathcal{G}^c}\sum_{t\in\mathcal{T}} c^{curt}_g\,S^{base}\,p^c_{g,t}, \tag{12}
```

```math
\min\; F + F^{curt}. \tag{13}
```

### Branch flows (W-space)

The four flows are the SOCWR/ACR expressions, reused verbatim
(`add_socwr_branch_expressions!`). The branch angle difference is tied to the bus angles.
QC has no separate ``\delta^{\min}\le\theta_i-\theta_j\le\delta^{\max}`` constraint. That
limit is enforced by the bound (5) on ``\mathrm{td}`` together with the exact coupling:

```math
\theta_{f(\ell),t} - \theta_{t(\ell),t} = \mathrm{td}_{\ell,t},
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{14}
```

Here ``\mathrm{td}`` is the raw angle difference ``\theta_{f(\ell)}-\theta_{t(\ell)}``,
with no phase-shift term; the transformer phase ``\Phi_\ell`` enters only through
``t_r, t_i`` in the flow coefficients below.

```math
P^f_{\ell,t} = \frac{(g_\ell+g^{fr}_\ell)}{tm_\ell}\,W_{ii}
  + \frac{(-g_\ell t_r + b_\ell t_i)}{tm_\ell}\,W_r
  + \frac{(-b_\ell t_r - g_\ell t_i)}{tm_\ell}\,W_i, \tag{15}
```

```math
Q^f_{\ell,t} = -\frac{(b_\ell+b^{fr}_\ell)}{tm_\ell}\,W_{ii}
  - \frac{(-b_\ell t_r - g_\ell t_i)}{tm_\ell}\,W_r
  + \frac{(-g_\ell t_r + b_\ell t_i)}{tm_\ell}\,W_i, \tag{16}
```

```math
P^t_{\ell,t} = (g_\ell+g^{to}_\ell)\,W_{jj}
  + \frac{(-g_\ell t_r - b_\ell t_i)}{tm_\ell}\,W_r
  + \frac{(b_\ell t_r - g_\ell t_i)}{tm_\ell}\,W_i, \tag{17}
```

```math
Q^t_{\ell,t} = -(b_\ell+b^{to}_\ell)\,W_{jj}
  - \frac{(-b_\ell t_r + g_\ell t_i)}{tm_\ell}\,W_r
  + \frac{(g_\ell t_r + b_\ell t_i)}{tm_\ell}\,W_i. \tag{18}
```

### Power balance

At each non-virtual bus (`non_virtual_bus_ids`), the shunt term rides on the lifted square
``W_{b,t}`` (not ``V_{b,t}^2``):

```math
\sum_{g\in\mathcal{G}(b)} p^g_{g,t}
  - \sum_{s\in\mathcal{S}(b)} P^s_{s,t}
  - P^d_{b,t}
  - \sum_{\ell:\,f(\ell)=b} P^f_{\ell,t}
  - \sum_{\ell:\,t(\ell)=b} P^t_{\ell,t}
  - g^{sh}_b\,W_{b,t} = 0, \qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{19}
```

```math
\sum_{g\in\mathcal{G}(b)} Q^g_{g,t}
  - \sum_{s\in\mathcal{S}(b)} Q^s_{s,t}
  - Q^d_{b,t}
  - \sum_{\ell:\,f(\ell)=b} Q^f_{\ell,t}
  - \sum_{\ell:\,t(\ell)=b} Q^t_{\ell,t}
  + b^{sh}_b\,W_{b,t} = 0, \qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{20}
```

### Reference bus

```math
\theta_{ref,t} = 0, \qquad t\in\mathcal{T}. \tag{21}
```

### Magnitude lifting (secant band)

For every bus (`add_qc_bus_variables!`): parabolic floor and secant ceiling (Constraint 1):

```math
V_{b,t}^2 \le W_{b,t} \le (V_b^{\min}+V_b^{\max})\,V_{b,t} - V_b^{\min}V_b^{\max},
  \qquad b\in\mathcal{N},\ t\in\mathcal{T}. \tag{22}
```

### Trigonometric envelopes

**Cosine** (Constraint 3), quadratic ceiling and secant floor, with
``m = \max(|\delta^{\min}_\ell|,|\delta^{\max}_\ell|)``:

```math
\mathrm{cs}_{\ell,t} \le 1 - \frac{1-\cos m}{m^2}\,\mathrm{td}_{\ell,t}^2, \tag{23}
```

```math
\mathrm{cs}_{\ell,t} \ge
  \frac{\cos\delta^{\min}_\ell - \cos\delta^{\max}_\ell}{\delta^{\min}_\ell - \delta^{\max}_\ell}
  \,(\mathrm{td}_{\ell,t} - \delta^{\min}_\ell) + \cos\delta^{\min}_\ell. \tag{24}
```

**Sine** (Constraint 4), three cases on the sign of the aperture, with the same ``m``.
Crossing zero (``\delta^{\min}<0<\delta^{\max}``): two tangents. All negative
(``\delta^{\max}\le 0``): secant ceiling, tangent floor. All positive
(``\delta^{\min}\ge 0``): tangent ceiling, secant floor.

```math
\begin{aligned}
\text{crosses }0:\quad
  &\mathrm{si} \le \cos\tfrac{m}{2}\big(\mathrm{td}-\tfrac{m}{2}\big)+\sin\tfrac{m}{2},
  &&\mathrm{si} \ge \cos\tfrac{m}{2}\big(\mathrm{td}+\tfrac{m}{2}\big)-\sin\tfrac{m}{2};\\[2pt]
\delta^{\max}\!\le 0:\quad
  &\mathrm{si} \le \tfrac{\sin\delta^{\min}-\sin\delta^{\max}}{\delta^{\min}-\delta^{\max}}
     (\mathrm{td}-\delta^{\min})+\sin\delta^{\min},
  &&\mathrm{si} \ge \cos\tfrac{m}{2}\big(\mathrm{td}+\tfrac{m}{2}\big)-\sin\tfrac{m}{2};\\[2pt]
\delta^{\min}\!\ge 0:\quad
  &\mathrm{si} \le \cos\tfrac{m}{2}\big(\mathrm{td}-\tfrac{m}{2}\big)+\sin\tfrac{m}{2},
  &&\mathrm{si} \ge \tfrac{\sin\delta^{\min}-\sin\delta^{\max}}{\delta^{\min}-\delta^{\max}}
     (\mathrm{td}-\delta^{\min})+\sin\delta^{\min}.
\end{aligned} \tag{25}
```

### McCormick products

For a product ``z=x\,y`` with ``x\in[\underline{x},\overline{x}]``,
``y\in[\underline{y},\overline{y}]`` the four envelope cuts are

```math
\begin{aligned}
z &\ge \underline{x}\,y + \underline{y}\,x - \underline{x}\,\underline{y}, &
z &\ge \overline{x}\,y + \overline{y}\,x - \overline{x}\,\overline{y},\\
z &\le \underline{x}\,y + \overline{y}\,x - \underline{x}\,\overline{y}, &
z &\le \overline{x}\,y + \underline{y}\,x - \overline{x}\,\underline{y}.
\end{aligned} \tag{26}
```

Instantiated three times (the recursive chain), in order, each per branch ``\ell`` and time
``t``:

| # | ``z`` | ``x`` (box) | ``y`` (box) |
|---|---|---|---|
| Product 1 | ``\mathrm{vv}_{\ell,t}`` | ``V_{i,t}\in[V_i^{\min},V_i^{\max}]`` | ``V_{j,t}\in[V_j^{\min},V_j^{\max}]`` |
| Product 2 | ``W_r`` | ``\mathrm{vv}_{\ell,t}\in[V_i^{\min}V_j^{\min},V_i^{\max}V_j^{\max}]`` | ``\mathrm{cs}_{\ell,t}\in[\mathrm{cs}^{\min},\mathrm{cs}^{\max}]`` |
| Product 3 | ``W_i`` | ``\mathrm{vv}_{\ell,t}`` (same box) | ``\mathrm{si}_{\ell,t}\in[\sin\delta^{\min},\sin\delta^{\max}]`` |

### Branch thermal limits

Posted as **quadratic** inequalities (QC is an NLP; unlike SOCWR these are not conic
constraints):

```math
(P^f_{\ell,t})^2 + (Q^f_{\ell,t})^2 \le \mathrm{rate}_{a,\ell}^2,\qquad
(P^t_{\ell,t})^2 + (Q^t_{\ell,t})^2 \le \mathrm{rate}_{a,\ell}^2,
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{27}
```

### Angle-difference tangent cuts

On the cross-blocks (Constraint 7), with the tangent argument clipped to
``\pm(\pi/2-10^{-6})``; each cut is posted only when the corresponding limit is not `NaN`:

```math
\tan(\delta^{\min}_\ell)\,W_r \le W_i \le \tan(\delta^{\max}_\ell)\,W_r,
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{28}
```

### Lifted nonlinear cuts (LNC)

The same two cuts as [SOCWR (Jabr conic relaxation)](@ref) (`add_qc_lnc_cuts!`), posted per
branch only when ``\delta^{\min}_\ell,\delta^{\max}_\ell`` are non-`NaN`, both lie in
``[-\pi/2,\pi/2]`` and ``\delta^{\min}_\ell<\delta^{\max}_\ell``. With
``\phi=\tfrac{\delta^{\max}+\delta^{\min}}{2}``, ``d=\tfrac{\delta^{\max}-\delta^{\min}}{2}``,
``s_f=V_f^{\min}+V_f^{\max}``, ``s_t=V_t^{\min}+V_t^{\max}``:

```math
s_f s_t(\cos\phi\,W_r+\sin\phi\,W_i)
  - V_t^{\max}\cos d\,s_t\,W_{ii} - V_f^{\max}\cos d\,s_f\,W_{jj}
  \ge V_f^{\max}V_t^{\max}\cos d\,(V_f^{\min}V_t^{\min}-V_f^{\max}V_t^{\max}), \tag{29}
```

```math
s_f s_t(\cos\phi\,W_r+\sin\phi\,W_i)
  - V_t^{\min}\cos d\,s_t\,W_{ii} - V_f^{\min}\cos d\,s_f\,W_{jj}
  \ge -V_f^{\min}V_t^{\min}\cos d\,(V_f^{\min}V_t^{\min}-V_f^{\max}V_t^{\max}). \tag{30}
```

### Current-magnitude tie (`ccm`)

The rotated-cone link and the algebraic closed form (Constraint 8), per branch and time:

```math
(P^f_{\ell,t})^2 + (Q^f_{\ell,t})^2 \le \frac{W_{ii}}{tm_\ell}\,\mathrm{ccm}_{\ell,t}, \tag{31}
```

```math
\mathrm{ccm}_{\ell,t} =
  (g_\ell^2+b_\ell^2)\Big(\frac{W_{ii}}{tm_\ell} + W_{jj}
    - \frac{2(t_r W_r + t_i W_i)}{tm_\ell}\Big)
  - (g^{fr}_\ell{}^2 + b^{fr}_\ell{}^2)\,\frac{W_{ii}}{tm_\ell}
  + 2\,(g^{fr}_\ell P^f_{\ell,t} - b^{fr}_\ell Q^f_{\ell,t}). \tag{32}
```

### Generation constraints

The active/reactive time envelopes (added only where strictly tighter than the scalar
bound) and the renewable curtailment equality, as in [AC Polar](@ref):

```math
p^g_{g,t} \le \overline{P}^{ts}_{g,t}\ (g\in\mathcal{G}\setminus\mathcal{G}^c),\qquad
p^g_{g,t} \ge \underline{P}^{ts}_{g,t}\ (g\in\mathcal{G}),\qquad
\underline{Q}^{ts}_{g,t} \le Q^g_{g,t} \le \overline{Q}^{ts}_{g,t}\ (g\in\mathcal{G}). \tag{33}
```

```math
p^g_{g,t} + p^c_{g,t} = \overline{P}^{ts}_{g,t}, \qquad g\in\mathcal{G}^c,\ t\in\mathcal{T}. \tag{34}
```

### Storage constraints

!!! note "QC storage is the SOCWR W-space model (not the AC Polar one)"
    QC borrows the SOCWR storage model. When a case carries storage, QC dispatches to
    `add_socwr_storage_constraints!` (`storage/dispatch.jl`), so storage is relaxed in
    W-space: non-simultaneity is a **binary big-M** (not the bilinear ``cd=0`` of
    [AC Polar](@ref)), and the inverter limit is a **rotated cone** on a squared-current
    variable ``\mathrm{ccm}^s`` (not the constant-``\overline{S}_s^2`` quadratic).

Energy balance (`_storage_energy_balance!`):

```math
E_{s,1} = E^0_s + \eta^c_s c_{s,1} - \frac{d_{s,1}}{\eta^d_s},\qquad
E_{s,t} = E_{s,t-1} + \eta^c_s c_{s,t} - \frac{d_{s,t}}{\eta^d_s}\ (t\ge 2). \tag{35}
```

Non-simultaneous charge/discharge via the binary ``z_{s,t}``:

```math
c_{s,t} \le z_{s,t}\,C_s^{\max}, \qquad
d_{s,t} \le (1-z_{s,t})\,D_s^{\max}. \tag{36}
```

Inverter limit (rotated cone) and converter AC coupling:

```math
(P^s_{s,t})^2 + (Q^s_{s,t})^2 \le W_{B(s),t}\,\mathrm{ccm}^s_{s,t}, \tag{37}
```

```math
P^s_{s,t} + (d_{s,t} - c_{s,t}) = P^{loss}_s + r_s\,\mathrm{ccm}^s_{s,t}, \qquad
Q^s_{s,t} = Q^{sc}_{s,t} + Q^{loss}_s + x_s\,\mathrm{ccm}^s_{s,t}. \tag{38}
```

### Complete model

```math
\min_{\substack{\theta,V,W,\,W_r,W_i,\mathrm{td},\mathrm{vv},\mathrm{cs},\mathrm{si},\mathrm{ccm},\\
  p^g,Q^g,p^c,\,E,c,d,P^s,Q^s,Q^{sc},\mathrm{ccm}^s,z}}
\quad F + F^{curt} \qquad (\text{Eq. } 13)
```

subject to:

| Block | Equations |
|---|---|
| Flows (coupling + W-space) | (14)–(18) |
| Balance | (19), (20) |
| Reference bus | (21) |
| Magnitude lifting | (22) |
| Trig envelopes | (23)–(25) |
| McCormick products | (26) ×3 |
| Thermal | (27) |
| Tangent / LNC cuts | (28)–(30) |
| Current tie | (31), (32) |
| Generation | (33), (34) |
| Storage (if present) | (35)–(38) |
| Bounds | (1)–(11) |

## What it captures / what it loses

!!! info "Philosophy"
    **Captures.** By opening the blocks, QC adds (i) the exact angle coupling
    ``\mathrm{td} = \theta_i - \theta_j``, which restores angle cycle consistency (the bus
    angles telescope to ``\sum\mathrm{td} = 0`` around any loop), and (ii) per-piece
    envelopes (secant/parabola on ``\mathrm{cs}``, tangents on ``\mathrm{si}``, recursive
    McCormick on the products, plus the squared-current tie ``\mathrm{ccm}``) in place of
    SOCWR's single voltage-product cone. Neither relaxation dominates the other in general.

    **Loses.** The link from ``\mathrm{td}`` to the blocks passes through *relaxed*
    envelopes, so the net only closes as far as the envelopes squeeze. QC remains a
    per-branch (per-bus-pair) relaxation: it does not impose the network-wide, rank-1
    condition that [SDPWRM](@ref) places on the whole voltage matrix.

## Lirion vs. the original

!!! warning
    Lirion implements the **recursive-McCormick (QCRM)** variant of the quadratic-convex
    relaxation of Hijazi, Coffrin and Van Hentenryck (2017), as provided in PowerModels
    [Coffrin2018](@cite). Two implementation notes:

    - The cross-blocks ``W_r, W_i`` (`R`, `T_var` in the code) are shared **per bus-pair**:
      parallel branches alias one product and each posts its McCormick chain onto that shared
      pair, which is a valid, tighter-or-equal relaxation.
    - The relaxation is assembled as an **NLP**, keeping the thermal limit (Constraint 6) as
      a quadratic constraint and the power-magnitude tie (Constraint 8) as a rotated-cone
      quadratic solved with Ipopt, rather than a purely conic or linear form. The LNC cuts
      (Constraint 9) are the lifted nonlinear cuts of Coffrin, Hijazi and Van Hentenryck
      (2015), the same block as the [SOCWR (Jabr conic relaxation)](@ref) model (also posted
      by [SDPWRM](@ref)).

## See also

- Neighbours: [SOCWR (Jabr conic relaxation)](@ref) · [SDPWRM](@ref) · [AC Rectangular (ACR)](@ref) · [AC Polar](@ref)
- Foundations: notation and network model in [Foundations](@ref)
- Reference implementation: constructor `qc()` (returns `QC`); source
  `src/modeling/branch/qc.jl`
- PowerModels equivalent: `QCRMPowerModel` [Coffrin2018](@cite)
