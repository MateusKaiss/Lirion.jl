# SDPWRM

The SDPWRM model is a convex relaxation written in **W-space**. The
[SOCWR (Jabr conic relaxation)](@ref) and [QC (Quadratic Convex)](@ref) relaxations tighten
the voltage-product blocks one branch at a time. SDPWRM gathers every block into a single
Hermitian matrix ``W = V\,V^{*}`` over **all** buses and imposes one condition on it,
positive semidefiniteness ``W \succeq 0``. That single condition already contains the whole
SOCWR relaxation as its ``2 \times 2`` slice and adds the loop couplings SOCWR misses, up to
the full network, so SDPWRM is never looser than SOCWR. It minimizes the same generation
cost as the exact [AC Polar](@ref) model, over the convex set of semidefinite ``W`` instead
of over genuine voltages.

!!! note "A convex relaxation"
    SDPWRM is a **relaxation**, not an approximation: it enlarges the feasible set (it drops
    the non-convex ``\operatorname{rank}(W) = 1`` requirement), so its optimum is a
    **guaranteed lower bound** on the true AC cost. A semidefinite ``W`` that is not rank one
    corresponds to no genuine voltage vector, so the bound can be loose, but it always errs
    low, never high.

## Derivation

### 1. The local tie does not close the cycles

In the [SOCWR (Jabr conic relaxation)](@ref) the tie between the voltage-product blocks is
**local**, holding branch by branch:

```math
W_r^2 + W_i^2 \;\le\; W_{ii}\,W_{jj}.
```

Each branch settles its own accounts without consulting its neighbours, and in a meshed
network that lets through a lie no single branch can catch. Each branch carries an angle
buried in its blocks, ``\theta_{ij} = \operatorname{atan2}(W_i, W_r)``, and nothing ties
these angles to one another. Take a triangle of three lines, each reporting a ``30^\circ``
difference. Branch by branch every one is impeccable, sitting exactly on the SOCWR cone.
Walk the loop, though, and the three differences add to
``30^\circ + 30^\circ + 30^\circ = 90^\circ``, when genuine bus angles would have to bring
them back to zero:

```math
(\theta_1 - \theta_2) + (\theta_2 - \theta_3) + (\theta_3 - \theta_1) = 0.
```

Three honest-looking branches, one impossible triangle. SOCWR never catches it, because it
never looks at more than one branch at a time. The worked example at the end of this page
carries the same triangle through in full.

The [QC (Quadratic Convex)](@ref) relaxation attacked this with the exact angle coupling
``td = \theta_i - \theta_j``, which brings the bus angles back so that they telescope, but
the link between ``td`` and the blocks runs through the (relaxed) envelopes, so it only
tightens as far as the envelopes squeeze. SDPWRM attacks from another direction:
instead of improving branch by branch, it imposes a single condition on the **entire
network**.

### 2. The central idea: the matrix W

Until now the blocks were indexed by branch (`wr`, `wi`) and by bus (`wii`). SDPWRM
organizes them into a matrix. Define the vector of all complex bus voltages,

```math
V = (V_1, V_2, \ldots, V_n),
```

and the matrix ``W`` as the **outer product** ``V V^{*}``, column times row. An inner
product compresses two vectors into a single number; the outer product does the reverse,
expanding them into an ``n \times n`` matrix:

```math
W = V\,V^{*}, \qquad\text{that is}\qquad W_{ij} = V_i\,V_j^{*}.
```

That entry ``V_i V_j^{*}`` is the object behind the blocks of the
[AC Rectangular (ACR)](@ref) model, ``V_i V_j^{*} = W_r + j\,W_i`` (the imaginary unit is
the plain ``j`` of the [AC Polar](@ref) page; bus indices are only subscripts, so the two
never collide). Hence:

- diagonal: ``W_{ii} = |V_i|^2``, the squared voltage magnitudes of every bus (the bus
  block `wii`);
- off-diagonal: ``W_{ij} = W_r + j\,W_i``, the cross blocks (`wr`, `wi`) of the pair
  ``ij``.

Every block of every branch and bus lives in a single object. The matrix has an entry for
every pair of buses, whether or not a line joins them.

Swapping the two indices conjugates the entry: ``W_{ji} = V_j V_i^{*}`` is the complex
conjugate of ``W_{ij} = V_i V_j^{*}``. A matrix whose mirror entries are conjugates like
this equals its own conjugate transpose, and is called **Hermitian**, the complex
counterpart of a symmetric matrix. Section 6 leans on exactly this structure when it
rebuilds ``W`` over the reals.

### 3. Rank and cycle consistency

The matrix lets us ask the one question SOCWR could not: do all these blocks come from a
single common ``V``? The tool that answers it is **rank**.

The rank of a matrix is the number of linearly independent columns, how many distinct
directions the columns actually occupy. A generic ``3 \times 3`` matrix has rank 3, three
columns pointing in distinct directions. If one column is a combination of the others, the
rank drops. The extreme case, rank 1, means every column points the same way: each is a
multiple of a single vector.

That extreme is our question answered yes. If ``W`` genuinely factors as ``V V^{*}``, every
column is a multiple of that one ``V``, so ``W`` has rank 1. A single generator exists, each
bus then carries one angle, and the differences telescope around any loop on their own:

```math
\operatorname{rank}(W) = 1 \;\Longrightarrow\; \text{one generator } V \;\Longrightarrow\; \text{every cycle closes}.
```

Rank 1 is the whole of genuineness. It locks the magnitudes and the angles together onto one
``V``. Cycle consistency is only its angle half. Keep that split in view. The relaxation is
about to drop rank 1, and the question becomes how much of the cycle structure survives
without it (Section 5).

### 4. Keep the convex half

The exact requirement is that ``W`` have the form ``V V^{*}`` for some ``V \neq 0``, which is
equivalent to two conditions, ``W`` **positive semidefinite** (meaning ``x^{*} W x \ge 0``
for every vector ``x``, written ``W \succeq 0``) and **rank one**:

```math
W = V\,V^{*}\ \text{ for some } V \neq 0 \quad\Longleftrightarrow\quad W \succeq 0 \ \text{ and }\ \operatorname{rank}(W) = 1.
```

In OPF the bus voltages are bounded away from zero, so ``V \neq 0`` always holds.

Any physical ``W`` passes the semidefinite test for free:

```math
x^{*} W x = x^{*} V V^{*} x = (V^{*}x)^{*}(V^{*}x) = |V^{*}x|^2 \;\ge\; 0,
```

a squared magnitude. So ``W \succeq 0`` is a **necessary** condition for ``W`` to come from
real voltages. Rank one, on the other hand, is non-convex. Imposing
``\operatorname{rank}(W) = 1`` would bring back the full difficulty of the AC problem.
SDPWRM keeps only the convex half of the pair:

```math
W \succeq 0, \qquad\text{i.e.}\qquad x^{*} W x \ge 0 \ \text{ for every vector } x,
```

with ``x^{*}`` the conjugate transpose (the vector laid down and its entries conjugated).

### 5. The ladder of principal minors

What does ``W \succeq 0`` actually impose? A **principal submatrix** is the block left when
you pick a set of indices and keep their rows and columns; its **minor** is that block's
determinant. The criterion is that ``W \succeq 0`` if and only if **all** principal minors
are ``\ge 0``. This unfolds the single constraint into a ladder (one rung per block size),
and each rung is an old acquaintance.

**Rung 1 (``1 \times 1`` blocks).** The diagonal entries: ``W_{ii} \ge 0``.

**Rung 2 (``2 \times 2`` blocks), the whole SOCWR.** The block on ``\{i,j\}`` is

```math
\begin{bmatrix} W_{ii} & W_r + j\,W_i \\ W_r - j\,W_i & W_{jj} \end{bmatrix},
```

and its determinant ``\ge 0`` gives ``W_{ii}W_{jj} - (W_r^2 + W_i^2) \ge 0``, that is
``W_r^2 + W_i^2 \le W_{ii}W_{jj}``. The SOCWR tie is the ``2 \times 2`` minor of the SDP. It
was not replaced, it was **absorbed**, and it follows for free that the SDP is never looser
than SOCWR.

**Rung 3 (``3 \times 3`` blocks), the cycle constraint.** Here is the model's novelty, read
off at the SOCWR cone boundary. Take the triangle ``\{1,2,3\}`` and set every magnitude to
1, so the minor speaks about the angles alone. Each branch then picks its own angle through
a unit-modulus block, ``W_{12} = e^{j\alpha}``, ``W_{23} = e^{j\beta}``,
``W_{31} = e^{j\gamma}``, with ``\alpha, \beta, \gamma`` free and not yet tied to one
another. The Hermitian submatrix is

```math
\begin{bmatrix}
1 & e^{j\alpha} & e^{-j\gamma} \\
e^{-j\alpha} & 1 & e^{j\beta} \\
e^{j\gamma} & e^{-j\beta} & 1
\end{bmatrix}.
```

Writing ``\xi = \alpha + \beta + \gamma``, the determinant is

```math
\det = e^{j\xi} + e^{-j\xi} - 2 = 2(\cos\xi - 1),
```

using Euler's identity ``e^{j\xi} + e^{-j\xi} = 2\cos\xi``. Since ``\cos\xi \le 1`` always,
the determinant can never be positive, while PSD demands ``\det \ge 0``. Squeezed between
the formula's natural ceiling and the constraint's floor, it survives only at ``\det = 0``,
which forces ``\cos\xi = 1``, that is ``\xi = 0`` once the per-branch apertures stay inside
``\pm 90^\circ``. So on this boundary slice **rung 3 is the cycle constraint**,
``\alpha + \beta + \gamma = 0``, and only now does a common generator exist,

```math
V = \begin{bmatrix} 1 \\ e^{-j\alpha} \\ e^{-j(\alpha+\beta)} \end{bmatrix},
\qquad \gamma = -(\alpha + \beta),
```

with the submatrix equal to ``V V^{*}`` and therefore rank one.

And the ladder continues. ``4 \times 4`` couples quadrilaterals, and so on up to the
``n \times n`` minor, the whole network. ``W \succeq 0`` is that entire ladder at once,
packed into a single convex constraint.

### 6. The real 2n × 2n lifting

``W`` is complex; SDP solvers speak real numbers. Split
``W = W^{\mathrm{re}} + j\,W^{\mathrm{im}}`` and assemble the real matrix

```math
M = \begin{bmatrix} W^{\mathrm{re}} & W^{\mathrm{im}} \\ -W^{\mathrm{im}} & W^{\mathrm{re}} \end{bmatrix} \succeq 0.
```

Because ``W`` is Hermitian, ``W^{\mathrm{re}}`` is symmetric and ``W^{\mathrm{im}}`` is
antisymmetric. The bottom-left block ``-W^{\mathrm{im}}`` equals
``(W^{\mathrm{im}})^{\top}``, so ``M`` comes out symmetric, the form an SDP solver accepts.
Testing ``M \succeq 0`` is exactly testing ``W \succeq 0``.

!!! note "Proof: M is symmetric, and the two quadratic forms agree"
    **Why ``M`` is symmetric.** Symmetry means the entry at ``(p,q)`` equals the entry at
    ``(q,p)``. Take an entry of the top-right block, row ``i``, column ``n+j``:

    ```math
    M[i,\, n+j] = W^{\mathrm{im}}_{ij}.
    ```

    The mirror entry, in the transposed position (row ``n+j``, column ``i``), lies in the
    bottom-left block ``-W^{\mathrm{im}}`` at position ``(j,i)``:

    ```math
    M[n+j,\, i] = (-W^{\mathrm{im}})_{ji} = -W^{\mathrm{im}}_{ji}.
    ```

    Since ``W^{\mathrm{im}}`` is antisymmetric, ``W^{\mathrm{im}}_{ji} = -W^{\mathrm{im}}_{ij}``, so

    ```math
    M[n+j,\, i] = -(-W^{\mathrm{im}}_{ij}) = +W^{\mathrm{im}}_{ij} = M[i,\, n+j].
    ```

    **Why ``x^{*} W x = v^{\top} M v``.** Write ``x = a - j\,b`` (so
    ``x^{*} = a^{\top} + j\,b^{\top}``), ``W = W^{\mathrm{re}} + j\,W^{\mathrm{im}}``, and
    ``v = (a;\, b)``.

    *Left side.* First,

    ```math
    W x = (W^{\mathrm{re}} + j\,W^{\mathrm{im}})(a - j\,b)
    = \underbrace{(W^{\mathrm{re}}a + W^{\mathrm{im}}b)}_{R} + j\underbrace{(W^{\mathrm{im}}a - W^{\mathrm{re}}b)}_{I},
    ```

    ```math
    x^{*} W x = (a^{\top} + j\,b^{\top})(R + j\,I)
    = \underbrace{(a^{\top}R - b^{\top}I)}_{\text{real part}} + j\underbrace{(a^{\top}I + b^{\top}R)}_{\text{imaginary part}}.
    ```

    The imaginary part dies. With
    ``a^{\top}I = a^{\top}W^{\mathrm{im}}a - a^{\top}W^{\mathrm{re}}b`` and
    ``b^{\top}R = b^{\top}W^{\mathrm{re}}a + b^{\top}W^{\mathrm{im}}b``, each surviving term
    is zero: ``a^{\top}W^{\mathrm{im}}a = 0`` and ``b^{\top}W^{\mathrm{im}}b = 0`` by
    antisymmetry (a number equal to its own negative, since
    ``a^{\top}W^{\mathrm{im}}a = a^{\top}(W^{\mathrm{im}})^{\top}a = -a^{\top}W^{\mathrm{im}}a``),
    and ``-a^{\top}W^{\mathrm{re}}b + b^{\top}W^{\mathrm{re}}a = 0`` because
    ``b^{\top}W^{\mathrm{re}}a = a^{\top}(W^{\mathrm{re}})^{\top}b = a^{\top}W^{\mathrm{re}}b``
    (transposition rule ``(XYZ)^{\top} = Z^{\top}Y^{\top}X^{\top}`` and
    ``W^{\mathrm{re}}`` symmetric). The real part is then

    ```math
    x^{*} W x = a^{\top}W^{\mathrm{re}}a + a^{\top}W^{\mathrm{im}}b - b^{\top}W^{\mathrm{im}}a + b^{\top}W^{\mathrm{re}}b.
    ```

    *Right side.* With ``M`` and ``v`` as above,

    ```math
    M v = \begin{bmatrix} W^{\mathrm{re}}a + W^{\mathrm{im}}b \\ -W^{\mathrm{im}}a + W^{\mathrm{re}}b \end{bmatrix},
    ```

    ```math
    v^{\top} M v = a^{\top}(W^{\mathrm{re}}a + W^{\mathrm{im}}b) + b^{\top}(-W^{\mathrm{im}}a + W^{\mathrm{re}}b)
    = a^{\top}W^{\mathrm{re}}a + a^{\top}W^{\mathrm{im}}b - b^{\top}W^{\mathrm{im}}a + b^{\top}W^{\mathrm{re}}b.
    ```

    Identical term by term. One is PSD if and only if the other is, so sending ``M`` to the
    solver imposes the condition on ``W``.

## SDPWRM Formulation

SDPWRM reuses the shared generation, balance and objective blocks unchanged. What is
specific to it is the pair of real voltage matrices `WR_sdp`/`WI_sdp` and the single
positive-semidefinite constraint that couples them. The W-space branch flows are produced
by the same `add_socwr_branch_expressions!` as SOCWR and QC, with `R`, `T_var` and
`v_magnitude_sq` aliasing matrix entries.

Without storage the problem is a pure semidefinite program, which Lirion solves with
Clarabel by default. The storage block adds a charge/discharge binary (Eq. 20), so a
storage-bearing instance becomes a mixed-integer SDP and needs a solver that handles both
the semidefinite cone and integer variables. The benchmark cases carry no storage.

### Sets

- ``\mathcal{N}``: buses; ``n = |\mathcal{N}|``. A fixed integer index
  ``\operatorname{idx}(b) \in \{1,\dots,n\}`` places each bus in the matrices; the matrices
  ``W^{\mathrm{re}}, W^{\mathrm{im}}`` span all ``n`` buses.
- ``\mathcal{N}^{nv} \subseteq \mathcal{N}``: non-virtual buses, on which the power balance
  is enforced (`non_virtual_bus_ids`).
- ``\mathcal{L}``: active branches (in-service lines and transformers).
- ``\mathcal{G}``: generators; ``\mathcal{G}^c \subseteq \mathcal{G}`` the curtailable subset.
- ``\mathcal{S}``: storage units.
- ``\mathcal{T} = \{1,\dots,T\}``: time periods.

### Indices and incidence functions

- ``b, i, j \in \mathcal{N}``; ``\ell \in \mathcal{L}``; ``g \in \mathcal{G}``;
  ``s \in \mathcal{S}``; ``t \in \mathcal{T}``.
- ``f(\ell),\ t(\ell)``: *from* and *to* buses of branch ``\ell``; their matrix
  rows/columns are ``i_\ell = \operatorname{idx}(f(\ell))`` and
  ``j_\ell = \operatorname{idx}(t(\ell))``.
- ``B(s)``: connection bus of storage ``s``.

### Parameters

**Branches.** ``g_\ell, b_\ell`` (series conductance/susceptance);
``g^{fr}_\ell, b^{fr}_\ell`` and ``g^{to}_\ell, b^{to}_\ell`` (shunt, *from* / *to* sides);
``k_\ell = 1/\tau_\ell`` the inverse tap ratio (`inv_tap`);
``\tau^{re}_\ell = \tau_\ell\cos\Phi_\ell`` (`tap_real`) and
``\tau^{im}_\ell = \tau_\ell\sin\Phi_\ell`` (`tap_imag`), the real/imaginary parts of the
complex tap; ``\overline{S}_\ell`` the apparent-flow limit (`rate_a`);
``\theta^{\min}_\ell, \theta^{\max}_\ell`` the angular-aperture limits
(`ang_min`, `ang_max`).

**Buses.** ``V_b^{\min}, V_b^{\max}`` voltage-magnitude limits; ``g^{sh}_b, b^{sh}_b``
shunt conductance/susceptance; ``P^d_{b,t}, Q^d_{b,t}`` active/reactive demand.

**Generators.** ``P_g^{\min}, P_g^{\max}`` and ``Q_g^{\min}, Q_g^{\max}`` static limits;
active/reactive time envelopes ``\overline{P}^{ts}_{g,t}, \underline{P}^{ts}_{g,t}`` and
``\overline{Q}^{ts}_{g,t}, \underline{Q}^{ts}_{g,t}``; cost coefficients ``c_{g,t,k}``;
unit curtailment cost ``c^{curt}_g``.

**Storage.** ``E_s^{\min}, E_s^{\max}`` energy limits; initial energy ``E^0_s``;
charge/discharge power caps ``C_s^{\max}, D_s^{\max}``; efficiencies ``\eta^c_s, \eta^d_s``;
converter resistance and reactance ``r_s, x_s``; fixed losses ``P^{loss}_s, Q^{loss}_s``;
inverter thermal limit ``\overline{S}_s``; reactive limits ``Q_s^{\min}, Q_s^{\max}``.

### Decision variables

**Voltage (W-space).** Two real matrices per period replace the magnitude/angle pair.
There is no voltage-angle variable, and hence no reference-bus equation, since ``W`` is
invariant to a global phase rotation.

```math
W^{\mathrm{re}}_t \in \mathbb{S}^{n}, \qquad
W^{\mathrm{im}}_t \in \mathbb{R}^{n\times n}, \qquad t\in\mathcal{T}. \tag{1}
```

``W^{\mathrm{re}}_t`` is declared **symmetric** in the code (`Symmetric`); ``W^{\mathrm{im}}_t``
is forced **antisymmetric** by the PSD embedding (Eq. 13). The diagonal
``W^{\mathrm{re}}_{t,ii}`` is the squared voltage magnitude of the bus at index ``i``,
exposed as the expression ``w_{b,t} :=`` `v_magnitude_sq`
``\,= W^{\mathrm{re}}_{t,\operatorname{idx}(b)\,\operatorname{idx}(b)}``. Both matrices carry a
symmetry inherited from trigonometry:

- ``W^{\mathrm{re}}`` is **symmetric**: ``W^{\mathrm{re}}_{ij} = v_i v_j \cos\theta_{ij}``;
  swapping ``i \leftrightarrow j`` commutes the magnitudes and sends
  ``\theta_{ji} = -\theta_{ij}``, but the cosine is even, so the value is unchanged.
- ``W^{\mathrm{im}}`` is **antisymmetric**: ``W^{\mathrm{im}}_{ij} = v_i v_j \sin\theta_{ij}``,
  and the sine is odd, so ``W^{\mathrm{im}}_{ji} = -W^{\mathrm{im}}_{ij}``; the diagonal is
  zero (``\sin 0 = 0``).

**Generation.**

```math
p^g_{g,t} \in [P_g^{\min}, P_g^{\max}], \qquad
Q^g_{g,t} \in [Q_g^{\min}, Q_g^{\max}], \qquad g\in\mathcal{G},\ t\in\mathcal{T}. \tag{2}
```

```math
p^c_{g,t} \in [0, \overline{P}^{ts}_{g,t}], \qquad g\in\mathcal{G}^c,\ t\in\mathcal{T}. \tag{3}
```

**Storage** (the W-space block shared with SOCWR).

```math
E_{s,t}\in[E_s^{\min},E_s^{\max}],\quad
c_{s,t}\in[0,C_s^{\max}],\quad
d_{s,t}\in[0,D_s^{\max}],\quad
P^s_{s,t}\in[-\overline{S}_s,\overline{S}_s], \tag{4}
```

```math
Q^s_{s,t}\in[Q_s^{\min},Q_s^{\max}],\quad
Q^{sc}_{s,t}\in[Q_s^{\min},Q_s^{\max}],\quad
z_{s,t}\in\{0,1\},\quad
\ell^{cm}_{s,t}\in\big[0,\ (\overline{S}_s/V^{\min})^2\big], \tag{5}
```

with ``z_{s,t}`` the charge/discharge selector and ``\ell^{cm}_{s,t}`` the squared
converter current magnitude (`ccms`); ``V^{\min}`` is a global voltage floor
(`model.ext[:v_min]`, default ``0.9``).

### Entry bounds

**Diagonal, from the voltage limits:**

```math
(V_i^{\min})^2 \;\le\; W^{\mathrm{re}}_{ii} \;\le\; (V_i^{\max})^2.
```

**Off-diagonal, from the voltage and angle limits**, and only for pairs joined by a branch
(a pair with no line carries no angle information; the entry exists, but is held only by the
PSD constraint). These are the same three corner cases derived from the ``W_r``, ``W_i``
bounds of the [QC (Quadratic Convex)](@ref) (the corner principle, with the sign of the
sine deciding the pairing).

If ``\theta^{\min} \ge 0``:

```math
V_i^{\min} V_j^{\min}\cos\theta^{\max} \;\le\; W^{\mathrm{re}}_{ij} \;\le\; V_i^{\max} V_j^{\max}\cos\theta^{\min},
```

```math
V_i^{\min} V_j^{\min}\sin\theta^{\min} \;\le\; W^{\mathrm{im}}_{ij} \;\le\; V_i^{\max} V_j^{\max}\sin\theta^{\max}.
```

If ``\theta^{\max} \le 0``:

```math
V_i^{\min} V_j^{\min}\cos\theta^{\min} \;\le\; W^{\mathrm{re}}_{ij} \;\le\; V_i^{\max} V_j^{\max}\cos\theta^{\max},
```

```math
V_i^{\max} V_j^{\max}\sin\theta^{\min} \;\le\; W^{\mathrm{im}}_{ij} \;\le\; V_i^{\min} V_j^{\min}\sin\theta^{\max}.
```

If the interval crosses zero:

```math
V_i^{\min} V_j^{\min}\min\!\left(\cos\theta^{\min},\ \cos\theta^{\max}\right) \;\le\; W^{\mathrm{re}}_{ij} \;\le\; V_i^{\max} V_j^{\max},
```

```math
V_i^{\max} V_j^{\max}\sin\theta^{\min} \;\le\; W^{\mathrm{im}}_{ij} \;\le\; V_i^{\max} V_j^{\max}\sin\theta^{\max}.
```

**Mirror.** The other triangle of the matrix shares these bounds through the symmetries:
``W^{\mathrm{re}}_{ji}`` inherits the same interval; ``W^{\mathrm{im}}_{ji} = -W^{\mathrm{im}}_{ij}``
inherits the reflected one (if ``W^{\mathrm{im}}_{ij} \in [a,b]`` then
``W^{\mathrm{im}}_{ji} \in [-b,-a]``).

**Parallel lines.** If two or three branches join the same pair of buses, each brings its
own limits but the matrix entry is a single one (indexing is by bus, not by branch). The
resolution is the **intersection** of the intervals (the smallest of the ceilings, the
largest of the floors). Each branch tightens the same entry, and the joint squeeze holds.

### Objective function

The objective is the shared generation cost, identical to the exact [AC Polar](@ref)
model: the fuel-cost polynomial evaluated on power in physical units (``S^{base}p^g``) plus,
when curtailable generators are present, a linear curtailment cost.

```math
F = \sum_{t\in\mathcal{T}}\sum_{g\in\mathcal{G}}\sum_{k=1}^{m_{g,t}}
c_{g,t,k}\,(S^{base}\,p^g_{g,t})^{m_{g,t}-k}, \qquad
F^{curt} = \sum_{g\in\mathcal{G}^c}\sum_{t\in\mathcal{T}} c^{curt}_g\,S^{base}\,p^c_{g,t},
```

```math
\min\; F + F^{curt}. \tag{6}
```

### Branch flows

For branch ``\ell`` write the aliased blocks: the two endpoint diagonals
``W^{ff}_{\ell,t}=W^{\mathrm{re}}_{t,i_\ell i_\ell}`` and
``W^{tt}_{\ell,t}=W^{\mathrm{re}}_{t,j_\ell j_\ell}``, and the off-diagonal pair
``W^{r}_{\ell,t}=W^{\mathrm{re}}_{t,i_\ell j_\ell}`` (= `R`) and
``W^{i}_{\ell,t}=W^{\mathrm{im}}_{t,i_\ell j_\ell}`` (= `T_var`). The four flows are
**linear** in these blocks:

```math
P^f_{\ell,t} = (g_\ell+g^{fr}_\ell)\,k_\ell^2\,W^{ff}_{\ell,t}
+ (-g_\ell\tau^{re}_\ell + b_\ell\tau^{im}_\ell)\,k_\ell^2\,W^{r}_{\ell,t}
+ (-b_\ell\tau^{re}_\ell - g_\ell\tau^{im}_\ell)\,k_\ell^2\,W^{i}_{\ell,t}. \tag{7}
```

```math
Q^f_{\ell,t} = -(b_\ell+b^{fr}_\ell)\,k_\ell^2\,W^{ff}_{\ell,t}
- (-b_\ell\tau^{re}_\ell - g_\ell\tau^{im}_\ell)\,k_\ell^2\,W^{r}_{\ell,t}
+ (-g_\ell\tau^{re}_\ell + b_\ell\tau^{im}_\ell)\,k_\ell^2\,W^{i}_{\ell,t}. \tag{8}
```

```math
P^t_{\ell,t} = (g_\ell+g^{to}_\ell)\,W^{tt}_{\ell,t}
+ (-g_\ell\tau^{re}_\ell - b_\ell\tau^{im}_\ell)\,k_\ell^2\,W^{r}_{\ell,t}
- (-b_\ell\tau^{re}_\ell + g_\ell\tau^{im}_\ell)\,k_\ell^2\,W^{i}_{\ell,t}. \tag{9}
```

```math
Q^t_{\ell,t} = -(b_\ell+b^{to}_\ell)\,W^{tt}_{\ell,t}
- (-b_\ell\tau^{re}_\ell + g_\ell\tau^{im}_\ell)\,k_\ell^2\,W^{r}_{\ell,t}
- (-g_\ell\tau^{re}_\ell - b_\ell\tau^{im}_\ell)\,k_\ell^2\,W^{i}_{\ell,t}. \tag{10}
```

Note the inverse tap ``k_\ell^2`` scales only the *from*-side self term (Eq. 7, 8) and all
cross terms; the *to*-side self term (Eq. 9, 10) carries no tap factor.

### Power balance

At each non-virtual bus, injections meet withdrawals; the shunt draws on the diagonal
block ``w_{b,t}=W^{\mathrm{re}}_{t,\operatorname{idx}(b)\,\operatorname{idx}(b)}``.

**Active power.**

```math
\sum_{g\in\mathcal{G}(b)} p^g_{g,t}
- \sum_{s\in\mathcal{S}(b)} P^s_{s,t}
- P^d_{b,t}
- \sum_{\ell:\,f(\ell)=b} P^f_{\ell,t}
- \sum_{\ell:\,t(\ell)=b} P^t_{\ell,t}
- g^{sh}_b\,w_{b,t} = 0, \qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{11}
```

**Reactive power.**

```math
\sum_{g\in\mathcal{G}(b)} Q^g_{g,t}
- \sum_{s\in\mathcal{S}(b)} Q^s_{s,t}
- Q^d_{b,t}
- \sum_{\ell:\,f(\ell)=b} Q^f_{\ell,t}
- \sum_{\ell:\,t(\ell)=b} Q^t_{\ell,t}
+ b^{sh}_b\,w_{b,t} = 0, \qquad b\in\mathcal{N}^{nv},\ t\in\mathcal{T}. \tag{12}
```

The two shunt terms carry opposite signs, as they do throughout Lirion: the conductance
draws active power (``-g^{sh}_b w_{b,t}``), while the stored susceptance
(`pf.shunt_susceptance`) enters reactive balance with a plus (``+b^{sh}_b w_{b,t}``), the
sign the code posts.

### Positive-semidefinite constraint

The single condition that defines the relaxation. The code assembles **one global**
``2n\times 2n`` real symmetric block over **all** ``n`` buses (no clique or chordal
decomposition) and passes it to the solver, per period:

```math
M_t = \begin{bmatrix} W^{\mathrm{re}}_t & W^{\mathrm{im}}_t \\
-W^{\mathrm{im}}_t & W^{\mathrm{re}}_t \end{bmatrix} \succeq 0, \qquad t\in\mathcal{T}. \tag{13}
```

Requiring ``M_t`` symmetric is what forces ``W^{\mathrm{im}}_t`` antisymmetric; testing
``M_t \succeq 0`` is exactly testing the Hermitian
``W_t = W^{\mathrm{re}}_t + j\,W^{\mathrm{im}}_t \succeq 0`` (Section 6 of the derivation).
This global block absorbs the per-pair ``2\times 2`` Jabr cone of
[SOCWR (Jabr conic relaxation)](@ref), which is therefore **not** posted separately here.

### Branch thermal limits

Two second-order cones per branch, on the from- and to-side flows:

```math
(\overline{S}_\ell,\ P^f_{\ell,t},\ Q^f_{\ell,t}) \in \text{SecondOrderCone}, \qquad
(\overline{S}_\ell,\ P^t_{\ell,t},\ Q^t_{\ell,t}) \in \text{SecondOrderCone},
\qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{14}
```

### Angle limits

A linear cut on the blocks (there is no ``\theta`` variable, so no raw angle-difference
constraint is imposed, unlike SOCWR). Each side is posted only where its own aperture limit
is set: the upper cut where ``\theta^{\max}_\ell`` is present, the lower cut where
``\theta^{\min}_\ell`` is present:

```math
\tan(\theta^{\min}_\ell)\,W^{r}_{\ell,t} \;\le\; W^{i}_{\ell,t} \;\le\;
\tan(\theta^{\max}_\ell)\,W^{r}_{\ell,t}, \qquad \ell\in\mathcal{L},\ t\in\mathcal{T}. \tag{15}
```

The tangent argument is clipped to ``\pm(\pi/2 - 10^{-6})`` (`_tan_clip`) so the slope stays
finite as the aperture approaches ``\pm 90^\circ``.

### Lifted nonlinear cuts (LNC)

The same two cuts per branch as [SOCWR (Jabr conic relaxation)](@ref) (Coffrin–Hijazi–Van
Hentenryck 2015), added only when the aperture is finite, lies within ``\pm\pi/2`` and
satisfies ``\theta^{\min}_\ell<\theta^{\max}_\ell``. With
``\phi=(\theta^{\max}_\ell+\theta^{\min}_\ell)/2``,
``d=(\theta^{\max}_\ell-\theta^{\min}_\ell)/2``, ``s_f=V^{\min}_f+V^{\max}_f`` and
``s_t=V^{\min}_t+V^{\max}_t``:

```math
s_f s_t\,(\cos\phi\,W^{r}_{\ell,t}+\sin\phi\,W^{i}_{\ell,t})
- V^{\max}_t\cos d\,s_t\,W^{ff}_{\ell,t}
- V^{\max}_f\cos d\,s_f\,W^{tt}_{\ell,t}
\;\ge\; V^{\max}_f V^{\max}_t\cos d\,(V^{\min}_f V^{\min}_t - V^{\max}_f V^{\max}_t). \tag{16}
```

```math
s_f s_t\,(\cos\phi\,W^{r}_{\ell,t}+\sin\phi\,W^{i}_{\ell,t})
- V^{\min}_t\cos d\,s_t\,W^{ff}_{\ell,t}
- V^{\min}_f\cos d\,s_f\,W^{tt}_{\ell,t}
\;\ge\; -V^{\min}_f V^{\min}_t\cos d\,(V^{\min}_f V^{\min}_t - V^{\max}_f V^{\max}_t). \tag{17}
```

They are kept because the PSD ladder knows the *structure* (the blocks of a common
generator) but not the voltage and angle *limits*. No minor carries that information, and
the LNC cuts inject it. Accumulation, as always: each piece cuts a slack the others do not
see.

### Generation constraints

Active/reactive time envelopes are added only where a period limit is strictly tighter than
the static bound. The per-period upper active bound is skipped for renewables, which take
the curtailment equality instead; the lower active bound and the reactive envelope apply to
every generator.

```math
\underline{P}^{ts}_{g,t} \le p^g_{g,t} \le \overline{P}^{ts}_{g,t}\ (g\notin\mathcal{G}^c),
\quad
\underline{Q}^{ts}_{g,t} \le Q^g_{g,t} \le \overline{Q}^{ts}_{g,t},
\quad
p^g_{g,t}+p^c_{g,t}=\overline{P}^{ts}_{g,t}\ (g\in\mathcal{G}^c). \tag{18}
```

### Storage constraints (W-space)

SDPWRM has no storage file of its own. The dispatch sends storage to the SOCWR W-space
block `add_socwr_storage_constraints!`. Non-simultaneity is the binary big-M pair (Eq. 20),
not the bilinear ``c\,d = 0`` of AC Polar, and the converter-current bound rests on the
single global ``V^{\min}`` floor. **Energy balance:**

```math
E_{s,1} = E^0_s + \eta^c_s c_{s,1} - \frac{d_{s,1}}{\eta^d_s}, \qquad
E_{s,t} = E_{s,t-1} + \eta^c_s c_{s,t} - \frac{d_{s,t}}{\eta^d_s}\ (t\ge 2),
\qquad s\in\mathcal{S}. \tag{19}
```

**Charge/discharge exclusivity** (binary big-M by the ramp caps):

```math
c_{s,t} \le z_{s,t}\,C_s^{\max}, \qquad
d_{s,t} \le (1-z_{s,t})\,D_s^{\max}, \qquad s\in\mathcal{S},\ t\in\mathcal{T}. \tag{20}
```

**Converter thermal limit** (a rotated cone through the connection-bus diagonal
``w_{B(s),t}``) and the active/reactive **couplings** in the current-magnitude variable
``\ell^{cm}_{s,t}``:

```math
(P^s_{s,t})^2 + (Q^s_{s,t})^2 \;\le\; w_{B(s),t}\,\ell^{cm}_{s,t},
\qquad s\in\mathcal{S},\ t\in\mathcal{T}. \tag{21}
```

```math
P^s_{s,t} + (d_{s,t}-c_{s,t}) = P^{loss}_s + r_s\,\ell^{cm}_{s,t}, \qquad
Q^s_{s,t} = Q^{sc}_{s,t} + Q^{loss}_s + x_s\,\ell^{cm}_{s,t},
\qquad s\in\mathcal{S},\ t\in\mathcal{T}. \tag{22}
```

### Complete model

```math
\min_{\substack{W^{\mathrm{re}},\,W^{\mathrm{im}},\,p^g,\,Q^g,\,p^c,\\
E,\,c,\,d,\,P^s,\,Q^s,\,Q^{sc},\,z,\,\ell^{cm}}}
\quad F + F^{curt} \qquad (\text{Eq. } 6)
```

subject to:

| Block | Equations |
|---|---|
| W-space flows | (7)–(10) |
| Power balance | (11), (12) |
| PSD coupling | (13) |
| Thermal | (14) |
| Angle cut | (15) |
| LNC cuts | (16), (17) |
| Generation | (18) |
| Storage | (19)–(22) |
| Matrix-entry and variable bounds | (1)–(5) + Entry bounds |

## What it captures / what it loses

!!! info "Philosophy"
    **Captures:** a single convex tie that couples every loop at once. The ladder of
    principal minors runs from the ``2 \times 2`` SOCWR cone through the ``3 \times 3`` loop
    block up to the whole ``n \times n`` network, so SDPWRM is never looser than SOCWR and
    can be strictly tighter wherever a loop is present. What the tie does not do is force the
    loop angles to sum to zero at every feasible point. It pins a cycle exactly only when its
    branches sit on the SOCWR cone boundary (the unit-modulus case of rung 3), and in full
    generality only in the rank-one limit of the next bullet. A feasible ``W \succeq 0`` in
    the cone interior can still carry loop angles that do not close.

    **Loses:** the rank-one requirement. A feasible ``W \succeq 0`` need not factor as
    ``V V^{*}``, so it may match no genuine set of bus voltages. The optimum is a lower bound
    on the true AC cost, not a point the network can actually run.

## Lirion vs. the original

!!! warning
    The semidefinite program is handled on the real symmetric lifting ``M`` (the
    ``2n \times 2n`` matrix of Section 6), which is the standard way to pass a Hermitian PSD
    constraint to a real SDP solver. Lirion assembles one dense ``2n \times 2n`` block over
    all ``n`` buses and does not exploit chordal sparsity, the decomposition into smaller
    overlapping cliques that large-scale SDP-OPF codes lean on. The dense form stays faithful
    to the whole matrix and keeps the model simple, at the cost of how large a network solves
    in reasonable time. The convex semidefinite relaxation of the OPF is due to
    Bai, Wei, Fujisawa and Wang (2008), who first cast the problem over a semidefinite
    voltage matrix. Lavaei and Low (2012) later gave conditions under which the relaxation is
    exact (zero duality gap).

## Worked example

The guiding example runs throughout: three buses in a triangle (branches 1-2, 2-3, 3-1).

```@raw html
<figure style="text-align:center;margin:1.6rem 0">
<svg viewBox="0 0 680 360" role="img" aria-label="Two three-bus triangles. Left: each branch carries a 30 degree difference and the loop sums to 90 degrees, impossible for real bus angles. Right: genuine bus angles 0, minus 30 and minus 60 degrees give differences 30, 30 and minus 60 that sum to 0 and close the loop." style="max-width:640px;width:100%;height:auto;color:inherit" xmlns="http://www.w3.org/2000/svg" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round" font-family="Georgia, 'Times New Roman', serif">
  <defs>
    <marker id="lrn-ah" markerUnits="userSpaceOnUse" viewBox="0 0 12 12" refX="10" refY="6" markerWidth="12" markerHeight="12" orient="auto">
      <path d="M2,2.5 L11,6 L2,9.5 L4.5,6 Z" fill="currentColor" stroke="none"/>
    </marker>
  </defs>

  <line x1="340" y1="96" x2="340" y2="300" stroke-width="1.2" stroke-dasharray="3 6" stroke-opacity="0.22"/>

  <text x="170" y="32" text-anchor="middle" fill="currentColor" stroke="none" font-size="16.5">Three 30&#176; differences</text>
  <polyline points="170,92 125,172 80,252" marker-mid="url(#lrn-ah)"/>
  <polyline points="80,252 170,252 260,252" marker-mid="url(#lrn-ah)"/>
  <polyline points="260,252 215,172 170,92" marker-mid="url(#lrn-ah)"/>
  <circle cx="170" cy="92" r="6" fill="currentColor" stroke="none"/>
  <circle cx="80" cy="252" r="6" fill="currentColor" stroke="none"/>
  <circle cx="260" cy="252" r="6" fill="currentColor" stroke="none"/>
  <text x="170" y="76" text-anchor="middle" fill="currentColor" stroke="none" font-size="18">1</text>
  <text x="64" y="270" text-anchor="middle" fill="currentColor" stroke="none" font-size="18">2</text>
  <text x="276" y="270" text-anchor="middle" fill="currentColor" stroke="none" font-size="18">3</text>
  <text x="100" y="166" text-anchor="middle" fill="currentColor" stroke="none" font-size="15">30&#176;</text>
  <text x="240" y="166" text-anchor="middle" fill="currentColor" stroke="none" font-size="15">30&#176;</text>
  <text x="170" y="278" text-anchor="middle" fill="currentColor" stroke="none" font-size="15">30&#176;</text>
  <text x="170" y="316" text-anchor="middle" fill="currentColor" stroke="none" font-size="17">30&#176; + 30&#176; + 30&#176; = 90&#176;</text>
  <text x="170" y="338" text-anchor="middle" fill="currentColor" stroke="none" font-size="13.5" fill-opacity="0.72">the loop does not close</text>

  <text x="510" y="32" text-anchor="middle" fill="currentColor" stroke="none" font-size="16.5">Angles 0&#176;, &#8722;30&#176;, &#8722;60&#176;</text>
  <polyline points="510,92 465,172 420,252" marker-mid="url(#lrn-ah)"/>
  <polyline points="420,252 510,252 600,252" marker-mid="url(#lrn-ah)"/>
  <polyline points="600,252 555,172 510,92" marker-mid="url(#lrn-ah)"/>
  <circle cx="510" cy="92" r="6" fill="currentColor" stroke="none"/>
  <circle cx="420" cy="252" r="6" fill="currentColor" stroke="none"/>
  <circle cx="600" cy="252" r="6" fill="currentColor" stroke="none"/>
  <text x="510" y="76" text-anchor="middle" fill="currentColor" stroke="none" font-size="18">1</text>
  <text x="404" y="270" text-anchor="middle" fill="currentColor" stroke="none" font-size="18">2</text>
  <text x="616" y="270" text-anchor="middle" fill="currentColor" stroke="none" font-size="18">3</text>
  <text x="440" y="166" text-anchor="middle" fill="currentColor" stroke="none" font-size="15">30&#176;</text>
  <text x="582" y="166" text-anchor="middle" fill="currentColor" stroke="none" font-size="15">&#8722;60&#176;</text>
  <text x="510" y="278" text-anchor="middle" fill="currentColor" stroke="none" font-size="15">30&#176;</text>
  <text x="510" y="316" text-anchor="middle" fill="currentColor" stroke="none" font-size="17">30&#176; + 30&#176; &#8722; 60&#176; = 0&#176;</text>
  <text x="510" y="338" text-anchor="middle" fill="currentColor" stroke="none" font-size="13.5" fill-opacity="0.72">the loop closes</text>
</svg>
<figcaption style="font-size:0.9em;opacity:0.72;margin-top:0.5rem">
Both pictures show the same three buses and branches. On the left every branch carries a 30&#176; difference and the loop adds up to 90&#176;, which no bus angles can produce; on the right the genuine angles 0&#176;, &#8722;30&#176;, &#8722;60&#176; give differences 30&#176;, 30&#176;, &#8722;60&#176; that close the loop at 0&#176;. SOCWR checks each branch on its own and never sees the left case is impossible, while the SDP 3&#215;3 minor looks at all three together and rules it out.
</figcaption>
</figure>
```

!!! example "The lie SOCWR cannot see"
    All three buses have magnitude ``v = 1``. Suppose the optimizer picks the same blocks on
    every branch:

    | branch | ``W_r`` | ``W_i`` | implicit angle ``\theta_{ij} = \operatorname{atan2}(W_i, W_r)`` |
    |---|---|---|---|
    | 1-2 | ``\cos 30^\circ = 0.866`` | ``\sin 30^\circ = 0.5`` | ``30^\circ`` |
    | 2-3 | ``0.866`` | ``0.5`` | ``30^\circ`` |
    | 3-1 | ``0.866`` | ``0.5`` | ``30^\circ`` |

    Test each branch against the SOCWR cone:
    ``W_r^2 + W_i^2 = 0.866^2 + 0.5^2 = 1 = W_{ii}W_{jj}``. All three pass, with perfect
    equality (they look like genuine voltages!). But sum the angles around the loop. For any
    genuine bus angles the differences telescope to zero,
    ``(\theta_1 - \theta_2) + (\theta_2 - \theta_3) + (\theta_3 - \theta_1) = 0``; here they
    would sum to ``30^\circ + 30^\circ + 30^\circ = 90^\circ \ne 0``. No three bus angles
    produce three ``30^\circ`` differences around a triangle. SOCWR accepted a point that no
    physical network realizes, and no per-branch squeeze can bar it, because each branch in
    isolation is impeccable.

    Hand the same blocks to the SDP and it refuses them. With every magnitude ``1`` the three
    cross-entries are ``W_{12} = W_{23} = W_{31} = e^{j30^\circ}``, so the ``3\times 3``
    matrix is exactly the unit-modulus triangle of rung 3 with
    ``\alpha = \beta = \gamma = 30^\circ`` and ``\xi = 90^\circ``. Its determinant is
    ``2(\cos 90^\circ - 1) = -2``, which is negative, so the matrix is not positive
    semidefinite. The SDP throws out the very point SOCWR waved through.

!!! example "The same triangle, genuine voltages"
    Take genuine voltages ``V_1 = 1``, ``V_2 = e^{-j30^\circ}``, ``V_3 = e^{-j60^\circ}``
    (bus angles ``0^\circ``, ``-30^\circ``, ``-60^\circ``, so the differences are
    ``30^\circ`` on branches 1-2 and 2-3 and ``-60^\circ`` on 3-1, closing the loop):

    ```math
    V = \begin{bmatrix} 1 \\ e^{-j30^\circ} \\ e^{-j60^\circ} \end{bmatrix},
    \qquad
    V^{*} = \begin{bmatrix} 1 & e^{j30^\circ} & e^{j60^\circ} \end{bmatrix},
    ```

    ```math
    W = V\,V^{*}
    = \begin{bmatrix}
    V_1 V_1^{*} & V_1 V_2^{*} & V_1 V_3^{*} \\
    V_2 V_1^{*} & V_2 V_2^{*} & V_2 V_3^{*} \\
    V_3 V_1^{*} & V_3 V_2^{*} & V_3 V_3^{*}
    \end{bmatrix}
    = \begin{bmatrix}
    1 & e^{j30^\circ} & e^{j60^\circ} \\
    e^{-j30^\circ} & 1 & e^{j30^\circ} \\
    e^{-j60^\circ} & e^{-j30^\circ} & 1
    \end{bmatrix}.
    ```

    Because this ``W`` is built as ``V V^{*}`` it is rank one, so the loop closes on its own.
    This is the honest counterpart of the lie above: the same triangle, but with angles that
    come from three real bus phasors.

!!! example "A fully worked three-bus instance"
    Data: ``V_1 = 1.00\angle 0^\circ``, ``V_2 = 0.95\angle{-10^\circ}``,
    ``V_3 = 1.05\angle{-25^\circ}``. Angle differences and magnitude products:

    | pair | ``v_i v_j`` | ``\theta_{ij} = \theta_i - \theta_j`` | ``\cos\theta_{ij}`` | ``\sin\theta_{ij}`` |
    |---|---|---|---|---|
    | 1-2 | ``1.00 \times 0.95 = 0.95`` | ``0 - (-10) = 10^\circ`` | ``0.98481`` | ``0.17365`` |
    | 1-3 | ``1.00 \times 1.05 = 1.05`` | ``0 - (-25) = 25^\circ`` | ``0.90631`` | ``0.42262`` |
    | 2-3 | ``0.95 \times 1.05 = 0.9975`` | ``-10 - (-25) = 15^\circ`` | ``0.96593`` | ``0.25882`` |

    (The cycle checks out: ``10^\circ + 15^\circ + (-25^\circ) = 0``.)

    Build ``W^{\mathrm{re}}`` (with ``W^{\mathrm{re}}_{ij} = v_i v_j \cos\theta_{ij}`` and
    ``W^{\mathrm{re}}_{ii} = v_i^2``) and ``W^{\mathrm{im}}`` (with
    ``W^{\mathrm{im}}_{ij} = v_i v_j \sin\theta_{ij}`` and
    ``W^{\mathrm{im}}_{ii} = \sin 0 = 0``):

    ```math
    W^{\mathrm{re}} = \begin{bmatrix}
    1.00000 & 0.93557 & 0.95162 \\
    0.93557 & 0.90250 & 0.96351 \\
    0.95162 & 0.96351 & 1.10250
    \end{bmatrix} \ (\text{symmetric}),
    \qquad
    W^{\mathrm{im}} = \begin{bmatrix}
    0 & 0.16497 & 0.44375 \\
    -0.16497 & 0 & 0.25817 \\
    -0.44375 & -0.25817 & 0
    \end{bmatrix} \ (\text{antisymmetric}).
    ```

    The two together give the complex ``W = W^{\mathrm{re}} + j\,W^{\mathrm{im}}``, entry by
    entry:

    ```math
    W = \begin{bmatrix}
    1.00000 & 0.93557 + j\,0.16497 & 0.95162 + j\,0.44375 \\
    0.93557 - j\,0.16497 & 0.90250 & 0.96351 + j\,0.25817 \\
    0.95162 - j\,0.44375 & 0.96351 - j\,0.25817 & 1.10250
    \end{bmatrix}.
    ```

    Finally the real ``6 \times 6`` matrix that goes to the solver,
    ``M = \begin{bmatrix} W^{\mathrm{re}} & W^{\mathrm{im}} \\ -W^{\mathrm{im}} & W^{\mathrm{re}} \end{bmatrix}``:

    ```math
    M = \left(\begin{array}{ccc|ccc}
    1.00000 & 0.93557 & 0.95162 & 0 & 0.16497 & 0.44375 \\
    0.93557 & 0.90250 & 0.96351 & -0.16497 & 0 & 0.25817 \\
    0.95162 & 0.96351 & 1.10250 & -0.44375 & -0.25817 & 0 \\ \hline
    0 & -0.16497 & -0.44375 & 1.00000 & 0.93557 & 0.95162 \\
    0.16497 & 0 & -0.25817 & 0.93557 & 0.90250 & 0.96351 \\
    0.44375 & 0.25817 & 0 & 0.95162 & 0.96351 & 1.10250
    \end{array}\right).
    ```

    The bottom-left block is ``W^{\mathrm{im}}`` with every sign flipped which, since
    ``W^{\mathrm{im}}`` is already antisymmetric, equals ``(W^{\mathrm{im}})^{\top}``. The
    whole ``M`` comes out symmetric, the format the SDP solver accepts.

## See also

- Neighbours: [SOCWR (Jabr conic relaxation)](@ref) · [QC (Quadratic Convex)](@ref)
- Exact reference: [AC Polar](@ref); the ``W``-space blocks come from the
  [AC Rectangular (ACR)](@ref) model.
- Reference implementation: struct `SDPWRM`, source `src/modeling/branch/sdpwrm.jl`;
  PowerModels equivalent `SDPWRMPowerModel`.
