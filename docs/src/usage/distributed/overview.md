# Running a distributed solve

A distributed solve uses the same `solve` you already know. You keep the model and the case, and you swap the `algorithm`: instead of `centralized()`, you pass `admm(...)` or `aladin(...)`. Lirion then splits the network into zones, each zone solves its own local OPF, and the zones coordinate on the boundary until they agree.

The distributed algorithms run on the DC, AC, and SOCWR models.

!!! important "Important"
    The runs below are simple examples. The way you configure a distributed algorithm changes the optimization problem itself, so reading the next pages first is strongly recommended. See [Coupling and sectioning](coupling_sectioning.md), [ADMM settings](@ref), and [ALADIN settings](@ref). For the theory behind the algorithms, see [ADMM](@ref) and [ALADIN](@ref).

## With ADMM

`admm()` defaults to `coupling = :hybrid` and `sectioning = :tie_line`.

```julia
using Lirion

result = solve("pglib_opf_case57_ieee.m";
    model = ac(),
    algorithm = admm(rho = 1e3),
)
```

## With ALADIN

`aladin()` defaults to `coupling = :state` and `sectioning = :tie_line`.

```julia
using Lirion

result = solve("pglib_opf_case57_ieee.m";
    model = dc(),
    algorithm = aladin(),
)
```

## Reading the result

A distributed run returns the same `LirionResult` as a centralized one, plus the per-zone tables and the convergence history. For the extra fields and how to read them, see the Decentralized section of [Results](@ref).

## Where to go next

- [Coupling and sectioning](coupling_sectioning.md) — how the network is split and what the zones agree on.
- [ADMM settings](@ref) — penalties, adaptive `rho`, and stopping tolerances.
- [ALADIN settings](@ref) — proximal mode, stopping criterion, and the coordination step.
