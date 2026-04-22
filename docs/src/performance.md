# Performance & Threading

Lirion runs two very different workloads and each one has its own
optimal threading configuration. Getting this right can be the
difference between a distributed run taking **minutes** or **hours**.

When Lirion loads it prints the current setup:

```
[ Info: [Lirion] Threads: Julia=1 | BLAS=14
```

- **Julia threads**: used to solve zone subproblems in parallel
  (distributed algorithms only).
- **BLAS threads**: used by every solver for dense linear algebra
  (matrix factorization, Schur complement, etc.).


## The golden rule

**Never nest parallelism.** Either go *deep* (one solver, many BLAS
threads) **or** *wide* (many zones, one BLAS thread each). Doing both
(e.g., 28 zones × 14 BLAS threads on a 28-core machine) creates
thread contention and cache thrashing — you end up slower than
sequential.


## Centralized

A single solver instance; Julia does not parallelize at the zone
level. Give BLAS all the physical cores so matrix factorization is
fast.

### Recommended

- `JULIA_NUM_THREADS` = 1 (default)
- `BLAS` = physical core count (Julia's default on most systems)

### Why not hyperthreads for BLAS?

Hyperthreads share the floating-point/vector units of a physical
core. BLAS is FPU-bound, so two hyperthreads competing for the same
FPU give no speedup — only scheduling overhead. Julia's BLAS default
correctly picks the physical count (e.g. 14 on a 14c/28t CPU).

### Command

No environment variables needed — the default is already correct.
If you want to be explicit:

**Linux / macOS (bash / zsh):**

```bash
unset JULIA_NUM_THREADS
# (BLAS default is already physical cores)
```

**Windows (PowerShell):**

```powershell
Remove-Item Env:JULIA_NUM_THREADS -ErrorAction SilentlyContinue
```


## Distributed (ADMM / ALADIN)

Lirion solves N zone subproblems in parallel via Julia threads.
Each subproblem is small, so we want zone-level parallelism, not
BLAS-level. Pin BLAS / OMP / MKL to a single thread so they do not
spawn extra threads inside each zone solve.

### Recommended

- `JULIA_NUM_THREADS` = number of logical cores (one thread per zone)
- `OPENBLAS_NUM_THREADS` = 1
- `OMP_NUM_THREADS` = 1
- `MKL_NUM_THREADS` = 1

### Command

Set these **before** launching Julia.

**Linux (bash / zsh):**

```bash
export JULIA_NUM_THREADS=$(getconf _NPROCESSORS_ONLN)
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
julia --project
```

**macOS (bash / zsh):**

```bash
export JULIA_NUM_THREADS=$(sysctl -n hw.ncpu)
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
julia --project
```

**Windows (PowerShell):**

```powershell
$env:JULIA_NUM_THREADS   = $env:NUMBER_OF_PROCESSORS
$env:OPENBLAS_NUM_THREADS = "1"
$env:OMP_NUM_THREADS      = "1"
$env:MKL_NUM_THREADS      = "1"
julia --project
```

**Windows (cmd.exe):**

```cmd
set JULIA_NUM_THREADS=%NUMBER_OF_PROCESSORS%
set OPENBLAS_NUM_THREADS=1
set OMP_NUM_THREADS=1
set MKL_NUM_THREADS=1
julia --project
```

### Shortcut

`JULIA_NUM_THREADS` can also be set with the `-t` flag. `julia -t
auto` picks a sensible default based on the machine, cross-platform.
You still need the `*_NUM_THREADS=1` exports for BLAS / OMP / MKL.


## Verifying the configuration

After loading Lirion, the info line tells you what the current
session is using:

```julia
julia> using Lirion
[ Info: [Lirion] Threads: Julia=28 | BLAS=1
```

- `Julia=28 | BLAS=1` → correctly configured for distributed.
- `Julia=1  | BLAS=14` → correctly configured for centralized.

If Julia shows `BLAS=14` when you wanted `BLAS=1`, your
`OPENBLAS_NUM_THREADS=1` was not exported before launching Julia.


## What Lirion does for you

For distributed runs, Lirion wraps the solve call in a scoped helper
that **forces BLAS to 1** for the duration of the distributed
algorithm and **restores your previous value afterwards**. This is a
safety net — even if you forget to set `OPENBLAS_NUM_THREADS=1`, the
distributed run will not suffer from nested BLAS parallelism.

Lirion cannot set `JULIA_NUM_THREADS` for you — that is fixed at the
moment Julia starts. If you want many zones in parallel, you must
launch Julia with the environment variable or the `-t` flag set.
