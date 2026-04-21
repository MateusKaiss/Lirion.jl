# ============================================================
# paper_results.jl
# ============================================================
# Top-level driver that runs the three stages in sequence:
#   1. paper_results_centralized.jl   → cases/_centralized.toml
#   2. paper_results_admm.jl          → cases/_admm.toml
#   3. paper_results_aladin.jl        → cases/_aladin.toml
#
# Each stage is idempotent and writes its own TOML atomically after
# every solve. Re-running resumes from where the previous run stopped.
#
# Usage (after `using Lirion`):
#
#     julia> include(joinpath(pkgdir(Lirion), "paper", "paper_results.jl"))
#
# To run just one stage, include the corresponding file directly.
# ============================================================

const _PAPER_DIR = @__DIR__

include(joinpath(_PAPER_DIR, "paper_results_centralized.jl"))
include(joinpath(_PAPER_DIR, "paper_results_admm.jl"))
include(joinpath(_PAPER_DIR, "paper_results_aladin.jl"))
