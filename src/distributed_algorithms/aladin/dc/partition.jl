# ============================================================
# ALADIN-DC — Partition initialization
# ------------------------------------------------------------
# Phase 1 of canonical ALADIN uses the SAME JuMP model as the ADMM
# local subproblem (the algorithmic difference is purely in the
# objective, which is rebuilt every iteration by update_objective!).
# So we delegate partition setup to the ADMM init helpers and let
# build_aladin_zone_data extract the post-build sensitivity layout.
# ============================================================

@inline init_state_partition_aladin!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    rho::Float64,
) = init_state_partition!(partitions, sys, optimizer_factory, rho)

@inline init_power_partition_aladin!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    rho::Float64,
    interzonal_branches::Vector{Int},
) = init_power_partition!(partitions, sys, optimizer_factory, rho, interzonal_branches)

@inline init_hybrid_partition_aladin!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_theta::Float64,
    rho_p::Float64,
    interzonal_branches::Vector{Int},
) = init_hybrid_partition!(
    partitions,
    sys,
    optimizer_factory,
    rho_theta,
    rho_p,
    interzonal_branches,
)
