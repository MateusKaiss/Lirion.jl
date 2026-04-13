# =======================================================================
# ALADIN-DC — Per-variant workspace types (stubs for Milestone 1b)
# =======================================================================

# ---------------------------------------------------------------
# Partition type (DC under ALADIN)
# ---------------------------------------------------------------

mutable struct PartitionDCAladin <: AbstractPartition
    # Placeholder — real fields added in 1c.
    PartitionDCAladin() = new()
end

# ---------------------------------------------------------------
# Variant-specific workspaces
# ---------------------------------------------------------------

mutable struct DCStateWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Any
    iter::Int
    r_primal::Float64
    r_dual::Float64
    eps_primal::Float64
    eps_dual::Float64

    DCStateWorkspaceALADIN() = new(nothing, 0, Inf, Inf, 0.0, 0.0)
end

mutable struct DCPowerWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Any
    iter::Int
    r_primal::Float64
    r_dual::Float64
    eps_primal::Float64
    eps_dual::Float64

    DCPowerWorkspaceALADIN() = new(nothing, 0, Inf, Inf, 0.0, 0.0)
end

mutable struct DCHybridWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Any
    iter::Int
    r_primal::Float64
    r_dual::Float64
    eps_primal::Float64
    eps_dual::Float64

    DCHybridWorkspaceALADIN() = new(nothing, 0, Inf, Inf, 0.0, 0.0)
end

# ---------------------------------------------------------------
# Workspace builder — dispatches on ALADIN{Variant, DC}
# ---------------------------------------------------------------

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    tag::ALADIN{State, DC},
)
    return error("ALADIN{State,DC} workspace not yet implemented (Milestone 1c)")
end

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    tag::ALADIN{Power, DC},
)
    return error("ALADIN{Power,DC} workspace not yet implemented (Milestone 1d)")
end

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    tag::ALADIN{Hybrid, DC},
)
    return error("ALADIN{Hybrid,DC} workspace not yet implemented (Milestone 1d)")
end
