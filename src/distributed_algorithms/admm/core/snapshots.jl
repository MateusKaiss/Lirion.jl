# ============================================================
# Buffer preparation
# ============================================================

@inline _bus_dims(part) = (length(part.bus.id), length(part.zone_shared_bus_ids))

# ------------------------------------------------------------
# DC
# ------------------------------------------------------------

function prepare_hotpaths!(part::PartitionDC)
    n_bus, n_shared = _bus_dims(part)
    part.theta_buf = zeros(n_bus)
    part.theta_shared_buf = zeros(n_shared)
    return nothing
end

# ------------------------------------------------------------
# AC
# ------------------------------------------------------------

function prepare_hotpaths!(part::PartitionAC)
    n_bus, n_shared = _bus_dims(part)
    part.theta_buf = zeros(n_bus)
    part.v_mag_buf = ones(n_bus)
    part.theta_shared_buf = zeros(n_shared)
    part.v_mag_shared_buf = ones(n_shared)
    return nothing
end

# ------------------------------------------------------------
# SOCWR
# ------------------------------------------------------------

function prepare_hotpaths!(part::PartitionSOCWR)
    n_bus, n_shared = _bus_dims(part)
    part.v2_buf = zeros(n_bus)
    part.v2_shared_buf = zeros(n_shared)
    return nothing
end

# ============================================================
# JuMP value helpers
# ============================================================

@inline function _val_or_start(
    model::JuMP.Model,
    x::JuMP.VariableRef,
    default::Float64,
    has_vals::Bool,
)
    if has_vals
        return JuMP.value(x)
    end
    sv = JuMP.start_value(x)
    return sv === nothing ? default : sv
end

@inline _val_at(
    model::JuMP.Model,
    bus_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2},
    b::Int,
    t::Int,
    default::Float64,
    has_vals::Bool,
) = _val_or_start(model, bus_var[b, t], default, has_vals)

# ============================================================
# Snapshot kernels
# ============================================================

@inline function _snapshot_bus_var!(
    buf::Vector{Float64},
    model::JuMP.Model,
    bus_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2},
    bus_ids::Vector{Int},
    default::Float64,
)
    has_vals = JuMP.has_values(model)
    @inbounds for i in eachindex(bus_ids)
        buf[i] = _val_at(model, bus_var, bus_ids[i], 1, default, has_vals)
    end
    return buf
end

@inline function _snapshot_shared_var!(
    buf::Vector{Float64},
    model::JuMP.Model,
    bus_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2},
    shared_bus_ids::Vector{Int},
    default::Float64,
)
    has_vals = JuMP.has_values(model)
    @inbounds for i in eachindex(shared_bus_ids)
        buf[i] = _val_at(
            model,
            bus_var,
            shared_bus_ids[i],
            1,
            default,
            has_vals,
        )
    end
    return buf
end

# ============================================================
# Full-bus snapshots
# ============================================================

# ------------------------------------------------------------
# DC
# ------------------------------------------------------------

@inline snapshot_theta!(part::PartitionDC) =
    _snapshot_bus_var!(part.theta_buf, part.model, part.theta_var, part.bus.id, 0.0)

# ------------------------------------------------------------
# AC
# ------------------------------------------------------------

@inline snapshot_theta!(part::PartitionAC) =
    _snapshot_bus_var!(part.theta_buf, part.model, part.theta_var, part.bus.id, 0.0)

@inline snapshot_v_mag!(part::PartitionAC) =
    _snapshot_bus_var!(part.v_mag_buf, part.model, part.v_mag_var, part.bus.id, 1.0)

# ------------------------------------------------------------
# SOCWR
# ------------------------------------------------------------

@inline snapshot_v2!(part::PartitionSOCWR) =
    _snapshot_bus_var!(part.v2_buf, part.model, part.v2_var, part.bus.id, 1.0)

# ============================================================
# Shared-bus snapshots
# ============================================================

# ------------------------------------------------------------
# DC
# ------------------------------------------------------------

@inline snapshot_theta_shared!(part::PartitionDC) =
    _snapshot_shared_var!(
        part.theta_shared_buf, part.model, part.theta_var,
        part.zone_shared_bus_ids, 0.0,
    )

# ------------------------------------------------------------
# AC
# ------------------------------------------------------------

@inline snapshot_theta_shared!(part::PartitionAC) =
    _snapshot_shared_var!(
        part.theta_shared_buf, part.model, part.theta_var,
        part.zone_shared_bus_ids, 0.0,
    )

@inline snapshot_v_mag_shared!(part::PartitionAC) =
    _snapshot_shared_var!(
        part.v_mag_shared_buf, part.model, part.v_mag_var,
        part.zone_shared_bus_ids, 1.0,
    )

# ------------------------------------------------------------
# SOCWR
# ------------------------------------------------------------

@inline snapshot_v2_shared!(part::PartitionSOCWR) =
    _snapshot_shared_var!(
        part.v2_shared_buf, part.model, part.v2_var,
        part.zone_shared_bus_ids, 1.0,
    )
