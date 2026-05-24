# ============================================================
# SystemData
# ============================================================

struct SystemData
    baseMVA::Float64
    T::Int

    bus::BusData
    branch::BranchData
    gen::GenData
    storage::StorageData

    bus_ts::BusTimeseriesData
    gencost_ts::GencostTimeseriesData

    powerflow::PowerFlowData
    objective::ObjectiveData

    case_raw::Dict{String, Any}
end

# ============================================================
# Distributed system
# ============================================================

struct DistributedSystem
    zones::Vector{Int}
    parts::Vector{SystemData}

    zone_branch_idx_in_full::Vector{Vector{Int}}
    branch_idx_full_to_zone::Vector{Vector{Int}}
    zone_shared_bus_ids::Vector{Vector{Int}}

    full_shared_bus_ids::Vector{Int}
    interzonal_branches::Vector{Int}
    branch_pairs::Vector{NTuple{4, Int}}
    bus_pairs::Vector{NTuple{3, Int}}
end

# ============================================================
# Utils
# ============================================================

function _raw_snapshot(case_raw::CaseRaw)
    d = Dict{String, Any}()
    d["name"] = case_raw.name
    d["baseMVA"] = case_raw.baseMVA
    d["T"] = case_raw.T
    d["maps"] = case_raw.maps
    return d
end

# ============================================================
# Builders
# ============================================================

function build_system_data(case_raw::CaseRaw)::SystemData
    bus_data = build_bus_data(case_raw)
    branch_data = build_branch_data(case_raw)
    gen_data = build_gen_data(case_raw)
    storage_data = build_storage_data(case_raw)

    bus_ts_data = build_bus_timeseries_data(case_raw)
    gencost_ts_data = build_gencost_timeseries_data(case_raw)

    powerflow_data = build_power_flow_data(case_raw)
    objective_data = build_objective_data(case_raw)

    return SystemData(
        case_raw.baseMVA,
        case_raw.T,
        bus_data,
        branch_data,
        gen_data,
        storage_data,
        bus_ts_data,
        gencost_ts_data,
        powerflow_data,
        objective_data,
        _raw_snapshot(case_raw),
    )
end

function build_distributed_system(
    case_raw::CaseRaw,
    sectioning::Sectioning,
)::DistributedSystem
    partitioning = prepare_distributed_partitions(case_raw, sectioning)
    interzonal_branches = _find_interzonal_branches(case_raw, sectioning)

    zones = partitioning.zone_ids
    parts_raw = partitioning.parts
    n_zones = length(zones)

    parts = Vector{SystemData}(undef, n_zones)
    zone_branch_idx_in_full = Vector{Vector{Int}}(undef, n_zones)
    branch_idx_full_to_zone = Vector{Vector{Int}}(undef, n_zones)
    zone_shared_bus_ids = Vector{Vector{Int}}(undef, n_zones)

    @inbounds for i in 1:n_zones
        zone_raw = parts_raw[i]

        parts[i] = build_system_data(zone_raw)

        zone_branch_idx_in_full[i] = zone_raw.dist.zone_branch_idx_in_full
        branch_idx_full_to_zone[i] = zone_raw.dist.branch_idx_full_to_zone
        zone_shared_bus_ids[i] = zone_raw.dist.zone_shared_bus_ids
    end

    full_shared_bus_ids = case_raw.dist.full_shared_bus_ids
    branch_pairs = case_raw.dist.branch_pairs
    bus_pairs = case_raw.dist.bus_pairs

    return DistributedSystem(
        zones,
        parts,
        zone_branch_idx_in_full,
        branch_idx_full_to_zone,
        zone_shared_bus_ids,
        full_shared_bus_ids,
        interzonal_branches,
        branch_pairs,
        bus_pairs,
    )
end
