@inline _wrap_angle_diff(Δ::Float64) = atan(sin(Δ), cos(Δ))

@inline _elapsed_s(t0_ns::Integer) = (time_ns() - t0_ns) / 1e9

@inline function _zone_cut_edges(part, interzonal_branches::AbstractVector{Int})
    m = part.branch_idx_full_to_zone
    return [e for e in interzonal_branches if m[e] != 0]
end

@inline function _admm_totals(partitions)
    parts = partitions.parts
    total = 0.0
    @inbounds @simd for k in eachindex(parts)
        total += parts[k].generation_cost
    end
    return total
end

@inline stop_if_converged(ws::AbstractADMMWorkspace) =
    (ws.r_primal < ws.eps_primal) && (ws.r_dual < ws.eps_dual)

@inline function _setup_partition_model!(
    part::AbstractPartition,
    sys_zone::SystemData,
    ::Type{M},
    optimizer_factory,
) where {M <: OPFModel}
    prepare_hotpaths!(part)
    model = build_model(sys_zone, M; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    return model
end
