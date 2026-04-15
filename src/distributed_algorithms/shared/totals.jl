@inline function _total_generation_cost(partitions)
    parts = partitions.parts
    total = 0.0
    @inbounds @simd for k in eachindex(parts)
        total += parts[k].generation_cost
    end
    return total
end

@inline function _finalize(partitions, t0_ns::Integer, iter_done::Integer)
    total_cost = _total_generation_cost(partitions)
    elapsed_s = _elapsed_s(t0_ns)
    return (;
        total_cost,
        partitions,
        iter_done = Int(iter_done),
        elapsed_s,
    )
end
