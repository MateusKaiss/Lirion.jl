@inline function _total_generation_cost(partitions)
    parts = partitions.parts
    total = 0.0
    @inbounds @simd for k in eachindex(parts)
        total += parts[k].generation_cost
    end
    return total
end
