# ============================================================
# Zone connectivity repair (pure graph, no Metis dependency)
# ============================================================

"""
    _zone_connected_components(zone, adj) -> Dict{Int,Vector{Vector{Int}}}

Group nodes by zone and, within each zone, into connected components of
the subgraph induced by that zone. Each node is visited once (O(V+E)).
Returns `zone_id => [component, ...]`, each component a vector of node
indices.
"""
function _zone_connected_components(
    zone::Vector{Int},
    adj::Vector{Vector{Int}},
)
    n = length(zone)
    seen = falses(n)
    comps_by_zone = Dict{Int,Vector{Vector{Int}}}()
    stack = Int[]
    for s in 1:n
        seen[s] && continue
        z = zone[s]
        comp = Int[]
        empty!(stack)
        push!(stack, s)
        seen[s] = true
        while !isempty(stack)
            u = pop!(stack)
            push!(comp, u)
            for v in adj[u]
                if !seen[v] && zone[v] == z
                    seen[v] = true
                    push!(stack, v)
                end
            end
        end
        push!(get!(comps_by_zone, z, Vector{Vector{Int}}()), comp)
    end
    return comps_by_zone
end

"""
    all_zones_connected(zone, adj) -> Bool

True iff every zone forms a single connected component on the graph.
"""
function all_zones_connected(zone::Vector{Int}, adj::Vector{Vector{Int}})
    for (_, comps) in _zone_connected_components(zone, adj)
        length(comps) <= 1 || return false
    end
    return true
end

"""
    enforce_zone_connectivity!(zone, adj) -> zone

Make every zone connected by a component-based, iterated repair:

- within each zone, keep the LARGEST component as the zone core;
- reassign every smaller component, as a whole, to the neighboring zone
  it has the most edges to (ties broken by smallest zone id);
- a component with no edge to any other zone is a genuine electrical
  island of the whole graph and is left in place;
- iterate to a fixed point.

Each whole-component move strictly decreases the global number of
connected components (the moved component merges into an adjacent
component of the target zone), so the loop terminates. The zone core is
never moved, so no zone is emptied and `k` is preserved. After this
call, `all_zones_connected` holds unless the underlying graph itself is
disconnected into more islands than there are zones to hold them.
"""
function enforce_zone_connectivity!(zone::Vector{Int}, adj::Vector{Vector{Int}})
    n = length(zone)
    # Global component count is >= number of zones and <= n, and drops by
    # >=1 per move, so n passes is a safe upper bound.
    for _pass in 1:(n + 1)
        comps_by_zone = _zone_connected_components(zone, adj)
        changed = false
        for (z, comps) in comps_by_zone
            length(comps) <= 1 && continue
            core_idx = argmax(map(length, comps))   # first largest = core
            for (ci, comp) in enumerate(comps)
                ci == core_idx && continue
                votes = Dict{Int,Int}()
                for u in comp
                    for v in adj[u]
                        zv = zone[v]
                        zv == z && continue
                        votes[zv] = get(votes, zv, 0) + 1
                    end
                end
                isempty(votes) && continue   # genuine island: cannot attach
                best = 0
                bestv = -1
                for zn in sort!(collect(keys(votes)))  # smallest zone id on ties
                    if votes[zn] > bestv
                        bestv = votes[zn]
                        best = zn
                    end
                end
                for u in comp
                    zone[u] = best
                end
                changed = true
            end
        end
        changed || break
    end
    return zone
end
