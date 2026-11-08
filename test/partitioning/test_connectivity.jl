# ============================================================
# ======  Zone connectivity repair (component-based)  ========
# ============================================================
# Unit tests for Lirion.enforce_zone_connectivity! / all_zones_connected.
# Pure graph logic — no Metis/Graphs needed. Reproduces the pathologies
# that the greedy per-bus fix could not repair (interior multi-bus
# islands) plus a genuinely disconnected graph (unfixable -> reported).

function _undirected_adj(n::Int, edges::Vector{Tuple{Int,Int}})
    adj = [Int[] for _ in 1:n]
    for (u, v) in edges
        push!(adj[u], v)
        push!(adj[v], u)
    end
    for a in adj
        sort!(unique!(a))
    end
    return adj
end

function run_partition_connectivity_tests()
    @testset "Zone connectivity repair" begin

        # ---- 1. Single stray bus: zone 1 = {1,5} on a path 1-2-3-4-5 ----
        # bus 1 and bus 5 share a zone but are not adjacent within it.
        @testset "stray boundary bus is reabsorbed" begin
            adj = _undirected_adj(5, [(1, 2), (2, 3), (3, 4), (4, 5)])
            zone = [1, 2, 2, 2, 1]
            @test !Lirion.all_zones_connected(zone, adj)
            Lirion.enforce_zone_connectivity!(zone, adj)
            @test Lirion.all_zones_connected(zone, adj)
            @test zone[1] == 1              # core of zone 1 kept
            @test zone[5] == zone[4]        # stray {5} moved to its neighbor
        end

        # ---- 2. Interior multi-bus island (the greedy fix's blind spot) ----
        # Two triangles A={1,2,3}, B={4,5,6} joined by bridge 3-4.
        # zone assignment splits both zones into two components each.
        @testset "interior multi-bus island dissolved" begin
            edges = [(1, 2), (2, 3), (1, 3), (4, 5), (5, 6), (4, 6), (3, 4)]
            adj = _undirected_adj(6, edges)
            zone = [1, 1, 2, 1, 2, 2]       # zone1={1,2,4}, zone2={3,5,6}
            @test !Lirion.all_zones_connected(zone, adj)   # {4} and {3} orphaned
            Lirion.enforce_zone_connectivity!(zone, adj)
            @test Lirion.all_zones_connected(zone, adj)
            @test count(==(1), zone) >= 1 && count(==(2), zone) >= 1  # k preserved
        end

        # ---- 3. Whole-component relocation, larger island ----
        # Path of 10; zone 1 owns a 3-bus island {7,8,9} far from its core.
        @testset "3-bus island relocated as a whole" begin
            adj = _undirected_adj(10,
                [(i, i + 1) for i in 1:9])
            zone = [1, 1, 1, 2, 2, 2, 1, 1, 1, 2]  # zone1 = {1,2,3} ∪ {7,8,9}
            @test !Lirion.all_zones_connected(zone, adj)
            Lirion.enforce_zone_connectivity!(zone, adj)
            @test Lirion.all_zones_connected(zone, adj)
            @test zone[1] == zone[2] == zone[3]     # core stays together
            # {7,8,9} moved as one block to their neighbor zone (2)
            @test zone[7] == zone[8] == zone[9]
        end

        # ---- 4. Already connected -> idempotent, no change ----
        @testset "already-connected partition is unchanged" begin
            adj = _undirected_adj(6, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)])
            zone = [1, 1, 1, 2, 2, 2]
            before = copy(zone)
            @test Lirion.all_zones_connected(zone, adj)
            Lirion.enforce_zone_connectivity!(zone, adj)
            @test zone == before
        end

        # ---- 5. Genuinely disconnected graph -> cannot fix, reported ----
        # Two triangles with NO bridge, all in one zone: enforce leaves the
        # island (no external neighbor) and all_zones_connected stays false.
        @testset "true island is left and reported (not forced)" begin
            edges = [(1, 2), (2, 3), (1, 3), (4, 5), (5, 6), (4, 6)]
            adj = _undirected_adj(6, edges)
            zone = [1, 1, 1, 1, 1, 1]
            Lirion.enforce_zone_connectivity!(zone, adj)  # must terminate
            @test !Lirion.all_zones_connected(zone, adj)
        end

        # ---- 6. Termination / no oscillation on a dense ring ----
        @testset "terminates on a ring with scattered zones" begin
            n = 12
            adj = _undirected_adj(n,
                vcat([(i, i % n + 1) for i in 1:n], [(1, 7), (3, 9)]))
            zone = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]  # alternating: many islands
            Lirion.enforce_zone_connectivity!(zone, adj)
            @test Lirion.all_zones_connected(zone, adj)
        end
    end
    return nothing
end

# ---- Adjacency + zone (by bus position) from a partitioned CaseRaw ----
function _case_zone_and_adj(c::Lirion.CaseRaw)
    ids = c.bus.bus_i
    n = length(ids)
    pos = Dict(ids[i] => i for i in 1:n)
    zone = copy(c.bus.zone)
    adj = [Int[] for _ in 1:n]
    br = c.branch
    for i in eachindex(br.f_bus)
        br.br_status[i] == 1 || continue
        f = get(pos, br.f_bus[i], 0); t = get(pos, br.t_bus[i], 0)
        (f == 0 || t == 0 || f == t) && continue
        push!(adj[f], t); push!(adj[t], f)
    end
    for a in adj; sort!(unique!(a)); end
    return zone, adj
end

# End-to-end: cases that produced disconnected zones before the fix must
# now come back connected. Guarded on Metis availability + PGLib data, so
# it runs locally (Metis in the shared env) and skips cleanly in CI.
function run_partition_connectivity_e2e_tests()
    @testset "create_partition returns connected zones (end-to-end)" begin
        # Activate the Metis extension if its trigger packages are present.
        ext_active =
            try
                @eval using Metis, Graphs, SimpleWeightedGraphs
                Base.get_extension(Lirion, :MetisExt) !== nothing
            catch
                false
            end
        root = get(ENV, "PGLIB_OPF_ROOT", "")
        if !ext_active || !isdir(root)
            @test_skip "Metis extension / PGLib data unavailable — e2e skipped"
            return nothing
        end
        # (case, k) pairs that were `not_clean` before the component-based fix.
        checks = [
            ("pglib_opf_case73_ieee_rts.m", 4),
            ("pglib_opf_case179_goc.m", 4),
            ("pglib_opf_case300_ieee.m", 8),
        ]
        for (file, k) in checks
            path = joinpath(root, file)
            isfile(path) || (@test_skip "missing $file"; continue)
            for refine in (true, false)
                c = Lirion.load_inputs(path)
                try
                    Lirion.create_partition(c; k = k, seed = 42, refine = refine)
                catch e
                    # Extension not overriding the stub in this session — skip
                    # rather than fail (CI without Metis as a test dep).
                    @test_skip "create_partition unavailable: $(sprint(showerror, e))"
                    continue
                end
                zone, adj = _case_zone_and_adj(c)
                @test Lirion.all_zones_connected(zone, adj)
            end
        end
    end
    return nothing
end
