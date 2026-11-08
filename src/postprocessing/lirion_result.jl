# ============================================================
# LirionResult — user-facing result container
# ============================================================

using DataFrames

"""
    LirionResult

Result of a Lirion.solve call.

Fields (user-facing):
- `algorithm`       : `:centralized`, `:admm`, or `:aladin`
- `base_model`      : `:DC`, `:AC`, `:ACR`, `:SOCWR`, ...
- `coupling`        : `:none` (centralized) or `:state` / `:power` / `:hybrid`
- `sectioning`      : `:none` / `:tie_line` / `:dummy_bus` / `:split_bus`
- `objective`       : total generation cost
- `status`          : termination status (`:optimal`, `:local_optimal`,
                       `:converged`, `:iteration_limit`, `:stalled` — ALADIN
                       hit max_iter with a frozen iterate — etc.)
- `solve_time`      : wall-clock elapsed time in seconds
- `iter_done`       : iteration count (`nothing` for centralized)

Tables (DataFrames, mirror the input matpower schema):
- `bus`     : bus_i, zone, pd, qd, …, va, vm (if applicable), vm_sq (SOCWR)
- `gen`     : gen_bus, status, pg, qg (if applicable), curtailment
- `branch`  : from, to, status, pf, qf, pt, qt (flows at both ends)

Distributed extras (`nothing` for centralized):
- `zones`     : `Vector{NamedTuple}` with one entry per zone, each holding
                `(id, bus, gen, branch, cost)`
- `consensus` : DataFrame of boundary variables in long format — one row
                per (variable, element, zone), with columns `(variable,
                bus_id, branch_id, zone, value, consensus, deviation)`

Internals (advanced use):
- `model`       : JuMP model (centralized) or `nothing`
- `system`      : `SystemData` (centralized) or `DistributedSystem`
- `partitions`  : `Partitions` (distributed) or `nothing`
- `diagnostics` : iteration-by-iteration history (distributed)

Dump to CSV with `Lirion.to_csvs(result, "out_dir/")`.
"""
struct LirionResult
    algorithm::Symbol
    base_model::Symbol
    coupling::Symbol
    sectioning::Symbol
    objective::Float64
    status::Symbol
    solve_time::Float64
    iter_done::Union{Nothing, Int}

    bus::DataFrame
    gen::DataFrame
    branch::DataFrame

    zones::Union{Nothing, Vector{NamedTuple}}
    consensus::Union{Nothing, DataFrame}

    model::Any
    system::Any
    partitions::Any
    diagnostics::AbstractRunDiagnostics
end

# ------------------------------------------------------------
# Helpers to build DataFrames from a solved JuMP model + SystemData
# ------------------------------------------------------------

_safe_value(v) =
    try
        JuMP.value(v)
    catch
        NaN
    end

"""Build a DataFrame with one row per bus holding the optimal state.

`bus` must expose `.id`, `.zone`, `.pd`, `.qd` (SystemData.bus or PartitionDC.bus style)."""
function _build_bus_df(model::JuMP.Model, bus)
    ids = bus.id
    df = DataFrame(;
        bus_i = ids,
        zone = copy(bus.zone),
        pd = copy(bus.pd),
        qd = copy(bus.qd),
    )

    if haskey(model, :theta)
        df.va = [_safe_value(model[:theta][b, 1]) for b in ids]
    elseif haskey(model, :v_real) && haskey(model, :v_imag)
        df.va = [
            atan(_safe_value(model[:v_imag][b, 1]), _safe_value(model[:v_real][b, 1]))
            for b in ids
        ]
    end

    if haskey(model, :v_magnitude)
        df.vm = [_safe_value(model[:v_magnitude][b, 1]) for b in ids]
    elseif haskey(model, :v_real) && haskey(model, :v_imag)
        df.vm = [
            hypot(_safe_value(model[:v_real][b, 1]), _safe_value(model[:v_imag][b, 1]))
            for b in ids
        ]
    elseif haskey(model, :phi)
        df.vm = [1.0 + _safe_value(model[:phi][b, 1]) for b in ids]
    end

    if haskey(model, :v_magnitude_sq)
        df.vm_sq = [_safe_value(model[:v_magnitude_sq][b, 1]) for b in ids]
    end

    return df
end

"""Build a DataFrame with one row per (in-service) generator.

`gen` must expose `.n`, `.status`, `.gen_bus` (SystemData.gen-style). Returns
an empty DataFrame when the zone has no generators or the model doesn't expose
`p_generation` (purely passive zone)."""
function _build_gen_df(model::JuMP.Model, gen)
    gen === nothing && return DataFrame()
    gen.n == 0 && return DataFrame()
    haskey(model, :p_generation) || return DataFrame()
    pg_var = model[:p_generation]
    # Filter by generators actually present in the JuMP model, not just by the
    # dataset's status flag. Preprocessing can introduce mismatches.
    ids = Int[]
    for g in 1:gen.n
        gen.status[g] == 1 || continue
        _pg_has_key(pg_var, g) || continue
        push!(ids, g)
    end
    isempty(ids) && return DataFrame()
    df = DataFrame(;
        gen_id = ids,
        gen_bus = gen.bus[ids],
        status = gen.status[ids],
        pg = [_safe_value(pg_var[g, 1]) for g in ids],
    )
    if haskey(model, :q_generation)
        qg_var = model[:q_generation]
        df.qg = [_pg_has_key(qg_var, g) ? _safe_value(qg_var[g, 1]) : NaN for g in ids]
    end
    if haskey(model, :curtailment)
        cv = model[:curtailment]
        df.curtailment = [_pg_has_key(cv, g) ? _safe_value(cv[g, 1]) : 0.0 for g in ids]
    end
    return df
end

# Support both DenseAxisArray (hasprimary key lookup via `haskey`) and
# SparseAxisArray (supports `haskey` on the key tuple).
@inline function _pg_has_key(container, g::Int)
    try
        if container isa JuMP.Containers.DenseAxisArray
            return g in axes(container, 1)
        elseif container isa JuMP.Containers.SparseAxisArray
            return haskey(container, (g, 1))
        end
        return true
    catch
        return false
    end
end

"""Build a DataFrame with one row per (in-service) branch.

When `zone_branch_idx_in_full` is supplied (distributed partitions), the
returned DataFrame also carries a `branch_id_global` column — the index of
each row in the ORIGINAL pre-partition case. This is what the top-level
dedupe uses to merge duplicates (shared tie-lines) across zones."""
function _build_branch_df(model::JuMP.Model, br; zone_branch_idx_in_full = nothing)
    br.n == 0 && return DataFrame()
    haskey(model, :p_flow_branches) || return DataFrame()
    ids = [l for l in 1:br.n if br.status[l] == 1]
    df = DataFrame(;
        branch_id = ids,
        from = br.from[ids],
        to = br.to[ids],
        status = br.status[ids],
        pf = [_safe_value(model[:p_flow_branches][l, 1]) for l in ids],
        pt = [_safe_value(model[:p_reverse_flow_branches][l, 1]) for l in ids],
    )
    if haskey(model, :q_flow_branches)
        df.qf = [_safe_value(model[:q_flow_branches][l, 1]) for l in ids]
    end
    if haskey(model, :q_reverse_flow_branches)
        df.qt = [_safe_value(model[:q_reverse_flow_branches][l, 1]) for l in ids]
    end
    if zone_branch_idx_in_full !== nothing
        df.branch_id_global = [zone_branch_idx_in_full[l] for l in ids]
    end
    return df
end

# ------------------------------------------------------------
# Status inference
# ------------------------------------------------------------

function _infer_status(model::JuMP.Model)::Symbol
    st = try
        JuMP.termination_status(model)
    catch
        return :unknown
    end
    st == MOI.OPTIMAL && return :optimal
    st == MOI.LOCALLY_SOLVED && return :local_optimal
    st == MOI.ALMOST_OPTIMAL && return :almost_optimal
    st == MOI.ITERATION_LIMIT && return :iteration_limit
    st == MOI.TIME_LIMIT && return :time_limit
    st == MOI.INFEASIBLE && return :infeasible
    return Symbol(lowercase(string(st)))
end

# ------------------------------------------------------------
# Centralized result builder (drop-in replacement for build_result_dict)
# ------------------------------------------------------------

function build_lirion_result_centralized(
    model::JuMP.Model,
    sys::SystemData,
    model_type::Type{<:OPFModel},
    t0_ns::Integer,
)::LirionResult
    bus_df = _build_bus_df(model, sys.bus)
    gen_df = _build_gen_df(model, sys.gen)
    branch_df = _build_branch_df(model, sys.branch)
    obj = JuMP.objective_value(model)
    elapsed_s = (time_ns() - t0_ns) / 1e9

    return LirionResult(
        :centralized,
        Symbol(nameof(model_type)),
        :none,
        :none,
        Float64(obj),
        _infer_status(model),
        Float64(elapsed_s),
        nothing,
        bus_df,
        gen_df,
        branch_df,
        nothing,
        nothing,
        model,
        sys,
        nothing,
        NoDiagnostics(),
    )
end

# ------------------------------------------------------------
# Distributed zone slicing — per-partition DataFrames
# ------------------------------------------------------------

function _build_zone_result(part::AbstractPartition, zone_sys::SystemData)::NamedTuple
    model = part.local_model::JuMP.Model
    bus_df = _build_bus_df(model, part.bus)
    gen_df = _build_gen_df(model, zone_sys.gen)
    # Pass the global mapping so the zone DF can expose `branch_id_global`.
    branch_df = _build_branch_df(
        model, part.branch;
        zone_branch_idx_in_full = part.zone_branch_idx_in_full,
    )
    return (
        id = part.id,
        bus = bus_df,
        gen = gen_df,
        branch = branch_df,
        cost = part.generation_cost,
    )
end

# ------------------------------------------------------------
# Consensus DataFrame (distributed)
# ------------------------------------------------------------

"""Build a DataFrame summarising boundary consensus residuals."""
function _build_consensus_df(ws::Any, partitions, sys)::DataFrame
    rows = NamedTuple[]

    is_aladin =
        hasproperty(ws, :zones) &&
        hasproperty(ws, :bus_pairs) && !isempty(ws.bus_pairs) &&
        !(eltype(ws.bus_pairs) <: NTuple)

    if is_aladin
        _push_bus_consensus_rows!(rows, ws, partitions)
    elseif sys !== nothing && !isempty(sys.full_shared_bus_ids)
        _push_admm_bus_consensus_rows!(rows, partitions, sys)
    end

    has_edge_pairs =
        hasproperty(ws, :interzonal_line_pairs) && !isempty(ws.interzonal_line_pairs) &&
        !(eltype(ws.interzonal_line_pairs) <: NTuple)
    if has_edge_pairs
        _push_edge_consensus_rows!(rows, ws, partitions)
    elseif hasproperty(ws, :interzonal) && !isempty(ws.interzonal)
        _push_admm_interzonal_rows!(rows, ws, partitions)
    end

    # DummyBus half-branch consensus. ALADIN uses `dummy_line_pairs`
    # (Vector of struct); ADMM uses `branch_pairs::Vector{NTuple{4,Int}}`.
    if hasproperty(ws, :dummy_line_pairs) && !isempty(ws.dummy_line_pairs)
        _push_dummy_consensus_rows!(rows, ws, partitions)
    elseif hasproperty(ws, :branch_pairs) && !isempty(ws.branch_pairs) &&
           eltype(ws.branch_pairs) <: NTuple
        _push_admm_dummy_rows!(rows, ws, partitions)
    end

    # SplitBus p_inj / q_inj consensus (ALADIN uses `split_bus_pairs`;
    # ADMM uses `bus_pairs::Vector{NTuple{3,Int}}` plus `p_inj_vars`).
    if hasproperty(ws, :split_bus_pairs) && !isempty(ws.split_bus_pairs)
        _push_split_bus_consensus_rows!(rows, ws, partitions)
    elseif hasproperty(ws, :bus_pairs) && !isempty(ws.bus_pairs) &&
           eltype(ws.bus_pairs) <: NTuple
        _push_admm_split_bus_consensus_rows!(rows, ws, partitions)
    end

    isempty(rows) && return DataFrame(;
        variable = Symbol[],
        bus_id = Union{Int, Missing}[],
        branch_id = Union{Int, Missing}[],
        zone = Int[],
        value = Float64[],
        consensus = Float64[],
        deviation = Float64[],
    )
    # Aggregate in long format: one row per (variable, element, zone), with
    # `consensus` = mean across zones sharing the element, and
    # `deviation = value - consensus`.
    df = DataFrame(rows)
    gdf = DataFrames.groupby(df, [:variable, :bus_id, :branch_id]; sort = false)
    out = DataFrames.transform(
        gdf,
        :value => (v -> fill(sum(v) / length(v), length(v))) => :consensus,
    )
    out.deviation = out.value .- out.consensus
    return out
end

# Emit a single reading: "zone k saw `value` for variable `var` at this
# physical element (bus_id OR branch_id, whichever applies)".
@inline function _push_bus_reading!(
    rows, variable::Symbol, bus_id::Int, zone::Int, value::Float64,
)
    push!(
        rows,
        (
            variable = variable,
            bus_id = bus_id,
            branch_id = missing,
            zone = zone,
            value = value,
        ),
    )
    return nothing
end

@inline function _push_branch_reading!(
    rows, variable::Symbol, branch_id::Int, zone::Int, value::Float64,
)
    push!(
        rows,
        (
            variable = variable,
            bus_id = missing,
            branch_id = branch_id,
            zone = zone,
            value = value,
        ),
    )
    return nothing
end

# --- ALADIN bus pairs (SharedBusPair / SharedBusPairAC / SharedBusPairSOCWR) ---
function _push_bus_consensus_rows!(rows, ws, partitions)
    for pair in ws.bus_pairs
        zi = pair.zone_i
        zj = pair.zone_j
        b_id = pair.bus_id
        if hasproperty(pair, :col_i)
            _push_bus_reading!(rows, :theta, b_id, zi, ws.zones[zi].y[pair.col_i])
            _push_bus_reading!(rows, :theta, b_id, zj, ws.zones[zj].y[pair.col_j])
        elseif hasproperty(pair, :theta_col_i)
            _push_bus_reading!(rows, :theta, b_id, zi, ws.zones[zi].y[pair.theta_col_i])
            _push_bus_reading!(rows, :theta, b_id, zj, ws.zones[zj].y[pair.theta_col_j])
            _push_bus_reading!(rows, :v_mag, b_id, zi, ws.zones[zi].y[pair.vm_col_i])
            _push_bus_reading!(rows, :v_mag, b_id, zj, ws.zones[zj].y[pair.vm_col_j])
        elseif hasproperty(pair, :v2_col_i)
            _push_bus_reading!(rows, :v2, b_id, zi, ws.zones[zi].y[pair.v2_col_i])
            _push_bus_reading!(rows, :v2, b_id, zj, ws.zones[zj].y[pair.v2_col_j])
        end
    end
    return nothing
end

function _push_edge_consensus_rows!(rows, ws, partitions)
    for pair in ws.interzonal_line_pairs
        zi = pair.zone_i
        zj = pair.zone_j
        mi = partitions.parts[zi].local_model
        mj = partitions.parts[zj].local_model
        le_i = pair.line_local_i
        le_j = pair.line_local_j
        gi = partitions.parts[zi].zone_edges[le_i]
        li = partitions.parts[zi].branch_idx_full_to_zone[gi]
        gj = partitions.parts[zj].zone_edges[le_j]
        lj = partitions.parts[zj].branch_idx_full_to_zone[gj]
        b_id = pair.line_global
        if haskey(mi, :p_flow_branches) && haskey(mj, :p_flow_branches)
            _push_branch_reading!(rows, :p_flow, b_id, zi,
                _safe_value(mi[:p_flow_branches][li, 1]))
            _push_branch_reading!(rows, :p_flow, b_id, zj,
                _safe_value(mj[:p_flow_branches][lj, 1]))
            if haskey(mi, :q_flow_branches) && haskey(mj, :q_flow_branches)
                _push_branch_reading!(rows, :q_flow, b_id, zi,
                    _safe_value(mi[:q_flow_branches][li, 1]))
                _push_branch_reading!(rows, :q_flow, b_id, zj,
                    _safe_value(mj[:q_flow_branches][lj, 1]))
            end
        end
    end
    return nothing
end

function _push_dummy_consensus_rows!(rows, ws, partitions)
    for pair in ws.dummy_line_pairs
        zi = pair.zone_i
        zj = pair.zone_j
        mi = partitions.parts[zi].local_model
        mj = partitions.parts[zj].local_model
        li = hasproperty(pair, :local_i) ? pair.local_i : nothing
        li === nothing && continue
        lj = pair.local_j
        dummy_bus_id = partitions.parts[zj].branch.from[lj]
        if haskey(mi, :p_flow_branches) && haskey(mj, :p_flow_branches)
            _push_bus_reading!(rows, :p_flow_dummy, dummy_bus_id, zi,
                _safe_value(mi[:p_flow_branches][li, 1]))
            _push_bus_reading!(rows, :p_flow_dummy, dummy_bus_id, zj,
                _safe_value(mj[:p_flow_branches][lj, 1]))
            if haskey(mi, :q_flow_branches) && haskey(mj, :q_flow_branches)
                _push_bus_reading!(rows, :q_flow_dummy, dummy_bus_id, zi,
                    _safe_value(mi[:q_flow_branches][li, 1]))
                _push_bus_reading!(rows, :q_flow_dummy, dummy_bus_id, zj,
                    _safe_value(mj[:q_flow_branches][lj, 1]))
            end
        end
    end
    return nothing
end

# ADMM: iterate sys.full_shared_bus_ids and collect per-zone θ/v/v² values
# from each partition's solved JuMP model.
function _push_admm_bus_consensus_rows!(rows, partitions, sys)
    shared = sys.full_shared_bus_ids
    isempty(shared) && return nothing

    parts = partitions.parts
    for b_id in shared
        hosts = Tuple{Int, JuMP.Model}[]
        for (zpos, part) in enumerate(parts)
            b_id in part.bus.id && push!(hosts, (zpos, part.local_model::JuMP.Model))
        end
        length(hosts) >= 2 || continue
        # Emit one reading per zone that hosts this shared bus.
        for (zone, m) in hosts
            haskey(m, :theta) && _push_bus_reading!(
                rows, :theta, b_id, zone, _safe_value(m[:theta][b_id, 1]))
            haskey(m, :v_magnitude) && _push_bus_reading!(
                rows, :v_mag, b_id, zone, _safe_value(m[:v_magnitude][b_id, 1]))
            haskey(m, :v_magnitude_sq) && _push_bus_reading!(
                rows, :v2, b_id, zone, _safe_value(m[:v_magnitude_sq][b_id, 1]))
        end
    end
    return nothing
end

function _push_admm_dummy_rows!(rows, ws, partitions)
    parts = partitions.parts
    for (ga, za, gb, zb) in ws.branch_pairs
        pa = parts[za]
        pb = parts[zb]
        la = pa.branch_idx_full_to_zone[ga]
        lb = pb.branch_idx_full_to_zone[gb]
        la == 0 && continue
        lb == 0 && continue
        mi = pa.local_model::JuMP.Model
        mj = pb.local_model::JuMP.Model
        dummy_bus_id = pb.branch.from[lb]
        if haskey(mi, :p_flow_branches) && haskey(mj, :p_flow_branches)
            _push_bus_reading!(rows, :p_flow_dummy, dummy_bus_id, za,
                _safe_value(mi[:p_flow_branches][la, 1]))
            _push_bus_reading!(rows, :p_flow_dummy, dummy_bus_id, zb,
                _safe_value(mj[:p_flow_branches][lb, 1]))
        end
        if haskey(mi, :q_flow_branches) && haskey(mj, :q_flow_branches)
            _push_bus_reading!(rows, :q_flow_dummy, dummy_bus_id, za,
                _safe_value(mi[:q_flow_branches][la, 1]))
            _push_bus_reading!(rows, :q_flow_dummy, dummy_bus_id, zb,
                _safe_value(mj[:q_flow_branches][lb, 1]))
        end
    end
    return nothing
end

function _push_admm_interzonal_rows!(rows, ws, partitions)
    parts = partitions.parts
    for e in ws.interzonal
        hosts = Tuple{Int, Int}[]
        for (zpos, part) in enumerate(parts)
            le = part.branch_idx_full_to_zone[e]
            le != 0 && push!(hosts, (zpos, le))
        end
        length(hosts) >= 2 || continue
        for (zone, local_idx) in hosts
            m = parts[zone].local_model::JuMP.Model
            haskey(m, :p_flow_branches) && _push_branch_reading!(
                rows, :p_flow, e, zone,
                _safe_value(m[:p_flow_branches][local_idx, 1]))
            haskey(m, :q_flow_branches) && _push_branch_reading!(
                rows, :q_flow, e, zone,
                _safe_value(m[:q_flow_branches][local_idx, 1]))
        end
    end
    return nothing
end

function _push_split_bus_consensus_rows!(rows, ws, partitions)
    hasproperty(ws, :zones) || return nothing
    for pair in ws.split_bus_pairs
        b_id = pair.bus_id
        has_q = hasproperty(pair, :q_inj_cols)
        for c in eachindex(pair.zones)
            zone = pair.zones[c]
            _push_bus_reading!(
                rows,
                :p_inj,
                b_id,
                zone,
                ws.zones[zone].y[pair.p_inj_cols[c]],
            )
            has_q && _push_bus_reading!(
                rows, :q_inj, b_id, zone, ws.zones[zone].y[pair.q_inj_cols[c]],
            )
        end
    end
    return nothing
end

function _push_admm_split_bus_consensus_rows!(rows, ws, partitions)
    hasproperty(ws, :p_inj_vars) || return nothing
    p_inj_vars = ws.p_inj_vars
    q_inj_vars = hasproperty(ws, :q_inj_vars) ? ws.q_inj_vars : nothing
    # One row per (split bus, zone copy); iterating the per-zone stores keeps a
    # bus shared by >2 zones from being counted once per pairwise edge.
    for zone in eachindex(p_inj_vars)
        for b_id in keys(p_inj_vars[zone])
            v = _admm_pinj_value(partitions.parts[zone], p_inj_vars[zone], b_id)
            v === nothing && continue
            _push_bus_reading!(rows, :p_inj, b_id, zone, v)
            if q_inj_vars !== nothing
                q = _admm_pinj_value(partitions.parts[zone], q_inj_vars[zone], b_id)
                q !== nothing && _push_bus_reading!(rows, :q_inj, b_id, zone, q)
            end
        end
    end
    return nothing
end

@inline function _admm_pinj_value(part, store, b_id)
    v = get(store, b_id, nothing)
    v === nothing && return nothing
    m = part.local_model::JuMP.Model
    JuMP.has_values(m) || return nothing
    return _safe_value(v)
end

# ------------------------------------------------------------
# Top-level table reconciliation (distributed)
# ------------------------------------------------------------
#
# Rules:
#   bus    — group by `bus_i`. `pd/qd` keep the owner's values (ghosts carry 0);
#            `va/vm/vm_sq` averaged across zones seeing this bus (= consensus
#             target at convergence).
#   gen    — each generator lives in exactly one zone (no overlap). Passthrough.
#   branch — see _merge_branch_rows below; depends on the sectioning.
# ------------------------------------------------------------

function _merge_bus_rows(zones_vec)
    isempty(zones_vec) && return DataFrame()
    df = vcat([z.bus for z in zones_vec]...; cols = :union)
    ids = sort!(unique(df.bus_i))
    rows = NamedTuple[]
    for b_id in ids
        sub = df[df.bus_i.==b_id, :]
        # Owner's row is the one with non-zero demand (ghosts have pd=qd=0
        # by construction in `create_partition`). Fall back to the first row.
        own_idx = findfirst(i -> sub.pd[i] != 0 || sub.qd[i] != 0, 1:nrow(sub))
        own = own_idx === nothing ? 1 : own_idx
        row = (; bus_i = b_id, zone = sub.zone[own], pd = sub.pd[own], qd = sub.qd[own])
        # Mean the consensus-carrying columns across zones that see the bus.
        for col in (:va, :vm, :vm_sq)
            if col in propertynames(sub)
                vals = collect(skipmissing(getproperty(sub, col)))
                isempty(vals) && continue
                row = merge(row, (; col => sum(vals) / length(vals)))
            end
        end
        push!(rows, row)
    end
    return DataFrame(rows)
end

function _merge_gen_rows(zones_vec)
    isempty(zones_vec) && return DataFrame()
    return vcat([z.gen for z in zones_vec]...; cols = :union)
end

# Merge branch rows across zones. Returns one row per original-network branch,
# with shared rows averaged. For DummyBus, halves are fused back to the
# original tie-line using `sys.branch_pairs` + the zone-local branch data.
function _merge_branch_rows(zones_vec, sys, config)
    isempty(zones_vec) && return DataFrame()
    raw = vcat([z.branch for z in zones_vec]...; cols = :union)
    if !(config.sectioning isa DummyBus)
        return _dedupe_branch_rows_by_global(raw)
    end
    # DummyBus: each tie-line became two halves. Pair them up via
    # `sys.branch_pairs = [(ga_expanded, za, gb_expanded, zb), …]`.
    return _fuse_dummybus_branch_rows(raw, zones_vec, sys)
end

function _dedupe_branch_rows_by_global(raw::DataFrame)
    "branch_id_global" in string.(propertynames(raw)) || return raw
    ids = sort!(unique(raw.branch_id_global))
    rows = NamedTuple[]
    for gid in ids
        sub = raw[raw.branch_id_global.==gid, :]
        row = (;
            branch_id_global = gid,
            from = sub.from[1],
            to = sub.to[1],
            status = sub.status[1],
        )
        # Average the flow columns across zones.
        for col in (:pf, :pt, :qf, :qt)
            if col in propertynames(sub)
                vals = collect(skipmissing(getproperty(sub, col)))
                isempty(vals) && continue
                row = merge(row, (; col => sum(vals) / length(vals)))
            end
        end
        push!(rows, row)
    end
    return DataFrame(rows)
end

function _fuse_dummybus_branch_rows(raw::DataFrame, zones_vec, sys)
    # Expanded-case global indices of the two halves per tie-line.
    pairs = sys.branch_pairs
    isempty(pairs) && return _dedupe_branch_rows_by_global(raw)
    # Collect "half" rows indexed by their branch_id_global (expanded-case idx).
    half_rows = Dict{Int, DataFrameRow}()
    for row in eachrow(raw)
        ismissing(row.branch_id_global) && continue
        half_rows[row.branch_id_global] = row
    end
    # Fuse pairs back to the original tie-line row. Keep any non-half rows
    # (internal branches) untouched.
    halves_set = Set{Int}()
    fused_rows = NamedTuple[]
    for (ga, za, gb, zb) in pairs
        push!(halves_set, ga)
        push!(halves_set, gb)
        ra = get(half_rows, ga, nothing)
        rb = get(half_rows, gb, nothing)
        (ra === nothing || rb === nothing) && continue
        from = ra.from       # original tie-line's `from`
        to = rb.to           # original tie-line's `to`
        nt = (;
            # Synthetic: we use half_a's expanded idx as the key so users can
            # still cross-reference. Could also add a `sectioning_notes` column.
            branch_id_global = ga,
            from = from,
            to = to,
            status = ra.status,
        )
        for col in (:pf, :pt, :qf, :qt)
            if col in propertynames(ra) && col in propertynames(rb)
                va = getproperty(ra, col)
                vb = getproperty(rb, col)
                if !ismissing(va) && !ismissing(vb)
                    # Halves carry the flow from each side — average gives the
                    # physical flow through the pre-partition tie-line.
                    nt = merge(nt, (; col => (va + vb) / 2))
                end
            end
        end
        push!(fused_rows, nt)
    end
    # Non-pair branches (internal) passthrough.
    for row in eachrow(raw)
        ismissing(row.branch_id_global) && continue
        row.branch_id_global in halves_set && continue
        nt = (; branch_id_global = row.branch_id_global,
            from = row.from, to = row.to, status = row.status)
        for col in (:pf, :pt, :qf, :qt)
            if col in propertynames(row)
                v = getproperty(row, col)
                ismissing(v) || (nt = merge(nt, (; col => v)))
            end
        end
        push!(fused_rows, nt)
    end
    return DataFrame(fused_rows)
end

# ------------------------------------------------------------
# Distributed result builder
# ------------------------------------------------------------

function build_lirion_result_distributed(
    partitions,
    ws,
    config,
    iter_done::Integer,
    t0_ns::Integer,
    sys,
)::LirionResult
    zones_vec = NamedTuple[]
    sizehint!(zones_vec, length(partitions.parts))
    for (i, part) in enumerate(partitions.parts)
        push!(zones_vec, _build_zone_result(part, sys.parts[i]))
    end

    # Top-level: reconcile per-zone views into a single physical-network view,
    # directly comparable to a centralized result. See `_merge_*` helpers.
    bus_df = _merge_bus_rows(zones_vec)
    gen_df = _merge_gen_rows(zones_vec)
    branch_df = _merge_branch_rows(zones_vec, sys, config)
    total_cost = sum(z.cost for z in zones_vec; init = 0.0)
    elapsed_s = (time_ns() - t0_ns) / 1e9

    cons_df = _build_consensus_df(ws, partitions, sys)
    diag = hasproperty(ws, :hist) ? ws.hist : NoDiagnostics()

    algo_sym =
        config isa ADMMConfig ? :admm :
        config isa ALADINConfig ? :aladin :
        :distributed
    max_iter = hasproperty(config, :max_iter) ? config.max_iter : 0
    # recorded terminal status wins; iter_done heuristic is the fallback
    ts = hasproperty(diag, :meta) ? get(diag.meta, :terminal_status, nothing) : nothing
    status_sym = ts isa Symbol ? ts :
        (iter_done >= max_iter ? :iteration_limit : :converged)

    return LirionResult(
        algo_sym,
        Symbol(nameof(config.base_model)),
        config.coupling,
        Symbol(lowercase(string(nameof(typeof(config.sectioning))))),
        Float64(total_cost),
        status_sym,
        Float64(elapsed_s),
        Int(iter_done),
        bus_df,
        gen_df,
        branch_df,
        zones_vec,
        cons_df,
        nothing,
        sys,
        partitions,
        diag,
    )
end

# ------------------------------------------------------------
# Pretty-printing
# ------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", r::LirionResult)
    println(io, "LirionResult")
    if r.algorithm === :centralized
        println(io, "  ├─ algorithm:  centralized")
        println(io, "  ├─ model:      ", r.base_model)
    else
        println(io, "  ├─ algorithm:  ", uppercase(string(r.algorithm)))
        println(io, "  ├─ model:      ", r.base_model, " / ", r.coupling,
            " / ", r.sectioning)
    end
    println(io, "  ├─ objective:  ", round(r.objective; digits = 4))
    println(io, "  ├─ status:     ", r.status)
    println(io, "  ├─ solve_time: ", round(r.solve_time; digits = 3), " s")
    r.iter_done === nothing || println(io, "  ├─ iter_done:  ", r.iter_done)
    println(io, "  ├─ bus:        ", size(r.bus, 1), " rows × ",
        size(r.bus, 2), " cols")
    println(io, "  ├─ gen:        ", size(r.gen, 1), " rows × ",
        size(r.gen, 2), " cols")
    print(io, "  ", r.zones === nothing ? "└─" : "├─", " branch:     ",
        size(r.branch, 1), " rows × ", size(r.branch, 2), " cols")
    if r.zones !== nothing
        println(io)
        println(io, "  ├─ zones:      ", length(r.zones), " zone(s)")
        if r.consensus !== nothing
            print(io, "  └─ consensus:  ", size(r.consensus, 1), " row(s)")
        else
            print(io, "  └─ consensus:  —")
        end
    end
    return nothing
end

Base.show(io::IO, ::LirionResult) = print(io, "LirionResult(...)")

# Alias old name so existing user code keeps working.
const LirionRunResult = LirionResult

# ------------------------------------------------------------
# CSV export
# ------------------------------------------------------------

"""
    to_csvs(result::LirionResult, outdir::AbstractString)

Dump the result's DataFrames as CSV files under `outdir` (created if missing).

Files written:
- `summary.csv`   : one row with algorithm / model / objective / status / time
- `bus.csv`, `gen.csv`, `branch.csv`   : merged across zones for distributed
- `consensus.csv` (distributed only)   : boundary-variable residuals
- `zone_<id>/{bus,gen,branch}.csv`     : per-zone tables (distributed)

Requires `CSV` to be loaded (`using CSV`). This function is defined via the
`CSVExt` extension.
"""
function to_csvs end
