# ============================================================
# ALADIN-AC — Zone-level data layout & builder
# ============================================================

# ============================================================
# Struct
# ============================================================

mutable struct AladinZoneDataAC <: AbstractAladinZoneData
    local_model::JuMP.Model

    all_vars::Vector{JuMP.VariableRef}
    n_vars::Int
    var_to_col::Dict{JuMP.VariableRef, Int}

    theta_cols::Vector{Int}
    vm_cols::Vector{Int}
    pg_cols::Vector{Int}
    qg_cols::Vector{Int}
    curt_cols::Vector{Int}

    pg_c2_sbase2::Vector{Float64}     # 2·c2·sbase²
    pg_c1_sbase::Vector{Float64}      # c1·sbase

    shared_theta_cols::Vector{Int}
    shared_vm_cols::Vector{Int}

    edge_from_theta_col::Vector{Int}
    edge_to_theta_col::Vector{Int}
    edge_from_vm_col::Vector{Int}
    edge_to_vm_col::Vector{Int}
    edge_g::Vector{Float64}
    edge_b::Vector{Float64}
    edge_b_sh::Vector{Float64}
    edge_inv_tap::Vector{Float64}
    edge_tap_angle::Vector{Float64}

    lb::Vector{Float64}
    ub::Vector{Float64}

    x_prev::Vector{Float64}
    y::Vector{Float64}

    g::Vector{Float64}
    hess::SparseArrays.SparseMatrixCSC{Float64, Int}

    active_jac::SparseArrays.SparseMatrixCSC{Float64, Int}

    delta::Vector{Float64}

    active_bound_cols::Vector{Int}
    quad_static::JuMP.QuadExpr
end

# ============================================================
# Variable / cost / topology collectors
# ============================================================

# ------------------------------------------------------------
# Column classes (θ / vm / pg / qg / curt) + cost coefficients
# ------------------------------------------------------------

function _collect_ac_cols(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    obj::ObjectiveData,
)
    theta_cols = Int[]
    vm_cols = Int[]
    pg_cols = Int[]
    qg_cols = Int[]
    curt_cols = Int[]

    pg_c2_sbase2 = Float64[]
    pg_c1_sbase = Float64[]

    sbase = obj.sbase

    if haskey(model, :theta)
        for v in values(model[:theta].data)
            col = get(var_to_col, v, 0)
            col == 0 && continue
            push!(theta_cols, col)
        end
    end

    if haskey(model, :v_magnitude)
        for v in values(model[:v_magnitude].data)
            col = get(var_to_col, v, 0)
            col == 0 && continue
            push!(vm_cols, col)
        end
    end

    if haskey(model, :p_generation)
        pg_var = model[:p_generation]
        for g in 1:obj.number_of_generators
            obj.generator_status[g] == 1 || continue
            v = pg_var[g, 1]
            col = get(var_to_col, v, 0)
            col == 0 && continue
            push!(pg_cols, col)
            coefs = obj.generator_costs[1][g]
            @assert length(coefs) <= 3 "ALADIN: generator cost must be at most quadratic (got degree $(length(coefs)-1))"
            c2 = (length(coefs) >= 3) ? coefs[end-2] : 0.0
            c1 = (length(coefs) >= 2) ? coefs[end-1] : 0.0
            push!(pg_c2_sbase2, 2.0 * c2 * sbase^2)
            push!(pg_c1_sbase, c1 * sbase)
        end
    end

    if haskey(model, :q_generation)
        qg_var = model[:q_generation]
        for g in 1:obj.number_of_generators
            obj.generator_status[g] == 1 || continue
            v = qg_var[g, 1]
            col = get(var_to_col, v, 0)
            col == 0 && continue
            push!(qg_cols, col)
        end
    end

    if haskey(model, :curtailment)
        for v in values(model[:curtailment].data)
            col = get(var_to_col, v, 0)
            col == 0 && continue
            push!(curt_cols, col)
        end
    end

    return (
        theta_cols = theta_cols,
        vm_cols = vm_cols,
        pg_cols = pg_cols,
        qg_cols = qg_cols,
        curt_cols = curt_cols,
        pg_c2_sbase2 = pg_c2_sbase2,
        pg_c1_sbase = pg_c1_sbase,
    )
end

# ------------------------------------------------------------
# Shared θ / vm columns
# ------------------------------------------------------------

function _collect_shared_ac_cols(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    shared_bus_ids::Vector{Int},
)
    theta_var = model[:theta]
    vm_var = model[:v_magnitude]
    theta_cols = Vector{Int}(undef, length(shared_bus_ids))
    vm_cols = Vector{Int}(undef, length(shared_bus_ids))
    @inbounds for i in eachindex(shared_bus_ids)
        b = shared_bus_ids[i]
        theta_cols[i] = var_to_col[theta_var[b, 1]]
        vm_cols[i] = var_to_col[vm_var[b, 1]]
    end
    return theta_cols, vm_cols
end

# ------------------------------------------------------------
# Interzonal edge layout
# ------------------------------------------------------------

function _collect_ac_edge_layout(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    part::PartitionAC;
    include_edges::Bool,
)
    include_edges || return (
        edge_from_theta_col = Int[],
        edge_to_theta_col = Int[],
        edge_from_vm_col = Int[],
        edge_to_vm_col = Int[],
        edge_g = Float64[],
        edge_b = Float64[],
        edge_b_sh = Float64[],
        edge_inv_tap = Float64[],
        edge_tap_angle = Float64[],
    )

    zone_edges = part.zone_edges
    n = length(zone_edges)

    from_theta = Vector{Int}(undef, n)
    to_theta = Vector{Int}(undef, n)
    from_vm = Vector{Int}(undef, n)
    to_vm = Vector{Int}(undef, n)
    g_vec = Vector{Float64}(undef, n)
    b_vec = Vector{Float64}(undef, n)
    b_sh_vec = Vector{Float64}(undef, n)
    inv_tap_vec = Vector{Float64}(undef, n)
    tap_angle_vec = Vector{Float64}(undef, n)

    theta_var = model[:theta]
    vm_var = model[:v_magnitude]
    branch = part.branch
    branch_map = part.branch_idx_full_to_zone

    @inbounds for (k, e) in enumerate(zone_edges)
        l = branch_map[e]
        from_b = branch.from[l]
        to_b = branch.to[l]

        from_theta[k] = var_to_col[theta_var[from_b, 1]]
        to_theta[k] = var_to_col[theta_var[to_b, 1]]
        from_vm[k] = var_to_col[vm_var[from_b, 1]]
        to_vm[k] = var_to_col[vm_var[to_b, 1]]

        g_vec[k] = branch.g[l]
        b_vec[k] = branch.b[l]
        b_sh_vec[k] = branch.b_sh[l]
        inv_tap_vec[k] = branch.inv_tap[l]
        tap_angle_vec[k] = branch.tap_angle[l]
    end

    return (
        edge_from_theta_col = from_theta,
        edge_to_theta_col = to_theta,
        edge_from_vm_col = from_vm,
        edge_to_vm_col = to_vm,
        edge_g = g_vec,
        edge_b = b_vec,
        edge_b_sh = b_sh_vec,
        edge_inv_tap = inv_tap_vec,
        edge_tap_angle = tap_angle_vec,
    )
end

# ============================================================
# Builder
# ============================================================

function build_aladin_zone_data(
    part::PartitionAC;
    include_edges::Bool = false,
)::AladinZoneDataAC
    obj = part.objective_params
    model = part.local_model::JuMP.Model

    all_vars, n_vars, var_to_col = _build_var_layout(model)
    cols = _collect_ac_cols(model, var_to_col, obj)
    shared_theta_cols, shared_vm_cols =
        _collect_shared_ac_cols(model, var_to_col, part.zone_shared_bus_ids)
    edata = _collect_ac_edge_layout(model, var_to_col, part; include_edges)
    lb, ub = _collect_bounds(all_vars)

    hess_initial = SparseArrays.spzeros(Float64, n_vars, n_vars)
    jac_initial = SparseArrays.spzeros(Float64, 0, n_vars)

    return AladinZoneDataAC(
        model,
        all_vars,
        n_vars,
        var_to_col,
        cols.theta_cols,
        cols.vm_cols,
        cols.pg_cols,
        cols.qg_cols,
        cols.curt_cols,
        cols.pg_c2_sbase2,
        cols.pg_c1_sbase,
        shared_theta_cols,
        shared_vm_cols,
        edata.edge_from_theta_col,
        edata.edge_to_theta_col,
        edata.edge_from_vm_col,
        edata.edge_to_vm_col,
        edata.edge_g,
        edata.edge_b,
        edata.edge_b_sh,
        edata.edge_inv_tap,
        edata.edge_tap_angle,
        lb,
        ub,
        _init_from_start_values(all_vars),  # x_prev ← JuMP warm-start
        zeros(Float64, n_vars),             # y
        zeros(Float64, n_vars),             # g
        hess_initial,
        jac_initial,
        zeros(Float64, n_vars),             # delta
        Int[],                              # active_bound_cols
        zero(JuMP.QuadExpr),                # quad_static
    )
end
