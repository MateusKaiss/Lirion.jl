# ============================================================
# ADMM-DC
# ============================================================

# ------------------------------------------------------------
# Partition initialization
# ------------------------------------------------------------

# ---------- dc bus partition ----------

function init_dc_bus_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
)
    prepare_hotpaths!(part)
    model = build_model(sys_zone, DC; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.theta_var = model[:theta]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end

    build_quad_bus_expr!(part, rho_theta)
    @objective(model, Min, part.obj_base + part.quad_bus_expr)

    return nothing
end

function init_dc_bus_partition!(
    partitions::Partitions{PartitionDC},
    sys::ADMMSystem,
    optimizer_factory,
    rho_theta::Float64,
)
    n_zones = length(partitions.parts)
    @assert length(sys.parts) == n_zones "ADMMSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n_zones
        init_dc_bus_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_theta,
        )
    end
    return nothing
end

# ---------- dc branch partition ----------

function init_dc_branch_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho::Float64,
    interzonal_branches::Vector{Int},
)
    prepare_hotpaths!(part)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)

    model = build_model(sys_zone, DC; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.theta_var = model[:theta]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end

    build_quad_branch_expr!(part, rho, part.zone_edges)
    @objective(model, Min, part.obj_base + part.quad_branch_expr)

    return nothing
end

function init_dc_branch_partition!(
    partitions::Partitions{PartitionDC},
    sys::ADMMSystem,
    optimizer_factory,
    rho::Float64,
    interzonal_branches::Vector{Int},
)
    n_zones = length(partitions.parts)
    @assert length(sys.parts) == n_zones "ADMMSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n_zones
        init_dc_branch_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho,
            interzonal_branches,
        )
    end
    return nothing
end

# ---------- dc hybrid partition ----------

function init_dc_hybrid_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    interzonal_branches::Vector{Int},
)
    prepare_hotpaths!(part)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)

    model = build_model(sys_zone, DC; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.theta_var = model[:theta]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end

    build_quad_bus_expr!(part, rho_theta)
    build_quad_branch_expr!(part, rho, part.zone_edges)

    @objective(model, Min, part.obj_base + part.quad_bus_expr + part.quad_branch_expr)

    return nothing
end

function init_dc_hybrid_partition!(
    partitions::Partitions{PartitionDC},
    sys::ADMMSystem,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    interzonal_branches::Vector{Int},
)
    n_zones = length(partitions.parts)
    @assert length(sys.parts) == n_zones "ADMMSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n_zones
        init_dc_hybrid_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_theta,
            rho,
            interzonal_branches,
        )
    end
    return nothing
end

# ============================================================
# RUNNERS (ADMM-DC)
# ============================================================

# # ------------------------------------------------------------
# # Runner: ADMM bus/DC
# # ------------------------------------------------------------

function run_admm(
    sys::ADMMSystem,
    config::DistributedOPFModel,
    ::ADMM{Bus, DCTag},
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == DC "ADMM bus/DC requer base_model == DC."

    shared_buses = sys.shared_buses
    max_iter = config.max_iter

    rho = Ref(config.rho_theta)
    prev_rho = Ref(rho[])
    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer(config)

    partitions = build_dc_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(shared_buses)
    n_shared_buses = length(shared_buses_sorted)


    dual_theta = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )
    dual_theta_buf = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )

    theta_consensus = zeros(Float64, n_shared_buses)
    theta_consensus_prev = similar(theta_consensus)

    hist = ADMMHistory()
    hist.meta[:kind] = :dc_bus
    hist.meta[:tag] = "admm_dc_bus"

    theta_avg_history = Dict{Int, Vector{Float64}}()

    results = Vector{NamedTuple{(:gen_cost_expr, :real_cost), Tuple{Float64, Float64}}}(
        undef, n_zones,
    )

    init_dc_bus_partition!(partitions, sys, optimizer_factory, rho[])

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\nIteração $iter (bus/DC)")

        refresh_rho_bus_dc!(partitions, rho, prev_rho)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            theta_dual_buf_z = dual_theta_buf[zone_id]
            theta_dual_buf_z .= dual_theta[zone_id]

            lin = build_linear_bus_term(
                part,
                shared_buses_sorted,
                shared_buses_to_pos,
                theta_consensus,
                theta_dual_buf_z,
                rho[],
            )

            @objective(model, Min, part.obj_base + part.quad_bus_expr + lin)
            optimize!(model)

            results[i] = (
                gen_cost_expr = value(model[:generation_cost]),
                real_cost = compute_real_generation_cost(model, part.objective_params),
            )
        end

        @inbounds for (i, zone_id) in enumerate(zone_ids)
            part = get_partition(partitions, zone_id)
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost
        end

        theta_consensus_prev .= theta_consensus

        compute_theta_consensus_vec!(
            theta_consensus,
            partitions,
            shared_buses_sorted,
            shared_buses_to_pos,
        )

        @inbounds for b in shared_buses_sorted
            j = shared_buses_to_pos[b]
            push!(get!(theta_avg_history, b, Float64[]), theta_consensus[j])
        end

        update_duals_theta_all!(
            dual_theta,
            partitions,
            zone_ids,
            theta_consensus,
            rho,
            shared_buses_sorted,
            shared_buses_to_pos,
        )

        r = primal_residual_bus_dc(
            partitions,
            theta_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
        )

        s = dual_residual_bus_dc(theta_consensus, theta_consensus_prev, rho[])

        eps_r = eps_primal_bus_dc(
            partitions,
            theta_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
            eps_abs,
            eps_rel,
        )

        eps_s = eps_dual_bus_dc_unscaled(
            dual_theta, partitions, zone_ids, shared_buses_to_pos, eps_abs, eps_rel,
        )

        # eps_r = 1e-4

        # eps_s = 1e-4
        _record_iter!(hist, r, s, eps_r, eps_s; rhos = Dict(:rho_theta => rho[]))

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if stop_if_converged!(r, s, eps_r, eps_s)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        maybe_adapt_rho!(adaptive, rho, r, s, config)
    end

    hist.heatmaps[:theta_avg] = theta_avg_history

    result = _finalize_admm_result(partitions, t0_ns, iter_done)
    _print_admm_summary!(
        "[ADMM-bus/DC]",
        result;
        rho_lines = ("rho_theta inicial: $(config.rho_theta) | rho_theta final: $(rho[])",),
    )

    return LirionRunResult(; result, model = nothing, system = sys, diagnostics = hist)
end

# ------------------------------------------------------------
# Runner: ADMM bus/DC  (DEBUG pesado)
# ------------------------------------------------------------
# ------------------------------------------------------------
# Runner: ADMM bus/DC  (DEBUG EXTREMO em TODA iteração)
# ------------------------------------------------------------
# ---------- prints focados em "virtual vs real" (SEM sys_zone) ----------


# @inline function _count_fixed_theta(part::PartitionDC)
#     θ = part.theta_var  # JuMP container (provavelmente DenseAxisArray 2D: (bus, t))

#     nfix = 0
#     fixed_buses = Int[]

#     # eixo 1 costuma ser os bus IDs (Vector{Int})
#     ax1 = axes(θ, 1)
#     ax2 = axes(θ, 2)

#     # Caso 1: 2D (bus, t)
#     if length(ax2) > 0
#         for b in ax1
#             any_fixed = false
#             for t in ax2
#                 v = θ[b, t]
#                 if JuMP.is_fixed(v)
#                     any_fixed = true
#                     break
#                 end
#             end
#             if any_fixed
#                 nfix += 1
#                 push!(fixed_buses, Int(b))
#             end
#         end
#         return nfix, fixed_buses
#     end

#     # Caso 2: 1D (só bus) – fallback seguro
#     for b in ax1
#         v = θ[b]
#         if JuMP.is_fixed(v)
#             nfix += 1
#             push!(fixed_buses, Int(b))
#         end
#     end
#     return nfix, fixed_buses
# end

# @inline function _is_virtual_bus(part::PartitionDC, b::Int)::Union{Bool, Nothing}
#     # BusData tipicamente tem:
#     #   part.bus.id_to_idx :: Dict{Int,Int}
#     #   part.bus.is_virtual :: Vector{Bool}
#     bus = part.bus
#     if !(hasproperty(bus, :id_to_idx) && hasproperty(bus, :is_virtual))
#         return nothing
#     end
#     id2i = getproperty(bus, :id_to_idx)
#     isv  = getproperty(bus, :is_virtual)
#     if !haskey(id2i, b)
#         return nothing
#     end
#     idx = id2i[b]
#     if idx < 1 || idx > length(isv)
#         return nothing
#     end
#     return isv[idx]
# end

# @inline function _is_kcl_active(part::PartitionDC, b::Int)::Union{Bool, Nothing}
#     # No teu active_power_balance!, KCL só é aplicada em !is_virtual
#     v = _is_virtual_bus(part, b)
#     v === nothing && return nothing
#     return !v
# end

# function _print_shared_bus_virtual_table(partitions, zone_ids, shared_buses_sorted)
#     println("  [DEBUG] shared buses: status por zona (virtual/KCL):")
#     for b in shared_buses_sorted
#         line = "    bus=$b:"
#         for zone_id in zone_ids
#             part = get_partition(partitions, zone_id)
#             v  = _is_virtual_bus(part, b)
#             k  = _is_kcl_active(part, b)
#             if v === nothing || k === nothing
#                 line *= "  z$zone_id[virt=? kcl=?]"
#             else
#                 line *= "  z$zone_id[virt=$(v) kcl=$(k)]"
#             end
#         end
#         println(line)
#     end

#     # alerta: shared bus virtual em pelo menos uma zona
#     bad = Int[]
#     for b in shared_buses_sorted
#         any_virtual = false
#         any_unknown = false
#         for zone_id in zone_ids
#             part = get_partition(partitions, zone_id)
#             v = _is_virtual_bus(part, b)
#             if v === nothing
#                 any_unknown = true
#             elseif v
#                 any_virtual = true
#             end
#         end
#         if any_virtual || any_unknown
#             push!(bad, b)
#         end
#     end

#     if !isempty(bad)
#         println("  [ALERTA] shared buses suspeitas (virtual em alguma zona ou status desconhecido):")
#         println("           " * join(bad, ", "))
#     else
#         println("  [OK] nenhum shared bus parece virtual nas zonas (via part.bus.is_virtual).")
#     end
# end

# function run_admm(
#     sys::ADMMSystem,
#     config::DistributedOPFModel,
#     ::ADMM{Bus, DCTag},
# )::LirionRunResult
#     t0_ns = time_ns()
#     @assert config.base_model == DC "ADMM bus/DC requer base_model == DC."

#     shared_buses = sys.shared_buses
#     max_iter = config.max_iter

#     rho = Ref(config.rho_theta)
#     prev_rho = Ref(rho[])
#     eps_abs = config.eps_abs
#     eps_rel = config.eps_rel
#     adaptive = _adaptive_on(config)

#     optimizer_factory = resolve_optimizer(config)

#     partitions = build_dc_partitions(sys)
#     zone_ids = partitions.zones
#     n_zones = length(zone_ids)

#     shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(shared_buses)
#     n_shared_buses = length(shared_buses_sorted)

#     @show shared_buses
#     # ss
#     dual_theta = Dict{Int, Vector{Float64}}(
#         zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
#     )
#     dual_theta_buf = Dict{Int, Vector{Float64}}(
#         zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
#     )

#     theta_consensus = zeros(Float64, n_shared_buses)
#     theta_consensus_prev = similar(theta_consensus)

#     hist = ADMMHistory()
#     hist.meta[:kind] = :dc_bus
#     hist.meta[:tag]  = "admm_dc_bus"

#     theta_avg_history = Dict{Int, Vector{Float64}}()

#     results = Vector{NamedTuple{(:gen_cost_expr, :real_cost), Tuple{Float64, Float64}}}(
#         undef, n_zones,
#     )

#     # ============================================================
#     # Debug helpers (100% locais ao runner)  -- SEM sys_zone
#     # ============================================================

#     @inline function _θ_shared_vec(part::PartitionDC)
#         return snapshot_theta_shared!(part)  # θ nos shared (ordem part.shared_buses_zone)
#     end

#     @inline function _zone_gauge_stats(part::PartitionDC)
#         θs = _θ_shared_vec(part)
#         if isempty(θs)
#             return (mean = 0.0, span = 0.0, minv = 0.0, maxv = 0.0, n = 0)
#         end
#         mn = minimum(θs)
#         mx = maximum(θs)
#         return (mean = sum(θs) / length(θs), span = mx - mn, minv = mn, maxv = mx, n = length(θs))
#     end

#     @inline function _zone_primal_r(part::PartitionDC, z::Vector{Float64})
#         θs = _θ_shared_vec(part)
#         ids = part.shared_buses_zone
#         acc = 0.0
#         @inbounds for i in eachindex(ids)
#             b = ids[i]
#             j = shared_buses_to_pos[b]
#             d = θs[i] - z[j]
#             acc += d*d
#         end
#         return sqrt(acc)
#     end

#     @inline function _dual_norms(dual_theta::Dict{Int, Vector{Float64}})
#         maxn = 0.0
#         sumn = 0.0
#         for z in zone_ids
#             v = dual_theta[z]
#             nrm = sqrt(sum(vv*vv for vv in v))
#             maxn = max(maxn, nrm)
#             sumn += nrm
#         end
#         return maxn, sumn / max(1, length(zone_ids))
#     end

#     @inline function _safe_minmax(v::Vector{Float64})
#         isempty(v) && return (0.0, 0.0)
#         return (minimum(v), maximum(v))
#     end

#     @inline function _topk_abs_entries(v::Vector{Float64}, k::Int)
#         n = length(v)
#         k = min(k, n)
#         idxs = collect(1:n)
#         sort!(idxs, by = i -> abs(v[i]), rev = true)
#         return [(idxs[t], v[idxs[t]]) for t in 1:k]
#     end

#     @inline function _is_virtual_bus(part::PartitionDC, b::Int)::Union{Bool, Nothing}
#         bus = part.bus
#         if !(hasproperty(bus, :id_to_idx) && hasproperty(bus, :is_virtual))
#             return nothing
#         end
#         id2i = getproperty(bus, :id_to_idx)
#         isv  = getproperty(bus, :is_virtual)
#         if !haskey(id2i, b)
#             return nothing
#         end
#         idx = id2i[b]
#         if idx < 1 || idx > length(isv)
#             return nothing
#         end
#         return isv[idx]
#     end

#     @inline function _is_kcl_active(part::PartitionDC, b::Int)::Union{Bool, Nothing}
#         v = _is_virtual_bus(part, b)
#         v === nothing && return nothing
#         return !v
#     end

#     function _print_shared_bus_virtual_table(ps, zone_ids, shared_buses_sorted)
#         println("  [DEBUG] shared buses: status por zona (virtual/KCL):")
#         for b in shared_buses_sorted
#             line = "    bus=$b:"
#             for zone_id in zone_ids
#                 part = get_partition(ps, zone_id)
#                 v  = _is_virtual_bus(part, b)
#                 k  = _is_kcl_active(part, b)
#                 if v === nothing || k === nothing
#                     line *= "  z$zone_id[virt=? kcl=?]"
#                 else
#                     line *= "  z$zone_id[virt=$(v) kcl=$(k)]"
#                 end
#             end
#             println(line)
#         end

#         bad = Int[]
#         for b in shared_buses_sorted
#             any_virtual = false
#             any_unknown = false
#             for zone_id in zone_ids
#                 part = get_partition(ps, zone_id)
#                 v = _is_virtual_bus(part, b)
#                 if v === nothing
#                     any_unknown = true
#                 elseif v
#                     any_virtual = true
#                 end
#             end
#             if any_virtual || any_unknown
#                 push!(bad, b)
#             end
#         end

#         if !isempty(bad)
#             println("  [ALERTA] shared buses suspeitas (virtual em alguma zona ou status desconhecido):")
#             println("           " * join(bad, ", "))
#         else
#             println("  [OK] nenhum shared bus parece virtual nas zonas (via part.bus.is_virtual).")
#         end
#     end

#     # ============================================================
#     # INIT
#     # ============================================================

#     println("\n================ [ADMM-bus/DC] DEBUG INIT ================")
#     println("zones = $(zone_ids) | n_shared=$(n_shared_buses)")
#     println("shared_buses_sorted = $(shared_buses_sorted)")
#     println("  [DEBUG] shared buses: status por zona (virtual/KCL):")
#     _print_shared_bus_virtual_table(partitions, zone_ids, shared_buses_sorted)

#     init_dc_bus_partition!(partitions, sys, optimizer_factory, rho[])

#     iter_done = 0
#     for iter in 1:max_iter
#         iter_done = iter

#         rho_used_before = rho[]

#         println("\n================ [ADMM-bus/DC] ITER $iter ================")
#         println("rho_theta (start iter) = $(rho_used_before) | adaptive=$(adaptive)")

#         # rebuild do termo quadrático se rho mudou
#         refresh_rho_bus_dc!(partitions, rho, prev_rho)

#         if rho[] != rho_used_before
#             println("  [DEBUG] refresh_rho mudou rho antes dos solves: $rho_used_before -> $(rho[])")
#         end

#         rho_used = rho[]
#         println("rho_theta (used in solves) = $(rho_used)")

#         # (print periódico pra não poluir demais)
#         if iter == 1 || (iter % 250 == 0)
#             println("  [DEBUG] shared buses: status por zona (virtual/KCL) [periodic]:")
#             _print_shared_bus_virtual_table(partitions, zone_ids, shared_buses_sorted)
#         end

#         # ============================================================
#         # Solve local por zona
#         # ============================================================
#         Threads.@threads for i in 1:n_zones
#             zone_id = zone_ids[i]
#             part = get_partition(partitions, zone_id)
#             model = part.model::JuMP.Model

#             theta_dual_buf_z = dual_theta_buf[zone_id]
#             theta_dual_buf_z .= dual_theta[zone_id]

#             lin = build_linear_bus_term(
#                 part,
#                 shared_buses_sorted,
#                 shared_buses_to_pos,
#                 theta_consensus,
#                 theta_dual_buf_z,
#                 rho_used,
#             )

#             @objective(model, Min, part.obj_base + part.quad_bus_expr + lin)
            
#             optimize!(model)

#             # print(model)
#             if iter == 1
#                 println("  [DEBUG] θ fixed (âncoras) por zona:")
#                 for zone_id in zone_ids
#                     part = get_partition(partitions, zone_id)
#                     nfix, buses = _count_fixed_theta(part)
#                     println("    zone=$zone_id  n_fixed=$nfix  buses_fixed=$(buses)")
#                 end
#             end


#             _require_optimal!(model; iter=iter, zone_id=zone_id, tag="ADMM-bus/DC")

#             gc = try value(model[:generation_cost]) catch; NaN end
#             rc = try compute_real_generation_cost(model, part.objective_params) catch; NaN end

#             results[i] = (gen_cost_expr = gc, real_cost = rc)
#         end

#         @inbounds for (i, zone_id) in enumerate(zone_ids)
#             part = get_partition(partitions, zone_id)
#             out = results[i]
#             part.has_solution = true
#             part.generation_cost_model_expr = out.gen_cost_expr
#             part.real_generation_cost = out.real_cost
#         end

#         # custos por zona
#         println("  [DEBUG] custos por zona:")
#         @inbounds for (i, zone_id) in enumerate(zone_ids)
#             out = results[i]
#             println("    zone=$zone_id  gen_cost_expr=$(out.gen_cost_expr)  real_cost=$(out.real_cost)")
#         end

#         # ============================================================
#         # Consensus / duals / residuals
#         # ============================================================
#         theta_consensus_prev .= theta_consensus

#         compute_theta_consensus_vec!(
#             theta_consensus,
#             partitions,
#             shared_buses_sorted,
#             shared_buses_to_pos,
#         )

#         zmin, zmax = _safe_minmax(theta_consensus)
#         println("  [DEBUG] consenso z(theta): min=$zmin max=$zmax  range=$(zmax - zmin)  n_shared=$n_shared_buses")

#         # histórico por bus
#         @inbounds for b in shared_buses_sorted
#             j = shared_buses_to_pos[b]
#             push!(get!(theta_avg_history, b, Float64[]), theta_consensus[j])
#         end

#         # gauge stats por zona (SEMPRE)
#         println("  [DEBUG] gauge stats (θ_shared) por zona:")
#         for zone_id in zone_ids
#             part = get_partition(partitions, zone_id)
#             st = _zone_gauge_stats(part)
#             println("    zone=$zone_id shared_n=$(st.n) meanθ=$(st.mean) span=$(st.span) [min=$(st.minv), max=$(st.maxv)]")
#         end

#         # update duals
#         update_duals_theta_all!(
#             dual_theta,
#             partitions,
#             zone_ids,
#             theta_consensus,
#             rho,
#             shared_buses_sorted,
#             shared_buses_to_pos,
#         )

#         # residual primal
#         r = primal_residual_bus_dc(
#             partitions,
#             theta_consensus,
#             shared_buses_sorted,
#             shared_buses_to_pos,
#         )

#         # Boyd: s^k = rho * ||z^k - z^{k-1}||
#         dz_norm = sqrt(sum((theta_consensus[j] - theta_consensus_prev[j])^2 for j in eachindex(theta_consensus)))
#         s = rho_used * dz_norm

#         eps_r = eps_primal_bus_dc(
#             partitions,
#             theta_consensus,
#             shared_buses_sorted,
#             shared_buses_to_pos,
#             eps_abs,
#             eps_rel,
#         )

#         eps_s = eps_dual_bus_dc_unscaled(
#             dual_theta, partitions, zone_ids, shared_buses_to_pos, eps_abs, eps_rel,
#         )

#         println("  [DEBUG] dz_norm = $dz_norm  => s = rho*dz_norm = $(rho_used) * $dz_norm = $s")
#         println("  [DEBUG] residuals: ‖r‖=$r (eps_pri=$eps_r) | ‖s‖=$s (eps_du=$eps_s)")

#         # primal por zona + dual norms por zona (SEMPRE)
#         max_pri = -Inf
#         argmax_z = zone_ids[1]
#         println("  [DEBUG] r_z (por zona):")
#         for zone_id in zone_ids
#             part = get_partition(partitions, zone_id)
#             rz = _zone_primal_r(part, theta_consensus)
#             println("    zone=$zone_id  ‖r_z‖=$rz")
#             if rz > max_pri
#                 max_pri = rz
#                 argmax_z = zone_id
#             end
#         end
#         max_du, mean_du = _dual_norms(dual_theta)
#         println("  [DEBUG] pior zona (primal): zone=$argmax_z  max‖r_z‖=$max_pri")
#         println("  [DEBUG] dual norms: max||y_z||=$max_du  mean||y_z||=$mean_du")

#         # amostra dos maiores y
#         println("  [DEBUG] top-|y| por zona (top 5 entradas):")
#         for zone_id in zone_ids
#             v = dual_theta[zone_id]
#             top = _topk_abs_entries(v, 5)
#             msg = ""
#             for (idx, val) in top
#                 bus = shared_buses_sorted[idx]
#                 msg *= " (bus=$bus idx=$idx y=$(val))"
#             end
#             println("    zone=$zone_id:$msg")
#         end

#         # sanidade
#         if !(isfinite(r) && isfinite(s) && isfinite(eps_r) && isfinite(eps_s))
#             error("""
#             [ADMM-bus/DC FAIL-FAST] NaN/Inf detectado
#               iter=$iter
#               r=$r s=$s eps_r=$eps_r eps_s=$eps_s
#               rho=$rho_used
#               consenso[min,max]=[$zmin,$zmax]
#             """)
#         end

#         _record_iter!(hist, r, s, eps_r, eps_s; rhos = Dict(:rho_theta => rho_used))

#         if stop_if_converged!(r, s, eps_r, eps_s)
#             println("  [STOP] critério atingido em iter=$iter")
#             break
#         end

#         # Por que o adaptativo mudaria rho
#         if adaptive
#             if r > config.mu_incr * s
#                 println("  [DEBUG] adaptive_rho trigger: r > mu_incr*s  => tende a AUMENTAR rho")
#             elseif s > config.mu_decr * r
#                 println("  [DEBUG] adaptive_rho trigger: s > mu_decr*r  => tende a DIMINUIR rho")
#             else
#                 println("  [DEBUG] adaptive_rho: dentro da banda (sem pressão forte)")
#             end
#         end

#         rho_before = rho[]
#         maybe_adapt_rho!(adaptive, rho, r, s, config)
#         if rho[] != rho_before
#             println("  [DEBUG] rho atualizado no fim da iter: $rho_before -> $(rho[])")
#         end
#     end

#     hist.heatmaps[:theta_avg] = theta_avg_history

#     result = _finalize_admm_result(partitions, t0_ns, iter_done)
#     _print_admm_summary!(
#         "[ADMM-bus/DC]",
#         result;
#         rho_lines = ("rho_theta inicial: $(config.rho_theta) | rho_theta final: $(rho[])",),
#     )

#     return LirionRunResult(; result, model = nothing, system = sys, diagnostics = hist)
# end

















































































# ------------------------------------------------------------
# Runner: ADMM branch/DC
# ------------------------------------------------------------

function run_admm(
    sys::ADMMSystem,
    config::DistributedOPFModel,
    ::ADMM{Branch, DCTag},
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == DC "ADMM-branch/DC requer base_model == DC."

    interzonal_branches = sys.interzonal_branches
    max_iter = config.max_iter

    rho = Ref(config.rho_flow_p)
    prev_rho = Ref(rho[])
    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer(config)

    partitions = build_dc_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    flow_consensus_prev = Dict{Int, Float64}()

    duals_flow = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )

    hist = ADMMHistory()
    hist.meta[:kind] = :dc_branch
    hist.meta[:tag] = "admm_dc_branch"

    results = Vector{NamedTuple{(:gen_cost_expr, :real_cost), Tuple{Float64, Float64}}}(
        undef, n_zones,
    )

    init_dc_branch_partition!(
        partitions,
        sys,
        optimizer_factory,
        rho[],
        interzonal_branches,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-branch/DC] Iteração $iter")

        refresh_rho_branch_dc!(partitions, rho, prev_rho)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            lin = build_linear_branch_term(
                part,
                part.zone_edges,
                flow_consensus,
                duals_flow[zone_id],
                rho[],
            )

            @objective(model, Min, part.obj_base + part.quad_branch_expr + lin)
            optimize!(model)
            _require_optimal!(model; iter=iter, zone_id=zone_id, tag="ADMM-branch/DC")

            results[i] = (
                gen_cost_expr = value(model[:generation_cost]),
                real_cost = compute_real_generation_cost(model, part.objective_params),
            )
        end

        @inbounds for (i, zone_id) in enumerate(zone_ids)
            part = get_partition(partitions, zone_id)
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost
        end

        flow_consensus_prev = copy(flow_consensus)

        compute_flow_consensus!(flow_consensus, partitions, interzonal_branches)

        update_duals_flow!(
            duals_flow,
            partitions,
            flow_consensus,
            interzonal_branches,
            rho[],
        )

        r = primal_residual_branch_dc(partitions, flow_consensus, interzonal_branches)

        s = dual_residual_branch_dc(flow_consensus, flow_consensus_prev, rho[])

        eps_r = eps_primal_branch_dc(
            partitions,
            flow_consensus,
            interzonal_branches,
            eps_abs,
            eps_rel,
        )

        eps_s = eps_dual_branch_dc_unscaled(
            duals_flow,
            partitions,
            zone_ids,
            interzonal_branches,
            eps_abs,
            eps_rel,
        )
        

        _record_iter!(hist, r, s, eps_r, eps_s; rhos = Dict(:rho_flow => rho[]))

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if stop_if_converged!(r, s, eps_r, eps_s)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        maybe_adapt_rho!(adaptive, rho, r, s, config)
    end

    result = _finalize_admm_result(partitions, t0_ns, iter_done)
    _print_admm_summary!(
        "[ADMM-branch/DC]",
        result;
        rho_lines = ("rho_flow inicial: $(config.rho_flow_p) | rho_flow final: $(rho[])",),
    )

    return LirionRunResult(; result, model = nothing, system = sys, diagnostics = hist)
end

# ------------------------------------------------------------
# Runner: ADMM hybrid/DC
# ------------------------------------------------------------

function run_admm(
    sys::ADMMSystem,
    config::DistributedOPFModel,
    ::ADMM{Hybrid, DCTag},
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == DC "Hybrid ADMM/DC requer DC."

    shared_buses = sys.shared_buses
    interzonal_branches = sys.interzonal_branches
    max_iter = config.max_iter

    rho_theta = Ref(config.rho_theta)
    rho_flow = Ref(config.rho_flow_p)

    prev_rho_theta = Ref(rho_theta[])
    prev_rho_flow = Ref(rho_flow[])

    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer(config)

    partitions = build_dc_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    theta_consensus = zeros(Float64, n_shared_buses)
    theta_consensus_prev = similar(theta_consensus)

    flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    flow_consensus_prev = Dict{Int, Float64}()

    duals_theta = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )
    dual_theta_buf = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )
    duals_flow = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )

    hist = ADMMHistory()
    hist.meta[:kind] = :dc_hybrid
    hist.meta[:tag] = "admm_dc_hybrid"

    theta_avg_history = Dict{Int, Vector{Float64}}()

    results = Vector{NamedTuple{(:gen_cost_expr, :real_cost), Tuple{Float64, Float64}}}(
        undef, n_zones,
    )

    init_dc_hybrid_partition!(
        partitions,
        sys,
        optimizer_factory,
        rho_theta[],
        rho_flow[],
        interzonal_branches,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-hybrid/DC] Iteração $iter")

        refresh_rho_bus_dc!(partitions, rho_theta, prev_rho_theta)
        refresh_rho_branch_dc!(partitions, rho_flow, prev_rho_flow)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            theta_dual_buf_z = dual_theta_buf[zone_id]
            theta_dual_buf_z .= duals_theta[zone_id]

            lin_bus = build_linear_bus_term(
                part,
                shared_buses_sorted,
                shared_buses_to_pos,
                theta_consensus,
                theta_dual_buf_z,
                rho_theta[],
            )
            lin_branch = build_linear_branch_term(
                part,
                part.zone_edges,
                flow_consensus,
                duals_flow[zone_id],
                rho_flow[],
            )

            @objective(
                model,
                Min,
                part.obj_base + part.quad_bus_expr + part.quad_branch_expr + lin_bus +
                lin_branch
            )
            optimize!(model)
            _require_optimal!(model; iter=iter, zone_id=zone_id, tag="ADMM-hybrid/DC")

            results[i] = (
                gen_cost_expr = value(model[:generation_cost]),
                real_cost = compute_real_generation_cost(model, part.objective_params),
            )
        end

        @inbounds for (i, zone_id) in enumerate(zone_ids)
            part = get_partition(partitions, zone_id)
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost
        end

        theta_consensus_prev .= theta_consensus
        flow_consensus_prev = copy(flow_consensus)

        compute_theta_consensus_vec!(
            theta_consensus,
            partitions,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        @inbounds for b in shared_buses_sorted
            j = shared_buses_to_pos[b]
            push!(get!(theta_avg_history, b, Float64[]), theta_consensus[j])
        end

        compute_flow_consensus!(flow_consensus, partitions, interzonal_branches)

        update_duals_theta_all!(
            duals_theta,
            partitions,
            zone_ids,
            theta_consensus,
            rho_theta,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        update_duals_flow!(
            duals_flow,
            partitions,
            flow_consensus,
            interzonal_branches,
            rho_flow[],
        )

        r_bus = primal_residual_bus_dc(
            partitions,
            theta_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        r_br = primal_residual_branch_dc(partitions, flow_consensus, interzonal_branches)
        r_tot = hypot(r_bus, r_br)

        s_bus = dual_residual_bus_dc(theta_consensus, theta_consensus_prev, rho_theta[])
        s_br = dual_residual_branch_dc(flow_consensus, flow_consensus_prev, rho_flow[])
        s_tot = hypot(s_bus, s_br)

        eps_bus_r = eps_primal_bus_dc(
            partitions,
            theta_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
            eps_abs,
            eps_rel,
        )
        eps_br_r = eps_primal_branch_dc(
            partitions,
            flow_consensus,
            interzonal_branches,
            eps_abs,
            eps_rel,
        )
        eps_r_tot = hypot(eps_bus_r, eps_br_r)

        eps_bus_s = eps_dual_bus_dc_unscaled(
            duals_theta, partitions, zone_ids, shared_buses_to_pos, eps_abs, eps_rel,
        )
        eps_br_s = eps_dual_branch_dc_unscaled(
            duals_flow, partitions, zone_ids, interzonal_branches, eps_abs, eps_rel,
        )
        eps_s_tot = hypot(eps_bus_s, eps_br_s)

        _record_iter!(
            hist, r_tot, s_tot, eps_r_tot, eps_s_tot;
            rhos = Dict(:rho_theta => rho_theta[], :rho_flow => rho_flow[]),
        )

        _FAST_USE_LOGS && println("‖r‖ = $r_tot | eps_pri = $eps_r_tot")
        _FAST_USE_LOGS && println("‖s‖ = $s_tot | eps_du  = $eps_s_tot")

        if stop_if_converged!(r_tot, s_tot, eps_r_tot, eps_s_tot)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        maybe_adapt_rho!(
            adaptive,
            rho_theta, r_bus, s_bus,
            rho_flow, r_br, s_br,
            config,
        )
    end

    hist.heatmaps[:theta_avg] = theta_avg_history

    result = _finalize_admm_result(partitions, t0_ns, iter_done)
    _print_admm_summary!(
        "[ADMM-hybrid/DC]",
        result;
        rho_lines = (
            "rho_theta inicial: $(config.rho_theta) | rho_theta final: $(rho_theta[])",
            "rho_flow  inicial: $(config.rho_flow_p) | rho_flow  final: $(rho_flow[])",
        ),
    )

    return LirionRunResult(; result, model = nothing, system = sys, diagnostics = hist)
end





# ------------------------------------------------------------
# Runner: ADMM branch/DC  (COM DEBUGS DETALHADOS)
# ------------------------------------------------------------
# O que este runner imprime (quando debug=true):
#  - rho antes/depois do refresh e se mudou
#  - status do solver por zona + objective
#  - para alguns ramos interzonais: (p_z por zona), z(consenso), delta, sigma-suspeito (sinais opostos)
#  - norma e range dos duals por zona
#  - r, s, eps_r, eps_s
#
# Ajustes rápidos:
#   debug=true/false
#   dbg_every = imprime a cada dbg_every iterações (1 = todas)
#   dbg_edges_max = quantos ramos interzonais inspecionar
#   dbg_zones_max = quantas zonas imprimir na tabela
#   dbg_show_edges = lista opcional de ramos globais para sempre inspecionar (mesmo se não estiver entre os maiores)
# ------------------------------------------------------------

# function run_admm(
#     sys::ADMMSystem,
#     config::DistributedOPFModel,
#     ::ADMM{Branch, DCTag},
# )::LirionRunResult
#     t0_ns = time_ns()
#     @assert config.base_model == DC "ADMM-branch/DC requer base_model == DC."

#     interzonal_branches = sys.interzonal_branches
#     max_iter = config.max_iter

#     rho = Ref(config.rho_flow_p)
#     prev_rho = Ref(rho[])
#     eps_abs = config.eps_abs
#     eps_rel = config.eps_rel
#     adaptive = _adaptive_on(config)

#     optimizer_factory = resolve_optimizer(config)

#     partitions = build_dc_partitions(sys)
#     zone_ids = partitions.zones
#     n_zones = length(zone_ids)

#     flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
#     flow_consensus_prev = Dict{Int, Float64}()

#     duals_flow = Dict{Int, Dict{Int, Float64}}(
#         zone_id => Dict{Int, Float64}() for zone_id in zone_ids
#     )

#     hist = ADMMHistory()
#     hist.meta[:kind] = :dc_branch
#     hist.meta[:tag] = "admm_dc_branch"

#     results = Vector{NamedTuple{(:gen_cost_expr, :real_cost), Tuple{Float64, Float64}}}(
#         undef, n_zones,
#     )

#     init_dc_branch_partition!(
#         partitions,
#         sys,
#         optimizer_factory,
#         rho[],
#         interzonal_branches,
#     )

#     # --------------------------
#     # DEBUG knobs (altere aqui)
#     # --------------------------
#     debug::Bool = true
#     dbg_every::Int = 1
#     dbg_edges_max::Int = 8
#     dbg_zones_max::Int = 6
#     dbg_top_by_abs_p::Bool = true
#     dbg_show_edges::Vector{Int} = Int[]

#     # -----------------------------------------
#     # FAIL-FAST knobs (PARE NA HORA)
#     # -----------------------------------------
#     fail_fast::Bool = true
#     # pare se term != OPTIMAL/LOCALLY_SOLVED OU primal != FEASIBLE_POINT
#     strict_primal::Bool = true

#     # helpers locais
#     @inline function _safe_val(x, default::Float64=0.0)
#         return (x === nothing || !isfinite(x)) ? default : Float64(x)
#     end

#     @inline function _dc_flow_local(m::JuMP.Model, part::PartitionDC, e::Int)
#         l_local = get(part.branch_global_to_local, e, 0)
#         l_local == 0 && return nothing
#         return _dc_flow_val_model(m, part, l_local)
#     end

#     @inline function _ok_status(term, prim)::Bool
#         ok_term = (term == MOI.OPTIMAL) || (term == MOI.LOCALLY_SOLVED)
#         if !ok_term
#             return false
#         end
#         if strict_primal
#             return prim == MOI.FEASIBLE_POINT
#         else
#             # aceita UNKNOWN_RESULT_STATUS quando term é ok (às vezes acontece com alguns solvers)
#             return (prim == MOI.FEASIBLE_POINT) || (prim == MOI.UNKNOWN_RESULT_STATUS)
#         end
#     end

#     @inline function _mk_fail_msg(iter::Int, zone_id::Int, model::JuMP.Model)
#         term = try JuMP.termination_status(model) catch; nothing end
#         prim = try JuMP.primal_status(model) catch; nothing end
#         dual = try JuMP.dual_status(model) catch; nothing end
#         objv = try JuMP.objective_value(model) catch; NaN end
#         raw  = try JuMP.raw_status(model) catch; "" end
#         return """
#         [ADMM-branch/DC FAIL-FAST]
#         iter=$iter zone=$zone_id
#           termination_status = $term
#           primal_status       = $prim
#           dual_status         = $dual
#           objective_value     = $objv
#           raw_status          = $raw
#         """
#     end

#     function _dbg_pick_edges(
#         partitions::Partitions{PartitionDC},
#         zone_ids::Vector{Int},
#         interzonal_branches::Vector{Int},
#         dbg_edges_max::Int,
#         dbg_show_edges::Vector{Int},
#         top_by_abs_p::Bool,
#     )::Vector{Int}
#         chosen = Set{Int}()
#         for e in dbg_show_edges
#             if e in interzonal_branches
#                 push!(chosen, e)
#             end
#         end
#         if length(chosen) >= dbg_edges_max
#             return collect(chosen)[1:dbg_edges_max]
#         end
#         if !top_by_abs_p
#             for e in interzonal_branches
#                 push!(chosen, e)
#                 length(chosen) >= dbg_edges_max && break
#             end
#             return collect(chosen)
#         end
#         scores = Vector{Tuple{Float64,Int}}()
#         for e in interzonal_branches
#             acc = 0.0
#             cnt = 0
#             for z in zone_ids
#                 part = get_partition(partitions, z)
#                 m = part.model::JuMP.Model
#                 p = _dc_flow_local(m, part, e)
#                 p === nothing && continue
#                 acc += abs(p)
#                 cnt += 1
#             end
#             cnt == 0 && continue
#             push!(scores, (acc / cnt, e))
#         end
#         sort!(scores, by = x -> x[1], rev = true)
#         for (_, e) in scores
#             push!(chosen, e)
#             length(chosen) >= dbg_edges_max && break
#         end
#         return collect(chosen)
#     end

#     function _dbg_print_iter!(
#         iter::Int,
#         partitions::Partitions{PartitionDC},
#         zone_ids::Vector{Int},
#         interzonal_branches::Vector{Int},
#         flow_consensus::Dict{Int,Float64},
#         duals_flow::Dict{Int,Dict{Int,Float64}},
#         rho::Float64,
#         prev_rho_val::Float64,
#         r::Float64,
#         s::Float64,
#         eps_r::Float64,
#         eps_s::Float64,
#         results,
#     )
#         println("\n================ DEBUG ADMM-branch/DC ================")
#         println("iter = $iter")
#         println("rho = $rho | prev_rho = $prev_rho_val | changed = $(rho != prev_rho_val)")
#         println("r = $r (eps_r=$eps_r) | s = $s (eps_s=$eps_s)")

#         println("\n-- Zones (até $dbg_zones_max) --")
#         zshow = zone_ids[1:min(length(zone_ids), dbg_zones_max)]
#         for (k, z) in enumerate(zshow)
#             part = get_partition(partitions, z)
#             m = part.model::JuMP.Model
#             term = try JuMP.termination_status(m) catch; nothing end
#             prim = try JuMP.primal_status(m) catch; nothing end
#             objv = try JuMP.objective_value(m) catch; NaN end
#             gc   = try results[k].real_cost catch; NaN end
#             println("z=$z | term=$term | primal=$prim | obj=$objv | real_cost=$gc | #edges=$(length(part.zone_edges))")
#         end
#         if length(zone_ids) > length(zshow)
#             println("... ($(length(zone_ids)-length(zshow)) zonas omitidas)")
#         end

#         println("\n-- Dual stats (y) por zona: norm2 / min / max / nnz --")
#         for z in zshow
#             y = duals_flow[z]
#             if isempty(y)
#                 println("z=$z | empty")
#                 continue
#             end
#             acc = 0.0
#             y_min = Inf
#             y_max = -Inf
#             nnz = 0
#             for (_, v) in y
#                 if isfinite(v)
#                     acc += v*v
#                     y_min = min(y_min, v)
#                     y_max = max(y_max, v)
#                     nnz += 1
#                 end
#             end
#             println("z=$z | ‖y‖2=$(sqrt(acc)) | min=$y_min | max=$y_max | nnz=$nnz")
#         end

#         edges_to_check = _dbg_pick_edges(partitions, zone_ids, interzonal_branches,
#                                          dbg_edges_max, dbg_show_edges, dbg_top_by_abs_p)

#         println("\n-- Edges debug (até $dbg_edges_max) --")
#         println("Formato: e | z(consensus) | p_z por zona | (p_z - z) | sinais")
#         for e in edges_to_check
#             zc = get(flow_consensus, e, 0.0)
#             ps = Vector{Tuple{Int,Float64}}()
#             signs = Int[]
#             for z in zone_ids
#                 part = get_partition(partitions, z)
#                 m = part.model::JuMP.Model
#                 p = _dc_flow_local(m, part, e)
#                 p === nothing && continue
#                 push!(ps, (z, p))
#                 push!(signs, p > 0 ? 1 : (p < 0 ? -1 : 0))
#             end
#             has_pos = any(s -> s == 1, signs)
#             has_neg = any(s -> s == -1, signs)
#             sig_warn = (has_pos && has_neg) ? "⚠️(+/-)" : ""

#             buf = IOBuffer()
#             print(buf, "e=$e | z=$zc | ")
#             showk = min(length(ps), 6)
#             for i in 1:showk
#                 (zz, pv) = ps[i]
#                 dv = pv - zc
#                 print(buf, "z$zz:p=$(pv),Δ=$(dv)  ")
#             end
#             if length(ps) > showk
#                 print(buf, "...($(length(ps)-showk) mais) ")
#             end
#             print(buf, sig_warn)
#             println(String(take!(buf)))
#         end

#         println("======================================================\n")
#         return nothing
#     end

#     # --------------------------
#     # Thread-safe fail-fast
#     # --------------------------
#     fail_lock = Threads.SpinLock()
#     fail_flag = Ref(false)
#     fail_err  = Ref{Any}(nothing)

#     iter_done = 0
#     for iter in 1:max_iter
#         iter_done = iter
#         _FAST_USE_LOGS && println("\n[ADMM-branch/DC] Iteração $iter")

#         # reset fail-fast refs por iteração
#         fail_flag[] = false
#         fail_err[] = nothing

#         rho_before_refresh = rho[]
#         prev_rho_before = prev_rho[]
#         refresh_rho_branch_dc!(partitions, rho, prev_rho)

#         if debug && (iter % dbg_every == 0) && (rho[] != rho_before_refresh)
#             println("[DBG] rho mudou dentro de refresh_rho_branch_dc!: $rho_before_refresh -> $(rho[])")
#         end

#         Threads.@threads for i in 1:n_zones
#             # Se já falhou em outra thread, evita trabalho extra
#             if fail_fast && fail_flag[]
#                 continue
#             end

#             zone_id = zone_ids[i]
#             part = get_partition(partitions, zone_id)
#             model = part.model::JuMP.Model

#             try
#                 lin = build_linear_branch_term(
#                     part,
#                     part.zone_edges,
#                     flow_consensus,
#                     duals_flow[zone_id],
#                     rho[],
#                 )

#                 @objective(model, Min, part.obj_base + part.quad_branch_expr + lin)
#                 optimize!(model)

#                 term = JuMP.termination_status(model)
#                 prim = JuMP.primal_status(model)
#                 @info "OPTIMIZER" optimizer = JuMP.solver_name(model)

#                 if debug && (iter % dbg_every == 0)
#                     if !_ok_status(term, prim)
#                         println("[DBG] iter=$iter zone=$zone_id term=$term primal=$prim")
#                     end
#                 end

#                 if fail_fast && !_ok_status(term, prim)
#                     msg = _mk_fail_msg(iter, zone_id, model)
#                     lock(fail_lock)
#                     try
#                         if !fail_flag[]
#                             fail_flag[] = true
#                             fail_err[] = ErrorException(msg)
#                         end
#                     finally
#                         unlock(fail_lock)
#                     end
#                     continue
#                 end

#                 results[i] = (
#                     gen_cost_expr = value(model[:generation_cost]),
#                     real_cost = compute_real_generation_cost(model, part.objective_params),
#                 )
#             catch err
#                 if fail_fast
#                     lock(fail_lock)
#                     try
#                         if !fail_flag[]
#                             fail_flag[] = true
#                             fail_err[] = err
#                         end
#                     finally
#                         unlock(fail_lock)
#                     end
#                 else
#                     rethrow()
#                 end
#             end
#         end

#         # se alguma thread marcou falha, para IMEDIATAMENTE aqui
#         if fail_fast && fail_flag[]
#             # (opcional) imprime debug do estado do consenso antes de sair
#             if debug
#                 println("\n[DBG] Parando no iter=$iter por falha de zona.")
#             end
#             throw(fail_err[] isa Exception ? fail_err[] : ErrorException(string(fail_err[])))
#         end

#         @inbounds for (i, zone_id) in enumerate(zone_ids)
#             part = get_partition(partitions, zone_id)
#             out = results[i]
#             part.has_solution = true
#             part.generation_cost_model_expr = out.gen_cost_expr
#             part.real_generation_cost = out.real_cost
#         end

#         flow_consensus_prev = copy(flow_consensus)
#         compute_flow_consensus!(flow_consensus, partitions, interzonal_branches)

#         update_duals_flow!(
#             duals_flow,
#             partitions,
#             flow_consensus,
#             interzonal_branches,
#             rho[],
#         )

#         r = primal_residual_branch_dc(partitions, flow_consensus, interzonal_branches)
#         s = dual_residual_branch_dc(flow_consensus, flow_consensus_prev, rho[])

#         eps_r = eps_primal_branch_dc(
#             partitions,
#             flow_consensus,
#             interzonal_branches,
#             eps_abs,
#             eps_rel,
#         )

#         eps_s = eps_dual_branch_dc_unscaled(
#             duals_flow,
#             partitions,
#             zone_ids,
#             interzonal_branches,
#             eps_abs,
#             eps_rel,
#         )

#         _record_iter!(hist, r, s, eps_r, eps_s; rhos = Dict(:rho_flow => rho[]))

#         _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
#         _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

#         if debug && (iter % dbg_every == 0)
#             _dbg_print_iter!(
#                 iter,
#                 partitions,
#                 zone_ids,
#                 interzonal_branches,
#                 flow_consensus,
#                 duals_flow,
#                 rho[],
#                 prev_rho_before,
#                 r,
#                 s,
#                 eps_r,
#                 eps_s,
#                 results,
#             )
#         end

#         if stop_if_converged!(r, s, eps_r, eps_s)
#             _FAST_USE_LOGS && println("Stopping criterion met.")
#             break
#         end

#         maybe_adapt_rho!(adaptive, rho, r, s, config)
#     end

#     result = _finalize_admm_result(partitions, t0_ns, iter_done)
#     _print_admm_summary!(
#         "[ADMM-branch/DC]",
#         result;
#         rho_lines = ("rho_flow inicial: $(config.rho_flow_p) | rho_flow final: $(rho[])",),
#     )

#     return LirionRunResult(; result, model = nothing, system = sys, diagnostics = hist)
# end

