function resolve_optimizer_dc(config)
    opt = try
        getfield(config, :optimizer)
    catch
        nothing
    end

    if opt !== nothing
        if opt == COSMO.Optimizer
            return optimizer_with_attributes(
                COSMO.Optimizer,
                "eps_abs" => 1e-6,
                "eps_rel" => 1e-6,
                "max_iter" => 10_000,
                "verbose" => false,
                "decompose" => true,
            )
        else
            return opt
        end
    end

    eps_abs = try
        getfield(config, Symbol("ε_primal"))
    catch
        1e-6
    end
    eps_rel = try
        getfield(config, Symbol("ε_dual"))
    catch
        1e-6
    end

    return optimizer_with_attributes(
        OSQP.Optimizer,
        "eps_abs" => eps_abs,
        "eps_rel" => eps_rel,
        "max_iter" => 50_000,
        "verbose" => false,
        "adaptive_rho" => true,
        "linsys_solver" => "qdldl",
    )
end

####################
# PEQUENOS HELPERS #
####################


@inline function _val_or_start(m::JuMP.Model, x::JuMP.VariableRef)
    if JuMP.has_values(m)
        return JuMP.value(x)
    else
        sv = JuMP.start_value(x)
        return sv === nothing ? 0.0 : sv
    end
end

@inline function _theta_val(
    m::JuMP.Model,
    θ::JuMP.Containers.DenseAxisArray,
    b::Int,
    t::Int,
)
    return _val_or_start(m, θ[b, t])
end

@inline function snapshot_theta!(part::PartitionDC)
    # θ "completo" (todas as barras da zona)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.theta_buf
    bus_ids = part.bus.id
    @inbounds for i in eachindex(bus_ids)
        buf[i] = _theta_val(m, θ, bus_ids[i], 1)
    end
    return buf
end

@inline function snapshot_theta_shared!(part::PartitionDC)
    # θ só nas barras compartilhadas da zona
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.theta_shared_buf
    pos = part.shared_pos_local
    bus_ids = part.bus.id
    @inbounds for i in eachindex(pos)
        bidx = pos[i]
        buf[i] = _theta_val(m, θ, bus_ids[bidx], 1)
    end
    return buf
end

#############################
# FLUXO DC EM UM RAMO       #
#############################

@inline function _dc_flow_val_model(
    m::JuMP.Model,
    part::PartitionDC,
    l_local::Int;
    t::Int = 1,
)
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    f = part.branch.from[l_local]
    tb = part.branch.to[l_local]
    x = part.branch.x[l_local]
    @assert abs(x) > 1e-8 "x≈0 no ramo local $l_local (DC)."
    return (_theta_val(m, θ, f, t) - _theta_val(m, θ, tb, t)) / x
end

#############################
# EXPRESSÕES QUADRÁTICAS    #
#############################

function _build_quad_bus_expr!(part::PartitionDC, ρ::Float64)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    q = zero(JuMP.QuadExpr)

    # termo (ρ/2)*‖θ_shared‖² pra cada barra compartilhada
    for b in part.shared_buses_local
        q += 0.5 * ρ * θ[b, 1]^2
    end
    part.quad_bus_expr = q
    part.rho_bus_cached = ρ
    return q
end

function _build_quad_branch_expr!(
    part::PartitionDC,
    ρ::Float64,
    zone_edges::Vector{Int},
)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    q = zero(JuMP.QuadExpr)

    # termo (ρ/2)*‖p_flow(e)‖² para cada ramo de corte da zona
    for e in zone_edges
        l_local = get(part.edge_to_local, e, nothing)
        l_local === nothing && continue
        f = part.branch.from[l_local]
        t = part.branch.to[l_local]
        x = part.branch.x[l_local]
        @assert abs(x) > 1e-10 "x≈0 no ramo global $e (DC)."
        p = (θ[f, 1] - θ[t, 1]) / x
        q += 0.5 * ρ * p^2
    end
    part.quad_branch_expr = q
    part.rho_branch_cached = ρ
    return q
end

#############################
# EXPRESSÕES LINEARES       #
#############################

function _build_linear_bus_expr(
    part::PartitionDC,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    θ̄_vec::Vector{Float64},
    λ_vec_zone::Vector{Float64},
    ρ::Float64,
)
    m = part.model::JuMP.Model 
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr) #cria expressão linear vazia
    for b in part.shared_buses_local
        j = b2pos[b]
        lin += (λ_vec_zone[j] - ρ * θ̄_vec[j]) * θ[b, 1]
    end
    return lin
end

function _build_linear_branch_expr(
    part::PartitionDC,
    zone_edges::Vector{Int},
    z_flow::Dict{Int, Float64},
    λe_zone::Dict{Int, Float64},
    ρ::Float64,
)
    # soma_e ( λ_e - ρ z_e ) p_flow_e(θ)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)
    for e in zone_edges
        l_local = get(part.edge_to_local, e, nothing)
        l_local === nothing && continue
        f = part.branch.from[l_local]
        t = part.branch.to[l_local]
        x = part.branch.x[l_local]
        @assert abs(x) > 1e-10 "x≈0 no ramo global $e (DC)."
        coef = get(λe_zone, e, 0.0) - ρ * z_flow[e]
        lin += coef * (θ[f, 1] - θ[t, 1]) / x
    end
    return lin
end



function compute_theta_consensus_vec!(
    θ̄_vec::Vector{Float64},
    partitions::Dict{Int, PartitionDC},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    fill!(θ̄_vec, 0.0)
    counts = zeros(Int, length(shared_buses_sorted))

    # média das θ nas barras compartilhadas
    for (_, part) in partitions
        θs = snapshot_theta_shared!(part)
        lbs = part.shared_buses_local
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            θ̄_vec[j] += θs[i]
            counts[j] += 1
        end
    end

    @inbounds for j in eachindex(θ̄_vec)
        counts[j] != 0 && (θ̄_vec[j] /= counts[j])
    end
end

function compute_z_flow_consensus!(
    z_flow::Dict{Int, Float64},
    partitions::Dict{Int, PartitionDC},
    cut_edges::Vector{Int},
)
    sums = Dict{Int, Float64}()
    cnts = Dict{Int, Int}()

    # média dos fluxos p_e em cada ramo de corte
    for (_, part) in partitions
        m = part.model::JuMP.Model
        @inbounds for e in cut_edges
            l_local = get(part.edge_to_local, e, nothing)
            if l_local !== nothing
                p = _dc_flow_val_model(m, part, l_local)
                sums[e] = get(sums, e, 0.0) + p
                cnts[e] = get(cnts, e, 0) + 1
            end
        end
    end

    @inbounds for e in cut_edges
        n = get(cnts, e, 0)
        n > 0 && (z_flow[e] = sums[e] / n)
    end
end

#############################
# DUAIS EM θ                #
#############################

@inline function update_duals_theta_vec!(
    dual_vec_zone::Vector{Float64},
    part::PartitionDC,
    θ̄_vec::Vector{Float64},
    ρ::Float64,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    θs = snapshot_theta_shared!(part)
    lbs = part.shared_buses_local

    @inbounds for i in eachindex(lbs)
        j = b2pos[lbs[i]]
        dual_vec_zone[j] += ρ * (θs[i] - θ̄_vec[j])
    end
    return nothing
end

@inline function gather_zone_lambda!(
    λ_vec_zone::Vector{Float64},
    duals_dense_zone::Vector{Float64},
)
    λ_vec_zone .= duals_dense_zone
    return λ_vec_zone
end

#############################
# RESÍDUOS (BUS)            #
#############################

function primal_residual_bus_dc(
    partitions::Dict{Int, PartitionDC},
    θ̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    acc = 0.0
    for (_, part) in partitions
        θbuf = snapshot_theta!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_local
            j = b2pos[b]
            acc += (θbuf[pos[b]] - θ̄_vec[j])^2
        end
    end
    return sqrt(acc)
end

dual_residual_bus_dc(
    θ̄_vec::Vector{Float64},
    θ̄_prev_vec::Vector{Float64},
) = sqrt(sum((θ̄_vec[j] - θ̄_prev_vec[j])^2 for j in eachindex(θ̄_vec)))

function eps_primal_bus_dc(
    partitions::Dict{Int, PartitionDC},
    θ̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    N = max(1, length(shared_buses_sorted))
    norm_x_sq = 0.0
    for (_, part) in partitions
        θbuf = snapshot_theta!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_local
            norm_x_sq += (θbuf[pos[b]])^2
        end
    end
    norm_z_sq = sum((θ̄_vec[j])^2 for j in eachindex(θ̄_vec))

    return sqrt(N) * ϵ_abs + ϵ_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

function eps_dual_bus_dc(
    θ̄_vec::Vector{Float64},
    θ̄_prev_vec::Vector{Float64},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    N = max(1, length(θ̄_vec))
    return sqrt(N) * ϵ_abs +
           ϵ_rel * sqrt(sum((θ̄_vec[j] - θ̄_prev_vec[j])^2 for j in eachindex(θ̄_vec)))
end

#############################
# RESÍDUOS (BRANCH)         #
#############################


function primal_residual_branch_dc(
    partitions::Dict{Int, PartitionDC},
    z_flow::Dict{Int, Float64},
    cut_edges::Vector{Int},
)
    acc = 0.0
    for (_, part) in partitions
        m = part.model::JuMP.Model
        @inbounds for e in cut_edges
            l_local = get(part.edge_to_local, e, nothing)
            if l_local !== nothing
                acc += (_dc_flow_val_model(m, part, l_local) - z_flow[e])^2
            end
        end
    end
    return sqrt(acc)
end

dual_residual_branch_dc(
    z_flow::Dict{Int, Float64},
    z_flow_prev::Dict{Int, Float64},
) = sqrt(sum((z_flow[e] - get(z_flow_prev, e, 0.0))^2 for e in keys(z_flow)))

function eps_primal_branch_dc(
    partitions::Dict{Int, PartitionDC},
    z_flow::Dict{Int, Float64},
    cut_edges::Vector{Int},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    dim_pri = sum(length(_zone_cut_edges(part, cut_edges)) for (_, part) in partitions)
    return sqrt(max(1, dim_pri)) * ϵ_abs +
           ϵ_rel * max(sqrt(sum(z_flow[e]^2 for e in cut_edges)), 1e-12)
end

function eps_dual_branch_dc(
    z_flow::Dict{Int, Float64},
    z_flow_prev::Dict{Int, Float64},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    dim_du = max(1, length(z_flow))
    return sqrt(dim_du) * ϵ_abs +
           ϵ_rel * sqrt(sum((z_flow[e] - get(z_flow_prev, e, 0.0))^2 for e in keys(z_flow)))
end

##########################
# ρ ADAPTATIVO           #
##########################

function adaptive_rho_dc!(
    ρ::Base.RefValue{Float64},
    r::Float64,
    s::Float64;
    μ_incr = 3.0,
    μ_decr = 3.0,
    τ_incr = 1.05,
    τ_decr = 0.98,
)
    if r > μ_incr * s
        ρ[] *= τ_incr
    elseif s > μ_decr * r
        ρ[] *= τ_decr
    end
end

##########################
# DUAIS (FLOW)           #
##########################

function update_duals_flow!(
    duals::Dict{Int, Dict{Int, Float64}},
    partitions::Dict{Int, PartitionDC},
    z_flow::Dict{Int, Float64},
    cut_edges::Vector{Int},
    ρ::Float64,
)
    for (z, part) in partitions
        m = part.model::JuMP.Model
        @inbounds for e in cut_edges
            l_local = get(part.edge_to_local, e, nothing)
            if l_local !== nothing
                p_val = _dc_flow_val_model(m, part, l_local)
                duals[z][e] = get(duals[z], e, 0.0) + ρ * (p_val - z_flow[e])
            end
        end
    end
    return nothing
end


##########################
# FLAGS DE CONFIG        #
##########################

@inline _adaptive_on(config) =
    hasfield(typeof(config), :adaptive_rho) ? getfield(config, :adaptive_rho) : true

# =======================================================================
# Runners
# =======================================================================

##########################
# ADMM: BUS / DC         #
##########################

function run_admm(config::DistributedOPFModel, sys::ADMMSystem)
    t0_ns = time_ns()
    iter_done = 0
    @assert config.base_model == DC "run_admm (bus/DC) requer DC."

    partitions = build_dc_partitions(sys)
    shared_buses = sys.shared_buses

    for (_, part) in partitions
        prepare_hotpaths!(part)
    end
    shared_buses_sorted, b2pos = prepare_dense_maps(shared_buses)

    ρ = Ref(config.rho_theta)
    ϵ_abs = config.ε_primal
    ϵ_rel = config.ε_dual
    adaptive = _adaptive_on(config)

    duals_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    θ̄_vec = zeros(Float64, length(shared_buses_sorted))
    θ̄_prev_vec = similar(θ̄_vec)

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rho_vec = Float64[]
    θ_avg_history = Dict{Int, Vector{Float64}}()
    θ_local_history = Dict{Int, Dict{Int, Vector{Float64}}}()
    optimizer_factory = resolve_optimizer_dc(config)

    for (zone_id, part) in partitions
        sys_zone = sys.parts[zone_id]

        m = build_model(
            sys_zone,
            DC;
            optimizer = optimizer_factory,
        )

        part.model = m
        part.obj_base = m[:generation_cost]
        part.theta_var = m[:theta]
        part.has_solution = false

        _set_model_silent!(m)
        for I in keys(part.theta_var)
            set_start_value(part.theta_var[I], 0.0)
        end

        _build_quad_bus_expr!(part, ρ[])
        @objective(m, Min, part.obj_base + part.quad_bus_expr)
    end

    zone_ids = collect(keys(partitions))
    λ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )

    prev_rho = ρ[]

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\nIteração $k (bus/DC)")

        if ρ[] != prev_rho
            for (_, part) in partitions
                _build_quad_bus_expr!(part, ρ[])
            end
            prev_rho = ρ[]
        end

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)

        @threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            λz = gather_zone_lambda!(λ_zone_buf[z], duals_dense[z])

            lin =
                _build_linear_bus_expr(part, shared_buses_sorted, b2pos, θ̄_vec, λz, ρ[])

            m = part.model::JuMP.Model
            @objective(m, Min, part.obj_base + part.quad_bus_expr + lin)
            optimize!(m)

            gen_cost_expr = value(m[:generation_cost])
            real_cost = compute_real_generation_cost(m, part.objective_params)

            θs = snapshot_theta_shared!(part)
            θsnap = Dict{Int, Float64}()
            lbs = part.shared_buses_local
            @inbounds for i2 in eachindex(lbs)
                θsnap[lbs[i2]] = θs[i2]
            end
            results[i] =
                (gen_cost_expr = gen_cost_expr, real_cost = real_cost, θsnap = θsnap)
        end

        for (i, z) in enumerate(zone_ids)
            part = partitions[z]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost
            @inbounds for (b, θb) in out.θsnap

                # nível 1: zona
            # if !haskey(θ_local_history, z)
            #     θ_local_history[z] = Dict{Int, Vector{Float64}}()
            # end

            # # nível 2: barra
            # if !haskey(θ_local_history[z], b)
            #     θ_local_history[z][b] = Float64[]
            # end

            # # guarda o valor
            # push!(θ_local_history[z][b], θb)

                push!(
                    get!(
                        get!(θ_local_history, z, Dict{Int, Vector{Float64}}()),
                        b,
                        Float64[],
                    ),
                    θb,
                )
            end
        end

        θ̄_prev_vec .= θ̄_vec
        compute_theta_consensus_vec!(θ̄_vec, partitions, shared_buses_sorted, b2pos)
        @inbounds for b in shared_buses_sorted
            j = b2pos[b]
            push!(get!(θ_avg_history, b, Float64[]), θ̄_vec[j])
        end

        for z in zone_ids
            update_duals_theta_vec!(
                duals_dense[z],
                partitions[z],
                θ̄_vec,
                ρ[],
                shared_buses_sorted,
                b2pos,
            )
        end

        res_primal = primal_residual_bus_dc(partitions, θ̄_vec, shared_buses_sorted, b2pos)
        res_dual = dual_residual_bus_dc(θ̄_vec, θ̄_prev_vec)
        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        ϵ_pri =
            eps_primal_bus_dc(partitions, θ̄_vec, shared_buses_sorted, b2pos, ϵ_abs, ϵ_rel)
        ϵ_du = eps_dual_bus_dc(θ̄_vec, θ̄_prev_vec, ϵ_abs, ϵ_rel)
        push!(eps_primal_vec, ϵ_pri)
        push!(eps_dual_vec, ϵ_du)

        _FAST_USE_LOGS && println("‖r‖ = $res_primal | ε_pri = $ϵ_pri")
        _FAST_USE_LOGS && println("‖s‖ = $res_dual   | ε_du  = $ϵ_du")

        push!(rho_vec, ρ[])

        iter_done = k
        if res_primal < ϵ_pri && res_dual < ϵ_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_dc!(
                ρ,
                res_primal,
                res_dual;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
        end
    end

    save_dc_bus_plots(
        res_primal_vec,
        res_dual_vec,
        eps_primal_vec,
        eps_dual_vec,
        rho_vec,
        θ_avg_history;
        stem = "admm",
    )

    total_real_cost = sum(part.real_generation_cost for part in values(partitions))
    total_expr_cost = sum(part.generation_cost_model_expr for part in values(partitions))

    println("\nResumo (bus/DC):")
    println("→ Total (value(:generation_cost)) = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ inicial: $(config.rho_theta) | ρθ final: $(ρ[])")
    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (bus/DC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end

#################################
# ADMM: BRANCH/DC               #
#################################

function run_admm_branch(
    config::DistributedOPFModel,
    sys::ADMMSystem;
    use_legacy_rho::Bool = false,
    μ_incr::Float64 = 2.0,
    μ_decr::Float64 = 2.0,
    τ_incr::Float64 = 1.05,
    τ_decr::Float64 = 0.98,
)
    t0_ns = time_ns()
    iter_done = 0

    @assert config.base_model == DC "ADMM-branch/DC requer base_model == DC."

    partitions = build_dc_partitions(sys)
    cut_edges = sys.cut_edges

    for (_, part) in partitions
        prepare_hotpaths!(part)
        part.zone_edges = _zone_cut_edges(part, cut_edges)
    end

    ρ = Ref(config.rho_flow)
    ϵ_abs = config.ε_primal
    ϵ_rel = config.ε_dual
    adaptive = _adaptive_on(config)

    μ_incr = config.mu_incr
    μ_decr = config.mu_decr
    τ_incr = config.tau_incr
    τ_decr = config.tau_decr

    z_flow = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    z_flow_prev = Dict{Int, Float64}()
    duals =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rho_vec = Float64[]

    optimizer_factory = resolve_optimizer_dc(config)

    # build modelo base de cada zona (sem warm-start)
    for (zone_id, part) in partitions
        sys_zone = sys.parts[zone_id]

        m = build_model(
            sys_zone,
            DC;
            optimizer = optimizer_factory,
        )
        part.model = m
        part.obj_base = m[:generation_cost]
        part.theta_var = m[:theta]
        part.has_solution = false

        _set_model_silent!(m)
        for I in keys(part.theta_var)
            set_start_value(part.theta_var[I], 0.0)
        end

        _build_quad_branch_expr!(part, ρ[], part.zone_edges)
        @objective(m, Min, part.obj_base + part.quad_branch_expr)
    end

    zone_ids = collect(keys(partitions))
    prev_rho = ρ[]

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-branch/DC] Iteração $k")

        if ρ[] != prev_rho
            for (_, part) in partitions
                _build_quad_branch_expr!(part, ρ[], part.zone_edges)
            end
            prev_rho = ρ[]
        end

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)

        @threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]
            lin = _build_linear_branch_expr(part, part.zone_edges, z_flow, duals[z], ρ[])
            m = part.model::JuMP.Model
            @objective(m, Min, part.obj_base + part.quad_branch_expr + lin)
            optimize!(m)
            gen_cost_expr = value(m[:generation_cost])
            real_cost = compute_real_generation_cost(m, part.objective_params)
            results[i] = (gen_cost_expr = gen_cost_expr, real_cost = real_cost)
        end

        for (i, z) in enumerate(zone_ids)
            part = partitions[z]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost
        end

        z_flow_prev = copy(z_flow)
        compute_z_flow_consensus!(z_flow, partitions, cut_edges)
        update_duals_flow!(duals, partitions, z_flow, cut_edges, ρ[])

        res_primal = primal_residual_branch_dc(partitions, z_flow, cut_edges)
        res_dual = dual_residual_branch_dc(z_flow, z_flow_prev)
        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        ϵ_pri = eps_primal_branch_dc(partitions, z_flow, cut_edges, ϵ_abs, ϵ_rel)
        ϵ_du = eps_dual_branch_dc(z_flow, z_flow_prev, ϵ_abs, ϵ_rel)
        push!(eps_primal_vec, ϵ_pri)
        push!(eps_dual_vec, ϵ_du)

        _FAST_USE_LOGS && println("‖r‖ = $res_primal | ε_pri = $ϵ_pri")
        _FAST_USE_LOGS && println("‖s‖ = $res_dual   | ε_dual = $ϵ_du")

        push!(rho_vec, ρ[])

        iter_done = k
        if res_primal < ϵ_pri && res_dual < ϵ_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            if use_legacy_rho
                if res_primal > μ_incr * res_dual
                    ρ[] *= τ_incr
                elseif res_dual > μ_decr * res_primal
                    ρ[] *= τ_decr
                end
            else
                adaptive_rho_dc!(
                    ρ,
                    res_primal,
                    res_dual;
                    μ_incr = μ_incr,
                    μ_decr = μ_decr,
                    τ_incr = τ_incr,
                    τ_decr = τ_decr,
                )
            end
        end
    end

    save_dc_branch_plots(
        res_primal_vec,
        res_dual_vec,
        eps_primal_vec,
        eps_dual_vec,
        rho_vec;
        stem = "admm_branch",
    )

    total_real_cost = sum(part.real_generation_cost for part in values(partitions))
    total_expr_cost = sum(part.generation_cost_model_expr for part in values(partitions))

    println("\n[ADMM-branch/DC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρe_final = $(ρ[])")
    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (branch/DC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end

#################################
# ADMM: HÍBRIDO (BUS + BRANCH)  #
#################################

function run_admm_hybrid(config::DistributedOPFModel, sys::ADMMSystem)
    t0_ns = time_ns()
    iter_done = 0

    @assert config.base_model == DC "Hybrid ADMM/DC requer DC."

    partitions = build_dc_partitions(sys)
    shared_buses = sys.shared_buses
    cut_edges = sys.cut_edges

    for (_, part) in partitions
        prepare_hotpaths!(part)
        part.zone_edges = _zone_cut_edges(part, cut_edges)
    end

    shared_buses_sorted, b2pos = prepare_dense_maps(shared_buses)

    ρθ = Ref(config.rho_theta)
    ρe = Ref(config.rho_flow)
    ϵ_abs = config.ε_primal
    ϵ_rel = config.ε_dual
    adaptive = _adaptive_on(config)

    θ̄_vec = zeros(Float64, length(shared_buses_sorted))
    θ̄_prev_vec = similar(θ̄_vec)

    z_flow = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    z_flow_prev = Dict{Int, Float64}()

    duals_theta =
        Dict{Int, Vector{Float64}}(
            z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
        )
    duals_flow =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rhoθ_vec = Float64[]
    rhoe_vec = Float64[]
    θ_avg_history = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer_dc(config)

    # build base model de cada zona + termos quadráticos em θ_shared e p_flow de corte
    for (zone_id, part) in partitions
        sys_zone = sys.parts[zone_id]
        m = build_model(
            sys_zone,
            DC;
            optimizer = optimizer_factory,
        )
        part.model = m
        part.obj_base = m[:generation_cost]
        part.theta_var = m[:theta]
        part.has_solution = false
        _set_model_silent!(m)
        for I in keys(part.theta_var)
            set_start_value(part.theta_var[I], 0.0)
        end
        _build_quad_bus_expr!(part, ρθ[])
        _build_quad_branch_expr!(part, ρe[], part.zone_edges)
        @objective(m, Min, part.obj_base + part.quad_bus_expr + part.quad_branch_expr)
    end

    zone_ids = collect(keys(partitions))
    λθ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )

    prev_ρθ = ρθ[]
    prev_ρe = ρe[]

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-hybrid/DC] Iteração $k")

        if ρθ[] != prev_ρθ
            for (_, part) in partitions
                _build_quad_bus_expr!(part, ρθ[])
            end
            prev_ρθ = ρθ[]
        end
        if ρe[] != prev_ρe
            for (_, part) in partitions
                _build_quad_branch_expr!(part, ρe[], part.zone_edges)
            end
            prev_ρe = ρe[]
        end

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)

        @threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            λθ = gather_zone_lambda!(λθ_zone_buf[z], duals_theta[z])

            lin_bus = _build_linear_bus_expr(
                part,
                shared_buses_sorted,
                b2pos,
                θ̄_vec,
                λθ,
                ρθ[],
            )
            lin_branch = _build_linear_branch_expr(
                part,
                part.zone_edges,
                z_flow,
                duals_flow[z],
                ρe[],
            )

            m = part.model::JuMP.Model
            @objective(
                m,
                Min,
                part.obj_base + part.quad_bus_expr + part.quad_branch_expr +
                (lin_bus + lin_branch),
            )
            optimize!(m)

            gen_cost_expr = value(m[:generation_cost])
            real_cost = compute_real_generation_cost(m, part.objective_params)
            results[i] = (gen_cost_expr = gen_cost_expr, real_cost = real_cost)
        end

        for (i, z) in enumerate(zone_ids)
            part = partitions[z]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost
        end

        θ̄_prev_vec .= θ̄_vec
        z_flow_prev = copy(z_flow)

        compute_theta_consensus_vec!(θ̄_vec, partitions, shared_buses_sorted, b2pos)
        for b in shared_buses_sorted
            j = b2pos[b]
            push!(get!(θ_avg_history, b, Float64[]), θ̄_vec[j])
        end
        compute_z_flow_consensus!(z_flow, partitions, cut_edges)

        for z in zone_ids
            update_duals_theta_vec!(
                duals_theta[z],
                partitions[z],
                θ̄_vec,
                ρθ[],
                shared_buses_sorted,
                b2pos,
            )
        end
        update_duals_flow!(
            duals_flow,
            partitions,
            z_flow,
            cut_edges,
            ρe[],
        )

        res_pri_bus = primal_residual_bus_dc(partitions, θ̄_vec, shared_buses_sorted, b2pos)
        res_pri_br = primal_residual_branch_dc(partitions, z_flow, cut_edges)
        res_primal = hypot(res_pri_bus, res_pri_br)

        res_du_bus = dual_residual_bus_dc(θ̄_vec, θ̄_prev_vec)
        res_du_br = dual_residual_branch_dc(z_flow, z_flow_prev)
        res_dual = hypot(res_du_bus, res_du_br)

        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        ϵb_pri =
            eps_primal_bus_dc(partitions, θ̄_vec, shared_buses_sorted, b2pos, ϵ_abs, ϵ_rel)
        ϵe_pri = eps_primal_branch_dc(partitions, z_flow, cut_edges, ϵ_abs, ϵ_rel)
        push!(eps_primal_vec, hypot(ϵb_pri, ϵe_pri))

        ϵb_du = eps_dual_bus_dc(θ̄_vec, θ̄_prev_vec, ϵ_abs, ϵ_rel)
        ϵe_du = eps_dual_branch_dc(z_flow, z_flow_prev, ϵ_abs, ϵ_rel)
        push!(eps_dual_vec, hypot(ϵb_du, ϵe_du))

        _FAST_USE_LOGS &&
            println("‖r‖ = $(last(res_primal_vec)) | ε_pri = $(last(eps_primal_vec))")
        _FAST_USE_LOGS &&
            println("‖s‖ = $(last(res_dual_vec))   | ε_du  = $(last(eps_dual_vec))")

        push!(rhoθ_vec, ρθ[])
        push!(rhoe_vec, ρe[])

        iter_done = k
        if res_primal < last(eps_primal_vec) && res_dual < last(eps_dual_vec)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_dc!(
                ρθ,
                res_pri_bus,
                res_du_bus;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
            adaptive_rho_dc!(
                ρe,
                res_pri_br,
                res_du_br;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
        end
    end

    # Plots / resumo
    plot_admm_convergence(
        res_primal_vec,
        res_dual_vec,
        eps_primal_vec,
        eps_dual_vec,
        "plot_convergencia_admm_hybrid.html",
    )

    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoe_vec),
                    y = rhoe_vec,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρe",
                );
                p.layout.title = "Evolução de ρe";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρe";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rho_flow_hybrid.html",
    )

    plot_theta_avg_heatmap(θ_avg_history, "plot_theta_avg_heatmap_hybrid.html")

    total_real_cost = sum(part.real_generation_cost for part in values(partitions))
    total_expr_cost = sum(part.generation_cost_model_expr for part in values(partitions))

    println("\n[ADMM-hybrid/DC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ_final = $(ρθ[]) | ρe_final = $(ρe[])")
    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (hybrid/DC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end