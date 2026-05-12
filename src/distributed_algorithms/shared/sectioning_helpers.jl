# ============================================================
# Sectioning JuMP-side helpers (shared between ADMM and ALADIN)
# ============================================================

# ============================================================
# DC
# ============================================================

function _attach_p_inj_variables!(
    part::PartitionDC,
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_store::Dict{Int, JuMP.VariableRef},
)
    model = part.model::JuMP.Model
    pos = part.position
    zone_split_bus_ids = Int[]
    for (fbus, za, zb) in bus_pairs
        if za == pos || zb == pos
            push!(zone_split_bus_ids, fbus)
        end
    end
    isempty(zone_split_bus_ids) && return nothing

    bal_con = model[:active_power_balance]

    for b in zone_split_bus_ids
        v = @variable(model, base_name = "p_inj_$(b)")
        JuMP.set_start_value(v, 0.0)
        p_inj_store[b] = v
        JuMP.set_normalized_coefficient(bal_con[b, 1], v, 1.0)
    end
    return nothing
end

# ============================================================
# AC
# ============================================================

function _attach_pq_inj_variables_ac!(
    part::PartitionAC,
    sys_zone::SystemData,
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_store::Dict{Int, JuMP.VariableRef},
    q_inj_store::Dict{Int, JuMP.VariableRef},
)
    model = part.model::JuMP.Model
    pos = part.position
    pf = sys_zone.powerflow
    bus_idx_map = pf.bus_idx_map
    gs = pf.shunt_conductance
    bs = pf.shunt_susceptance

    zone_split_bus_ids = Int[]
    for (fbus, za, zb) in bus_pairs
        if za == pos || zb == pos
            push!(zone_split_bus_ids, fbus)
        end
    end
    isempty(zone_split_bus_ids) && return nothing

    bal_p = model[:active_power_balance]
    bal_q = model[:reactive_power_balance]
    v_mag = model[:v_magnitude]

    for b in zone_split_bus_ids
        vp = @variable(model, base_name = "p_inj_$(b)")
        JuMP.set_start_value(vp, 0.0)
        p_inj_store[b] = vp

        vq = @variable(model, base_name = "q_inj_$(b)")
        JuMP.set_start_value(vq, 0.0)
        q_inj_store[b] = vq

        JuMP.delete(model, bal_p[b, 1])
        @constraint(
            model,
            compute_active_power_balance(model, sys_zone, b, 1) -
            gs[bus_idx_map[b]] * v_mag[b, 1]^2 + vp == 0
        )

        JuMP.delete(model, bal_q[b, 1])
        @constraint(
            model,
            compute_reactive_power_balance(model, sys_zone, b, 1) +
            bs[bus_idx_map[b]] * v_mag[b, 1]^2 + vq == 0
        )
    end
    return nothing
end

# ============================================================
# SOCWR
# ============================================================

function _attach_pq_inj_variables_socwr!(
    part::PartitionSOCWR,
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_store::Dict{Int, JuMP.VariableRef},
    q_inj_store::Dict{Int, JuMP.VariableRef},
)
    model = part.model::JuMP.Model
    pos = part.position

    zone_split_bus_ids = Int[]
    for (fbus, za, zb) in bus_pairs
        if za == pos || zb == pos
            push!(zone_split_bus_ids, fbus)
        end
    end
    isempty(zone_split_bus_ids) && return nothing

    bal_p = model[:active_power_balance]
    bal_q = model[:reactive_power_balance]

    for b in zone_split_bus_ids
        vp = @variable(model, base_name = "p_inj_$(b)")
        JuMP.set_start_value(vp, 0.0)
        p_inj_store[b] = vp
        JuMP.set_normalized_coefficient(bal_p[b, 1], vp, 1.0)

        vq = @variable(model, base_name = "q_inj_$(b)")
        JuMP.set_start_value(vq, 0.0)
        q_inj_store[b] = vq
        JuMP.set_normalized_coefficient(bal_q[b, 1], vq, 1.0)
    end
    return nothing
end
