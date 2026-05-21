# ============================================================
# Result dictionary builder
# ============================================================

function build_result_dict(model::JuMP.Model, sys::SystemData)
    T = get(model.ext, :time_periods, sys.T)
    result = Dict{String, Any}()

    result["objective"] = objective_value(model)

    # -------------------------
    # Bus
    # -------------------------
    bus_ids = sys.bus.id

    # Voltage angles: ACR uses rectangular (vr, vi) and recovers θ = atan2(vi, vr).
    # Polar/DC/SOCWR expose `theta` directly.
    theta_matrix =
        haskey(model, :theta) ?
        [value(model[:theta][b, t]) for b in bus_ids, t in 1:T] :
        haskey(model, :v_real) && haskey(model, :v_imag) ?
        [
            atan(value(model[:v_imag][b, t]), value(model[:v_real][b, t]))
            for b in bus_ids, t in 1:T
        ] : nothing

    vm_matrix =
        haskey(model, :v_magnitude) ?
        [value(model[:v_magnitude][b, t]) for b in bus_ids, t in 1:T] :
        haskey(model, :v_real) && haskey(model, :v_imag) ?
        [
            hypot(value(model[:v_real][b, t]), value(model[:v_imag][b, t]))
            for b in bus_ids, t in 1:T
        ] :
        haskey(model, :phi) ?
        [1.0 + value(model[:phi][b, t]) for b in bus_ids, t in 1:T] :
        nothing

    vm_sq_matrix =
        haskey(model, :v_magnitude_sq) ?
        [value(model[:v_magnitude_sq][b, t]) for b in bus_ids, t in 1:T] :
        haskey(model, :v_real) && haskey(model, :v_imag) ?
        [
            value(model[:v_real][b, t])^2 + value(model[:v_imag][b, t])^2
            for b in bus_ids, t in 1:T
        ] : nothing

    result["bus"] = Dict(
        "theta" => theta_matrix,
        "vm" => vm_matrix,
        "vm_sq" => vm_sq_matrix,
    )

    # -------------------------
    # Gen (only in-service)
    # -------------------------
    if haskey(model, :p_generation) && sys.gen.n > 0
        gen_status = sys.gen.status
        gen_ids = [g for g in 1:length(gen_status) if gen_status[g] == 1]

        result["gen"] = Dict(
            "p" => [value(model[:p_generation][g, t]) for g in gen_ids, t in 1:T],
            "q" =>
                haskey(model, :q_generation) ?
                [value(model[:q_generation][g, t]) for g in gen_ids, t in 1:T] :
                nothing,
            "curtailment" =>
                haskey(model, :curtailment) ?
                [
                    try
                        value(model[:curtailment][g, t])
                    catch
                        0.0
                    end
                    for g in gen_ids, t in 1:T
                ] : nothing,
        )
    end

    # -------------------------
    # Branch (only in-service, local indexing)
    # -------------------------
    if haskey(model, :p_flow_branches) && sys.branch.n > 0
        br_status = sys.branch.status
        br_ids = [l for l in 1:length(br_status) if br_status[l] == 1]

        result["branch"] = Dict(
            "p_from" =>
                [value(model[:p_flow_branches][l, t]) for l in br_ids, t in 1:T],
            "p_to" => [
                value(model[:p_reverse_flow_branches][l, t]) for l in br_ids, t in 1:T
            ],
            "q_from" =>
                haskey(model, :q_flow_branches) ?
                [value(model[:q_flow_branches][l, t]) for l in br_ids, t in 1:T] :
                nothing,
            "q_to" =>
                haskey(model, :q_reverse_flow_branches) ?
                [value(model[:q_reverse_flow_branches][l, t]) for l in br_ids, t in 1:T] :
                nothing,
            "R" =>
                haskey(model, :R) ?
                [value(model[:R][l, t]) for l in br_ids, t in 1:T] :
                nothing,
            "T_var" =>
                haskey(model, :T_var) ?
                [value(model[:T_var][l, t]) for l in br_ids, t in 1:T] :
                nothing,
        )
    end

    # -------------------------
    # Storage
    # -------------------------
    if haskey(model, :storage) && sys.storage.n > 0
        nb = sys.storage.n
        result["storage"] = Dict(
            "storage" => [value(model[:storage][b, t]) for b in 1:nb, t in 1:T],
            "charge" => [value(model[:charge][b, t]) for b in 1:nb, t in 1:T],
            "discharge" => [value(model[:discharge][b, t]) for b in 1:nb, t in 1:T],
            "qsc" => [value(model[:qsc][b, t]) for b in 1:nb, t in 1:T],
            "ccms" =>
                haskey(model, :ccms) ?
                [value(model[:ccms][b, t]) for b in 1:nb, t in 1:T] : nothing,
            "p_injected" =>
                [value(model[:p_storage_injected][b, t]) for b in 1:nb, t in 1:T],
            "q_injected" =>
                [value(model[:q_storage_injected][b, t]) for b in 1:nb, t in 1:T],
        )
    end

    return result
end
