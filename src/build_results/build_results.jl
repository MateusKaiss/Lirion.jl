function build_result_dict(model::JuMP.Model, data::Dict)
    T = get(model.ext, :time_periods, 1)
    result = Dict{String, Any}()

    result["objective"] = objective_value(model)
    bus_ids = data["bus"]["bus_i"]
    result["bus"] = Dict(
        "theta" => [
            value(model[:theta][b, t]) for b in bus_ids, t in 1:T
        ],
        "vm" =>
            haskey(model, :v_magnitude) ?
            [value(model[:v_magnitude][b, t]) for b in bus_ids, t in 1:T] : nothing,
        "vm_sq" =>
            haskey(model, :v_magnitude_sq) ?
            [value(model[:v_magnitude_sq][b, t]) for b in bus_ids, t in 1:T] : nothing,
    )

    if haskey(model, :p_generation)
        gen_status = data["gen"]["gen_status"]
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
                ] : nothing)
    end

    if haskey(model, :p_flow_branches)
        branch_status = data["branch"]["br_status"]
        branch_ids = [g for g in 1:length(branch_status) if branch_status[g] == 1]
        result["branch"] = Dict(
            "p_from" => [value(model[:p_flow_branches][l, t]) for l in branch_ids, t in 1:T],
            "p_to"   => [value(model[:p_reverse_flow_branches][l, t]) for l in branch_ids, t in 1:T],
            "q_from" =>
                haskey(model, :q_flow_branches) ?
                [value(model[:q_flow_branches][l, t]) for l in branch_ids, t in 1:T] :
                nothing,
            "q_to"   =>
                haskey(model, :q_reverse_flow_branches) ?
                [value(model[:q_reverse_flow_branches][l, t]) for l in branch_ids, t in 1:T] :
                nothing,
            # --- extras só existem em SOCP ---
            "R" =>
                haskey(model, :R) ?
                [value(model[:R][l, t]) for l in branch_ids, t in 1:T] :
                nothing,
            "T_var" =>
                haskey(model, :T_var) ?
                [value(model[:T_var][l, t]) for l in branch_ids, t in 1:T] :
                nothing,
        )
    end

    if haskey(model, :storage)
        nb = length(data["storage"]["storage_bus"])
        result["storage"] = Dict(
            "storage" => [value(model[:storage][b, t]) for b in 1:nb, t in 1:T],
            "charge" => [value(model[:charge][b, t]) for b in 1:nb, t in 1:T],
            "discharge" => [value(model[:discharge][b, t]) for b in 1:nb, t in 1:T],
            "qsc" => [value(model[:qsc][b, t]) for b in 1:nb, t in 1:T],
            "ccms" =>
                haskey(model, :ccms) ?
                [value(model[:ccms][b, t]) for b in 1:nb, t in 1:T] : nothing,
            "p_injected" =>
                [value(model[:p_battery_injected][b, t]) for b in 1:nb, t in 1:T],
            "q_injected" =>
                [value(model[:q_battery_injected][b, t]) for b in 1:nb, t in 1:T],
        )
    end

    if haskey(model, :p_prosumer_injected)
        np = length(data["prosumer"]["bus"])
        result["prosumer"] = Dict(
            "p_injected" =>
                [value(model[:p_prosumer_injected][p, t]) for p in 1:np, t in 1:T],
            "q_injected" =>
                [value(model[:q_prosumer_injected][p, t]) for p in 1:np, t in 1:T],
            "pv_used" =>
                [value(model[:pv_generation_used][p, t]) for p in 1:np, t in 1:T],
            "curtailment" =>
                [value(model[:pv_curtailment][p, t]) for p in 1:np, t in 1:T],
        )
    end

    return result
end

function get_output_directory(case_path::AbstractString)
    case_dir = dirname(case_path)
    output_dir = joinpath(case_dir, "outputs")
    mkpath(output_dir)
    return output_dir
end

function save_subdict_to_csv(subdict::Dict, name::String, output_dir::String)
    df = DataFrame()
    for (k, v) in subdict
        if v !== nothing
            df[!, string(k)] = vec(v)
        else
            println("Field '$(k)' in section '$(name)' is `nothing` and was skipped.")
        end
    end
    CSV.write(joinpath(output_dir, "$name.csv"), df)
    return nothing
end

function save_result_to_csv(result_dict::Dict, data::Dict, case_path::AbstractString)
    output_dir = get_output_directory(case_path)

    for (section, subdict) in result_dict
        if subdict isa Dict
            save_subdict_to_csv(subdict, section, output_dir)
        elseif section == "objective"
            CSV.write(
                joinpath(output_dir, "objective.csv"),
                DataFrame(; objective = [subdict]),
            )
        end
    end

    return println("Resultados salvos em: $output_dir")
end
