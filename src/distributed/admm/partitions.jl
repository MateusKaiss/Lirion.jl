
function _set_model_silent!(m::JuMP.Model)
    _safe(k, v) = (
        try
            set_optimizer_attribute(m, k, v)
        catch
        end
    )

    try
        set_silent(m)
    catch
    end
    try
        set_optimizer_attribute(m, JuMP.MOI.Silent(), true)
    catch
    end

    solver = ""
    try
        backend = JuMP.backend(m)
        solver = lowercase(JuMP.MOI.get(backend, JuMP.MOI.SolverName()))
    catch
    end

    if occursin("ipopt", solver)
        _safe("print_level", 0)
        _safe("sb", "yes")
    elseif occursin("gurobi", solver)
        _safe("OutputFlag", 0)
    elseif occursin("cplex", solver)
        _safe("CPXPARAM_ScreenOutput", 0)
    elseif occursin("highs", solver)
        _safe("log_to_console", false)
        _safe("output_flag", false)
    elseif occursin("scs", solver)
        _safe("verbose", false)
    elseif occursin("clarabel", solver)
        _safe("verbose", false)
    elseif occursin("glpk", solver)
        _safe("msg_lev", 0)
    else
        _safe("verbose", false)
        _safe("output_flag", 0)
        _safe("print_level", 0)
    end
    return nothing
end

function _zone_cut_edges(subdata::Dict{String, Any}, cut_edges::Vector{Int})
    local_keep = Set(get(subdata, "branch_local_indices", Int[]))
    return [ℓ for ℓ in cut_edges if ℓ in local_keep]
end

# ---------- Helpers ----------
@inline function _edge_local(sub::Dict{String, Any}, e::Int)
    edgemap = get(sub, "edge_to_local", nothing)
    return edgemap === nothing ? findfirst(==(e), sub["branch_local_indices"]) :
           get(edgemap, e, nothing)
end

@inline function _bus_set(sub::Dict{String, Any})
    return get(sub, "bus_set", Set(Int[b for b in sub["bus"]["bus_i"]]))
end
