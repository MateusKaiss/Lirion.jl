
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

@inline function _zone_cut_edges(part::PartitionAC, cut_edges::Vector{Int})
    return [e for e in cut_edges if haskey(part.edge_to_local, e)]
end

@inline function _zone_cut_edges(part::PartitionDC, cut_edges::Vector{Int})
    return [e for e in cut_edges if haskey(part.edge_to_local, e)]
end

@inline function _zone_cut_edges(
    sub::Union{PartitionDC, PartitionAC, PartitionSOCP},
    cut_edges::Vector{Int},
)
    return [e for e in cut_edges if haskey(sub.edge_to_local, e)]
end


@inline _wrap_angle_diff(Δ::Float64) = atan(sin(Δ), cos(Δ))

@inline function _pretty_time(sec::Real)
    if sec < 60
        return "$(round(sec; digits = 3)) s"
    elseif sec < 3600
        return "$(round(sec/60; digits = 2)) min"
    else
        return "$(round(sec/3600; digits = 2)) h"
    end
end
