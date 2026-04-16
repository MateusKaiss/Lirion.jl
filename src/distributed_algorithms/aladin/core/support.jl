abstract type AbstractALADINWorkspace <: AbstractDistributedWorkspace end

mutable struct LambdaState
    lambda::Vector{Float64}
    lambda_prev::Vector{Float64}
    nu::Vector{Float64}
end

LambdaState(n::Int) = LambdaState(zeros(n), zeros(n), zeros(n))

@inline stop_if_converged(ws::AbstractALADINWorkspace) =
    (ws.r_primal < ws.eps_primal) && (ws.r_progress < ws.eps_progress)

@inline _line_search_step!(::AbstractALADINWorkspace) = nothing

function _init_from_start_values(all_vars::Vector{JuMP.VariableRef})
    n = length(all_vars)
    v = Vector{Float64}(undef, n)
    @inbounds for k in 1:n
        sv = JuMP.start_value(all_vars[k])
        v[k] = sv === nothing ? 0.0 : Float64(sv)
    end
    return v
end

function _resolve_qp_optimizer(user_opt)
    user_opt === nothing || return user_opt
    return JuMP.optimizer_with_attributes(
        OSQP.Optimizer,
        "eps_abs" => 1e-8,
        "eps_rel" => 1e-8,
        "max_iter" => 200_000,
        "verbose" => false,
        "polish" => true,
        "polish_refine_iter" => 5,
        "linsys_solver" => "qdldl",
    )
end
