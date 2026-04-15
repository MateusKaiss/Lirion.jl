@inline _admm_post_solve!(::AbstractADMMWorkspace, ::Integer, ::Any) = nothing

local_solve!(ws::AbstractADMMWorkspace) =
    _solve_zones_parallel!(ws, _admm_post_solve!)
