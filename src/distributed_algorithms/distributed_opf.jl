# ============================================================
# Dispatch
# ============================================================

# ------------------------------------------------------------
# Coupling - variant type
# ------------------------------------------------------------

@inline function _variant_type(::Val{c}) where {c}
    c === :state && return State
    c === :power && return Power
    c === :hybrid && return Hybrid
    return error("Unknown coupling: $(c). Expected :state, :power, or :hybrid.")
end

# ------------------------------------------------------------
# Algorithm display name
# ------------------------------------------------------------

_algo_name(::ADMMConfig) = "ADMM"
_algo_name(::ALADINConfig) = "ALADIN"

# ============================================================
# Utilities
# ============================================================

function _with_scoped_blas_threads(f, n::Integer = 1)
    prev = LinearAlgebra.BLAS.get_num_threads()
    try
        LinearAlgebra.BLAS.set_num_threads(n)
        return f()
    finally
        LinearAlgebra.BLAS.set_num_threads(prev)
    end
end

# ============================================================
# Entry point
# ============================================================

function solve_distributed(inputs::DistributedInputs)
    cfg = inputs.model_type
    sys = inputs.system

    @info "Running distributed OPF ($(_algo_name(cfg)))"
    @info "  coupling   = $(cfg.coupling)"
    @info "  sectioning = $(cfg.sectioning)"
    @info "  base_model = $(cfg.base_model)"
    @info "  optimizer  = $(inputs.optimizer)"

    optimizer_factory = resolve_optimizer(cfg, inputs.optimizer)

    return _with_scoped_blas_threads(1) do
        return run_algorithm(sys, cfg; optimizer_factory = optimizer_factory)
    end
end
