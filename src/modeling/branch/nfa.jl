# =======================================================================
# NFA (Network Flow Approximation)
# =======================================================================
# Simplest OPF model: pure active-power network flow.
# No voltages, no angles, no reactive power, no Ohm's law.
# Each branch has a single free active-power flow variable bounded by
# thermal limits; the reverse flow is -p_flow (lossless).

function add_nfa_branch_variables!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    @variable(model, p_flow_branches[l in active, t in 1:T])

    @inbounds for l in active, t in 1:T
        ra = br.rate_a[l]
        set_lower_bound(p_flow_branches[l, t], -ra)
        set_upper_bound(p_flow_branches[l, t], ra)
    end

    return nothing
end

function add_nfa_branch_expressions!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    @expression(model, p_reverse_flow_branches[l in active, t in 1:T],
        -model[:p_flow_branches][l, t]
    )

    return nothing
end
