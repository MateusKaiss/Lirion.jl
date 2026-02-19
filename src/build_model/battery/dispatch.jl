
# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================

function batteries!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
)
    bat = sys.battery
    if bat.n == 0
        @info "No battery information in SystemData (DC). Skipping storage modeling."
        return nothing
    end

    create_battery_variables!(model, bat)
    create_dc_battery_constraints!(model, bat)
    return nothing
end

function batteries!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{AC},
)
    bat = sys.battery
    if bat.n == 0
        @info "No battery information in SystemData (AC). Skipping storage modeling."
        return nothing
    end

    create_battery_variables!(model, bat)
    create_ac_battery_constraints!(model, bat)
    return nothing
end

function batteries!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{SOCP},
)
    bat = sys.battery
    if bat.n == 0
        @info "No battery information in SystemData (SOCP). Skipping storage modeling."
        return nothing
    end

    create_battery_variables!(model, bat)
    create_socp_battery_constraints!(model, bat)
    return nothing
end
