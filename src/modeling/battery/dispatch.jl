function batteries!(
    model::JuMP.Model,
    sys::SystemData,
    model_type::Type{<:OPFModel},
)
    bat = sys.battery

    if bat.n == 0
        @info "No battery information in SystemData. Skipping storage modeling."
        return nothing
    end

    add_battery_variables!(model, bat)
    add_battery_constraints!(model, bat, model_type)

    return nothing
end

# =======================================================================
# VARIABLES
# =======================================================================

function add_battery_variables!(model::JuMP.Model, bat::BatteryData)
    add_storage_variables!(model, bat)
    add_charge_discharge_variables!(model, bat)
    add_battery_power_variables!(model, bat)
    return nothing
end

# =======================================================================
# CONSTRAINTS
# =======================================================================

function add_battery_constraints!(model::JuMP.Model, bat::BatteryData, ::Type{DC})
    return add_dc_battery_constraints!(model, bat)
end

function add_battery_constraints!(model::JuMP.Model, bat::BatteryData, ::Type{AC})
    return add_ac_battery_constraints!(model, bat)
end

function add_battery_constraints!(model::JuMP.Model, bat::BatteryData, ::Type{ACR})
    return add_acr_battery_constraints!(model, bat)
end

function add_battery_constraints!(model::JuMP.Model, bat::BatteryData, ::Type{IVR})
    return add_acr_battery_constraints!(model, bat)
end

function add_battery_constraints!(model::JuMP.Model, bat::BatteryData, ::Type{LPAC})
    return add_dc_battery_constraints!(model, bat)
end

function add_battery_constraints!(model::JuMP.Model, bat::BatteryData, ::Type{SOCWR})
    return add_socwr_battery_constraints!(model, bat)
end

function add_battery_constraints!(model::JuMP.Model, bat::BatteryData, ::Type{SOCBF})
    return add_socwr_battery_constraints!(model, bat)
end

function add_battery_constraints!(model::JuMP.Model, bat::BatteryData, ::Type{BFA})
    return add_dc_battery_constraints!(model, bat)
end
