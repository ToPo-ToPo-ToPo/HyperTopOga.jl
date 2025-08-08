
#----------------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------------
function output_opt_result(x::Vector{Float64}, physics::AbstractPhysics, opt_settings::AbstractOptSettings, opt_filter::AbstractOptFilter)
    output_opt_result(x, physics, opt_settings, opt_filter, opt_settings.design_space_type)
end