

#--------------------------------------------------------------------------------------------------------
# 初期値に対する各目的関数の計算を行う
#--------------------------------------------------------------------------------------------------------
function initial_compute(opt::M, physics::P, opt_settings::OS, opt_filter::F, x::Vector{Float64}
    )::Vector{Float64} where {M<:AbstractMMA, P<:AbstractPhysics, OS<:AbstractOptSettings, F<:AbstractOptFilter}

    f0vals = Vector{Float64}(undef, length(opt.objective))
    
    for i in eachindex(opt.objective)
        f0vals[i] = compute_f(x, opt.objective[i], physics, opt_settings, opt_filter, opt.objective_df_methods[i])
    end
    
    return f0vals
end