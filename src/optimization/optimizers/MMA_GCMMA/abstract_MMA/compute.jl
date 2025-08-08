
#--------------------------------------------------------------------------------------------------------
# 目的関数、制約関数の計算+それぞれの感度の計算を行う
#--------------------------------------------------------------------------------------------------------
function compute(opt::M, physics::P, opt_settings::OS, opt_filter::F, opt_step::Int64, x::Vector{Float64}
    )::Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Matrix{Float64}} where {M<:AbstractMMA, P<:AbstractPhysics, OS<:AbstractOptSettings, F<:AbstractOptFilter}
   
    # 目的関数と感度の計算
    f0val = 0.0
    f0vals = Vector{Float64}(undef, length(opt.objective))
    df0dx = zeros(Float64, opt.n)
    for i in eachindex(opt.objective)
        # Compute function and gradient
        f, df = compute_f_and_grad(x, opt.objective[i], physics, opt_settings, opt_filter, opt.objective_df_methods[i])
        # Consider sensitivity filter
        df = filtering_grad(opt_filter, df, x)
        # Set Weight
        weight = opt.objective_weight[i](opt_step) / opt.objective0[i]
        # Update
        f0vals[i] = weight * f
        f0val += f0vals[i]
        # Update sensitivity
        df0dx .+= weight .* df
    end

    # 制約関数の計算
    fval = Vector{Float64}(undef, opt.m)
    dfdx = Matrix{Float64}(undef, opt.m, opt.n)
    for i in eachindex(opt.constraint)
        # Compute function and gradient
        f, df = compute_f_and_grad(x, opt.constraint[i], physics, opt_settings, opt_filter, opt.constraint_df_methods[i])
        # Consider sensitivity filter
        df = filtering_grad(opt_filter, df, x)
        # Set limit
        limit = opt.constraint_lim[i](opt_step)
        # Update
        fval[i] = f / limit - 1.0
        # Update sensitivity
        dfdx[i, :] .= df ./ limit
    end
   
    return f0val, f0vals, df0dx, fval, dfdx
end