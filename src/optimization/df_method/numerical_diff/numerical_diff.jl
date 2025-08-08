#-------------------------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------------------------
struct NumericalDiff <: AbstractDfMethod
    name::String
    #------------------------------------------------------------------------------
    #
    #------------------------------------------------------------------------------
    function NumericalDiff()
        return new("Numerical_diff")
    end 
end
#-------------------------------------------------------------------------------------------------
# 自動微分による感度解析
#-------------------------------------------------------------------------------------------------
function compute_f(
    x::Vector{Float64},
    eval::EF, 
    physics::P, 
    opt_settings::OptimizationSettings,
    opt_filter::OF,
    ::NumericalDiff
    )::Float64 where {EF<:AbstractEvaluationFunction, P<:AbstractPhysics, OF<:AbstractOptFilter}

    # 変数の更新(.=を使うと感度がズレるので=を使うこと)
    eval.x = x

    #
    return compute_eval_func(eval, physics, opt_settings, opt_filter)

end
#-------------------------------------------------------------------------------------------------
# 随伴変数法による感度解析
#-------------------------------------------------------------------------------------------------
function compute_f_and_grad(
    x::Vector{Float64}, 
    eval::EF, 
    physics::P, 
    opt_settings::OptimizationSettings,
    opt_filter::OF, 
    df_method::NumericalDiff
    )::Tuple{Float64, Vector{Float64}} where {EF<:AbstractEvaluationFunction, P<:AbstractPhysics, OF<:AbstractOptFilter}

    # 変数の更新(.=を使うと感度がズレるので=を使うこと)
    eval.x = x

    # Compute function
    func = compute_eval_func(eval, physics, opt_settings, opt_filter)

    # Compute gradient by adjoint variable method
    df = central_difference_method(x, eval, physics, opt_settings, opt_filter, 1.0e-05, df_method)
    
    # Update for filtering and Finish
    return func, df
end
#----------------------------------------------------------------------------
# 中央差分法にて関数の微分を行う
#----------------------------------------------------------------------------
function central_difference_method(
    x::Vector{Float64}, 
    eval::EF, 
    physics::P, 
    opt_settings::OptimizationSettings,
    opt_filter::OF,
    delta::Float64,
    ::NumericalDiff
    )::Vector{Float64} where {EF<:AbstractEvaluationFunction, P<:AbstractPhysics, OF<:AbstractOptFilter}
    
    # Init
    df = similar(x)
    # Main
    for i in eachindex(df)
        # 保存しておく
        x_org = x[i]
        
        # 更新
        x[i] += delta
        # 評価関数構造体へ入力
        eval.x = x
        #
        f_new = compute_eval_func(eval, physics, opt_settings, opt_filter)
        # 戻す
        x[i] = x_org
        
        # 更新
        x[i] -= delta
        # 評価関数構造体へ入力
        eval.x = x
        #
        f_old = compute_eval_func(eval, physics, opt_settings, opt_filter)
        # 戻す
        x[i]  = x_org

        # 感度の計算
        df[i] = (f_new - f_old) / (2.0*delta)
    end

    #
    return df

end