
#-------------------------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------------------------
struct AdjointDiff <: AbstractDfMethod 
    name::String
    #---------------------------------------------------
    #
    #---------------------------------------------------
    function AdjointDiff()
        return new("Adjoint_diff")
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
    ::AdjointDiff
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
    ::AdjointDiff
    )::Tuple{Float64, Vector{Float64}} where {
        EF<:AbstractEvaluationFunction, 
        P<:AbstractPhysics, 
        OF<:AbstractOptFilter
    }

    # 変数の更新(.=を使うと感度がズレるので=を使うこと)
    eval.x = x

    # Compute function
    func = compute_eval_func(eval, physics, opt_settings, opt_filter)

    # Compute gradient by adjoint variable method
    df = compute_grad_by_adjoint(eval, physics, opt_settings, opt_filter)
    
    # Update for filtering and Finish
    return func, df
end