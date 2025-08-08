
#-------------------------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------------------------
struct AutomaticDiff <: AbstractDfMethod 
    name::String
    #---------------------------------------------------
    #
    #---------------------------------------------------
    function AutomaticDiff()
        return new("Automatic_diff")
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
    ::AutomaticDiff
    )::Float64 where {
        EF<:AbstractEvaluationFunction, 
        P<:AbstractPhysics, 
        OF<:AbstractOptFilter
    }

    # フィルタリング + 変数の更新(.=を使うと感度がズレるので=を使うこと)
    eval.x = x

    #
    return compute_eval_func(eval, physics, opt_settings, opt_filter)

end
#-------------------------------------------------------------------------------------------------
# 自動微分による感度解析
#-------------------------------------------------------------------------------------------------
function compute_f_and_grad(
    x::Vector{Float64}, 
    eval::EF, 
    physics::P, 
    opt_settings::OptimizationSettings, 
    opt_filter::OF, ::AutomaticDiff
    )::Tuple{Float64, Vector{Float64}} where {
        EF<:AbstractEvaluationFunction, 
        P<:AbstractPhysics, 
        OF<:AbstractOptFilter
    }

    # Zero clear    
    deval = make_zero(eval)

    # フィルタリング + 変数の更新(.=を使うと感度がズレるので=を使うこと)
    eval.x = x
    
    # Compute function and gradient
    autodiff(
        Enzyme.set_runtime_activity(Enzyme.Reverse), 
        Const(compute_eval_func), 
        Duplicated(eval, deval), 
        Const(physics), 
        Const(opt_settings),
        Const(opt_filter)
    )
    
    # Update for filtering and Finish
    return eval.value, deval.x

end