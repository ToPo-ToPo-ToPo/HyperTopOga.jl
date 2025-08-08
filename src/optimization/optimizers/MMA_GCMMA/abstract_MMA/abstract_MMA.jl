
#--------------------------------------------------------------------------------------------------------
# MMAアルゴリズムの抽象型
#--------------------------------------------------------------------------------------------------------
abstract type AbstractMMA end
#--------------------------------------------------------------------------------------------------------
# 目的関数の設定
#--------------------------------------------------------------------------------------------------------
function set_objective_function!(opt::O, target::String, eval::E, weight::Float64, df_method::M=AutomaticDiff()
    ) where {O<:AbstractMMA, E<:AbstractEvaluationFunction, M<:AbstractDfMethod}

    # 目的関数の追加
    push!(opt.objective, eval)
    push!(opt.objective0, 1.0)
    push!(opt.objective_df_methods, df_method)

    # 最小化問題の場合の重み
    if target == "min"
        weight_func_min(x) = weight
        push!(opt.objective_weight, weight_func_min)

    # 最大化問題の場合の重み
    elseif target == "max"
        weight_func_max(x) = -weight
        push!(opt.objective_weight, weight_func_max)
    
    # どちらにも該当しない場合
    else
        error("You can only choose min or max for the objective function setting.")
    end
end
#--------------------------------------------------------------------------------------------------------
# 制約関数の設定
#--------------------------------------------------------------------------------------------------------
function add_inequality_constraint!(opt::O, eval::E, limit::Float64, df_method::M=AutomaticDiff()
    ) where {O<:AbstractMMA, E<:AbstractEvaluationFunction, M<:AbstractDfMethod}
    # 制約関数の追加
    push!(opt.constraint, eval)
    push!(opt.constraint_lim, x -> limit)
    push!(opt.constraint_df_methods, df_method)
    opt.m += 1
end
#--------------------------------------------------------------------------------------------------------
# 目的関数、制約関数の計算+それぞれの感度の計算を行う
#--------------------------------------------------------------------------------------------------------
function update_filter_status!(::O, filter::F, opt_step::Int64) where {O<:AbstractMMA, F<:AbstractOptFilter}
    update_status!(filter, opt_step)
end