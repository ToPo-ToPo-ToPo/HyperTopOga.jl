
#----------------------------------------------------------------------------
# 体積評価関数
#----------------------------------------------------------------------------
mutable struct Volume <: AbstractEvaluationFunction
    x::Vector{Float64}
    value::Float64
    target::Int64           # 対象の材料番号
    #--------------------------------------------------------
    # 内部コンストラクタ
    #--------------------------------------------------------
    function Volume(x0::Vector{Float64}, target::Int64)
        return new(x0, 0.0, target)
    end
end
#-----------------------------------------------------------------------------------
# インターフェース
#-----------------------------------------------------------------------------------
function compute_eval_func(
    eval::Volume, 
    physics::AbstractPhysics, 
    opt_settings::AbstractOptSettings, 
    opt_filter::AbstractOptFilter
    )::Float64

    # フィルタリング
    x_tilde = filtering(opt_filter, eval.x)
    # 体積を計算する
    return compute_eval_volume(x_tilde, eval, physics, opt_settings, opt_settings.design_space_type)
end
#-----------------------------------------------------------------------------------
# 全体の体積を計算 (Element base)
#-----------------------------------------------------------------------------------
function compute_eval_volume(
    x_tilde::Vector{Float64}, 
    eval::Volume, 
    physics::AbstractPhysics, 
    opt_settings::AbstractOptSettings, 
    ::ElementBase
    )::Float64
    # 
    volume = 0.0
    volume_d = 0.0
    @floop for (e, element) in enumerate(physics.elements)
        # 位置eにおける全ての設計変数のベクトルを取得
        xe = make_local_design_variables(x_tilde,  opt_settings.num_svalue_types, e)
        # 積分点ループ
        for ip in 1 : element.shape.num_eval_points
            # 積分点座標を取得
            coordinate = element.shape.eval_coordinates[ip, :]
            # 形状関数の微分
            dNdxi = make_dNdxi(element.shape, coordinate)
            # Jacobi行列の計算
            J = make_J(element.shape, dNdxi)
            #
            factor = integral_factor(element.shape, ip) * det(J)
            #
            @reduce volume += xe[eval.target] * factor
            @reduce volume_d += factor
        end
    end
    # Stock result
    eval.value = volume / volume_d
    #
    return eval.value
end
#-----------------------------------------------------------------------------------
# 体積制約の感度のインターフェース
#-----------------------------------------------------------------------------------
function compute_grad_by_automatic_adjoint(
    eval::Volume, 
    physics::AbstractPhysics, 
    opt_settings::AbstractOptSettings, 
    opt_filter::AbstractOptFilter
    )::Vector{Float64}

    # フィルタリング
    x_tilde = filtering(opt_filter, eval.x)
    # 感度を計算
    dfdx_tilde = compute_eval_volume_grad_by_adjoint(x_tilde, eval, physics, opt_settings, opt_settings.design_space_type)
    # フィルタリングの影響を考慮
    return filtering_df(opt_filter, dfdx_tilde, eval.x)
end
#-----------------------------------------------------------------------------------
# 体積制約の感度 (Element base)
#-----------------------------------------------------------------------------------
function compute_eval_volume_grad_by_adjoint(
    x_tilde::Vector{Float64},
    eval::Volume, 
    physics::AbstractPhysics, 
    opt_settings::AbstractOptSettings, 
    ::ElementBase
    )::Vector{Float64}
    # 
    df = zero(x_tilde)
    volume_d = 0.0
    @floop for (e, element) in enumerate(physics.elements)
        # 設計変数のindexと位置eにおける設計変数の組を取得
        index = make_index(opt_settings.num_svalue_types, eval.target, e)
        # 積分点ループ
        local temp = 0.0
        for ip in 1 : element.shape.num_eval_points
            # 積分点座標を取得
            coordinate = element.shape.eval_coordinates[ip, :]
            # 形状関数の微分
            dNdxi = make_dNdxi(element.shape, coordinate)
            # Jacobi行列の計算
            J = make_J(element.shape, dNdxi)
            #
            factor = integral_factor(element.shape, ip) * det(J)
            #
            temp += factor
            @reduce volume_d += factor
        end
        df[index] = temp
    end
    df ./= volume_d
    # 
    return df
end