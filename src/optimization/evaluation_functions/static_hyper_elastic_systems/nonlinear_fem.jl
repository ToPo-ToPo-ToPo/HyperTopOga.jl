
#-------------------------------------------------------------------------------------------------
# 非線形FEMにより、最終ステップの変位を求める
#-------------------------------------------------------------------------------------------------
function compute_nonlinear_fem!(
    eval::AbstractHyperElasticEvalFunc, 
    physics::StaticSolidMechanics, 
    element_states::Vector{ElementState}, 
    opt_settings::OptimizationSettings,
    opt_filter::AbstractOptFilter
    )::Nothing

    # 変位の初期化
    eval.U .= zero(eval.U)

    # vtkファイルに結果を書き出し
    if opt_settings.output_flag == true
        output_vtu(eval, physics, element_states, opt_settings, opt_filter, 0, opt_settings.design_space_type)
    end

    # 増分解析ループ
    for istep in 1 : eval.num_step
        # NR法により支配方程式の解を求める
        eval.U .= newton_rapson(eval.method, physics, element_states, istep, eval.U, opt_settings.output_flag)
        
        # vtkファイルに結果を書き出し
        if opt_settings.output_flag == true
            output_vtu(eval, physics, element_states, opt_settings, opt_filter, istep, opt_settings.design_space_type)
        end
    end

    return nothing

end