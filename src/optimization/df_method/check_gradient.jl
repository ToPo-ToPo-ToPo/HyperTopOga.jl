#------------------------------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------------------------------
function check_gradient(
    x::Vector{Float64}, 
    eval::AbstractEvaluationFunction, 
    physics::AbstractPhysics, 
    opt_settings::AbstractOptSettings,
    opt_filter::AbstractOptFilter,
    df_method1::AbstractDfMethod,
    df_method2::AbstractDfMethod
    )

    # Set
    num_svalue_types = opt_settings.num_svalue_types

    # 微分1
    println(df_method1.name)
    @CPUtime f0, df1 = compute_f_and_grad(x, eval, physics, opt_settings, opt_filter, df_method1)
    df1 .= df1 ./ f0

    # 結果を出力しないようにフラグを更新
    opt_settings.output_flag = false

    # 微分2
    println(df_method2.name)
    @CPUtime f0, df2 = compute_f_and_grad(x, eval, physics, opt_settings, opt_filter, df_method2)
    df2 .= df2 ./ f0

    # output
    sensitivity_check(x, df1, df2)
    println("評価関数値: ", f0) 

    # テキストファイルにデータを出力する
    output_gradient_data("output/check", df_method1.name, num_svalue_types, opt_filter.num_svalue_domain, df1)
    output_gradient_data("output/check", df_method2.name, num_svalue_types, opt_filter.num_svalue_domain, df2)

    # 最終的な最適化結果の出力
    #output_opt_result(x, physics, opt_settings, opt_filter, df1, df2, opt_settings.design_space_type)

end
#------------------------------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------------------------------
function check_gradient(
    x::Vector{Float64}, 
    eval::AbstractEvaluationFunction, 
    physics::AbstractPhysics, 
    opt_settings::AbstractOptSettings, 
    opt_filter::AbstractOptFilter,
    df_method1::AbstractDfMethod,
    df_method2::AbstractDfMethod,
    df_method3::AbstractDfMethod
    )

    # Set
    num_svalue_types = opt_settings.num_svalue_types

    # 自動微分 リバースモード
    println(df_method1.name)
    @CPUtime f0, df1 = compute_f_and_grad(x, eval, physics, opt_settings, opt_filter, df_method1)
    df1 .= df1 ./ f0

    # 解析的微分
    println(df_method2.name)
    @CPUtime f0, df2 = compute_f_and_grad(x, eval, physics, opt_settings, opt_filter, df_method2)
    df2 .= df2 ./ f0

    # 数値微分
    println(df_method3.name)
    @CPUtime f0, df3 = compute_f_and_grad(x, eval, physics, opt_settings, opt_filter, df_method3)
    df3 .= df3 ./ f0

    # output
    sensitivity_check(x, df1, df2, df3)
    println("評価関数値: ", f0) 

    # テキストファイルにデータを出力する
    output_gradient_data("output/check", df_method1.name, num_svalue_types, opt_filter.num_svalue_domain, df1)
    output_gradient_data("output/check", df_method2.name, num_svalue_types, opt_filter.num_svalue_domain, df2)
    output_gradient_data("output/check", df_method3.name, num_svalue_types, opt_filter.num_svalue_domain, df3)

    # 最終的な最適化結果の出力
    output_opt_result(x, physics, opt_settings, opt_filter, df1, df2, df3, opt_settings.design_space_type)

end

#------------------------------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------------------------------
function check_gradient_memory(
    x::Vector{Float64}, 
    eval::AbstractEvaluationFunction, 
    physics::AbstractPhysics, 
    opt_settings::AbstractOptSettings,
    opt_filter::AbstractOptFilter,
    df_method1::AbstractDfMethod
    )

    # Set
    num_svalue_types = opt_settings.num_svalue_types

    # 微分1
    println(df_method1.name)
    f0, df1 = compute_f_and_grad(x, eval, physics, opt_settings, opt_filter, df_method1)
    

end
#------------------------------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------------------------------
function check_gradient(
    x::Vector{Float64}, 
    eval::AbstractEvaluationFunction, 
    physics::AbstractPhysics, 
    opt_settings::AbstractOptSettings, 
    opt_filter::AbstractOptFilter,
    df_method1::AbstractDfMethod,
    df_method2::AbstractDfMethod,
    df_method3::AbstractDfMethod,
    df_method4::AbstractDfMethod
    )

    # Set
    num_svalue_types = opt_settings.num_svalue_types

    # 自動微分 リバースモード
    println(df_method1.name)
    @CPUtime f0, df1 = compute_f_and_grad(x, eval, physics, opt_settings, opt_filter, df_method1)
    df1 .= df1 ./ f0

    # 解析的微分
    println(df_method2.name)
    @CPUtime f0, df2 = compute_f_and_grad(x, eval, physics, opt_settings, opt_filter, df_method2)
    df2 .= df2 ./ f0

    # 解析的微分
    println(df_method3.name)
    @CPUtime f0, df3 = compute_f_and_grad(x, eval, physics, opt_settings, opt_filter, df_method3)
    df3 .= df3 ./ f0

    # 数値微分
    println(df_method4.name)
    @CPUtime f0, df4 = compute_f_and_grad(x, eval, physics, opt_settings, opt_filter, df_method4)
    df4 .= df4 ./ f0

    # テキストファイルにデータを出力する
    output_gradient_data("output/check", df_method1.name, num_svalue_types, opt_filter.num_svalue_domain, df1)
    output_gradient_data("output/check", df_method2.name, num_svalue_types, opt_filter.num_svalue_domain, df2)
    output_gradient_data("output/check", df_method3.name, num_svalue_types, opt_filter.num_svalue_domain, df3)

    # 最終的な最適化結果の出力
    output_opt_result(x, physics, opt_settings, opt_filter, df1, df2, df3, opt_settings.design_space_type)

end