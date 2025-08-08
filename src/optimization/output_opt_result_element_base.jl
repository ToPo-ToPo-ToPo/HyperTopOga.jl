
#----------------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------------
function output_opt_result(
    x::Vector{Float64}, 
    physics::AbstractPhysics, 
    opt_settings::AbstractOptSettings, 
    opt_filter::AbstractOptFilter, 
    ::ElementBase
    )
    
    # フィルタリング
    x_tilde = filtering(opt_filter, x)
    
    #
    vtk_name = physics.output_file_name * "_final"

    # vtkファイルに結果を書き出し
    vtk_datasets = Vector{VtkDataset}()

    # データの作成
    svalue = Vector{Float64}(undef, length(physics.elements))
    topology = similar(svalue)
    Threads.@threads for e in 1 : length(physics.elements)
        xe1 = make_local_design_variables(x, opt_settings.num_svalue_types, e)
        xe2 = make_local_design_variables(x_tilde, opt_settings.num_svalue_types, e)
        svalue[e] = compute(opt_settings.interpolate, xe1)
        topology[e] = compute(opt_settings.interpolate, xe2)
    end

    # 出力データに追加
    push!(vtk_datasets, VtkDataset("design_variables", "CellData", svalue))
    push!(vtk_datasets, VtkDataset("topology", "CellData", topology))
    
    # 個別の設計変数に対する結果を出力
    for isvalue in 1 : opt_settings.num_svalue_types
        # 初期化
        s_i = Vector{Float64}(undef, length(physics.elements))
        t_i = similar(s_i)
        Threads.@threads for e in 1 : length(physics.elements)
            #
            index = make_index(opt_settings.num_svalue_types, isvalue, e)
            # 設計変数とトポロジー
            s_i[e] = x[index]
            t_i[e] = x_tilde[index]
        end
        # 出力データに追加
        push!(vtk_datasets, VtkDataset("design_variables" * string(isvalue), "CellData", s_i))
        push!(vtk_datasets, VtkDataset("topology" * string(isvalue), "CellData", t_i))
    end
    
    #
    output_vtu(physics, vtk_datasets, vtk_name)

end
#----------------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------------
function output_opt_result(
    x::Vector{Float64}, 
    physics, 
    opt_settings::AbstractOptSettings, 
    opt_filter::AbstractOptFilter, 
    df1::Vector{Float64}, 
    df2::Vector{Float64},
    ::ElementBase
    )
    
    # フィルタリング
    x_tilde = filtering(opt_filter, x)
    
    #
    vtk_name = physics.output_file_name * "_final"

    # vtkファイルに結果を書き出し
    vtk_datasets = Vector{VtkDataset}()
    
    # 個別の設計変数に対する結果を出力
    for isvalue in 1 : opt_settings.num_svalue_types
        # 初期化
        df1_i = Vector{Float64}(undef, length(physics.elements))
        df2_i = similar(df1_i)
        s_i = similar(df1_i)
        t_i = similar(s_i)
         for e in 1 : length(physics.elements)
            index = make_index(opt_settings.num_svalue_types, isvalue, e)
            # 感度
            df1_i[e] = df1[index]
            df2_i[e] = df2[index]
            # 設計変数とトポロジー
            s_i[e] = x[index]
            t_i[e] = x_tilde[index]
        end
        # 出力データに追加
        push!(vtk_datasets, VtkDataset("df_type1_" * string(isvalue), "CellData", df1_i))
        push!(vtk_datasets, VtkDataset("df_type2_" * string(isvalue), "CellData", df2_i))
        push!(vtk_datasets, VtkDataset("design_variables" * string(isvalue), "CellData", s_i))
        push!(vtk_datasets, VtkDataset("topology" * string(isvalue), "CellData", t_i))
    end
    
    #
    output_vtu(physics, vtk_datasets, vtk_name)

end
#----------------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------------
function output_opt_result(
    x::Vector{Float64}, 
    physics::AbstractPhysics, 
    opt_settings::AbstractOptSettings, 
    opt_filter::AbstractOptFilter, 
    df_auto::Vector{Float64}, 
    df_adjoint::Vector{Float64}, 
    df_numeric::Vector{Float64},
    ::ElementBase
    )
    
    # フィルタリング
    x_tilde = filtering(opt_filter, x)
    
    #
    vtk_name = physics.output_file_name * "_final"

    # vtkファイルに結果を書き出し
    vtk_datasets = Vector{VtkDataset}()
    
    # 個別の設計変数に対する結果を出力
    for isvalue in 1 : opt_settings.num_svalue_types
        # 初期化
        df_auto_i = Vector{Float64}(undef, length(physics.elements))
        df_adjoint_i = similar(df_auto_i)
        df_numeric_i = similar(df_auto_i)
        s_i = similar(df_auto_i)
        t_i = similar(s_i)
         for e in 1 : length(physics.elements)
            index = make_index(opt_settings.num_svalue_types, isvalue, e)
            # 感度
            df_auto_i[e] = df_auto[index]
            df_adjoint_i[e] = df_adjoint[index]
            df_numeric_i[e] = df_numeric[index]
            # 設計変数とトポロジー
            s_i[e] = x[index]
            t_i[e] = x_tilde[index]
        end
        # 出力データに追加
        push!(vtk_datasets, VtkDataset("df_auto" * string(isvalue), "CellData", df_auto_i))
        push!(vtk_datasets, VtkDataset("df_adjoint" * string(isvalue), "CellData", df_adjoint_i))
        push!(vtk_datasets, VtkDataset("df_numeric" * string(isvalue), "CellData", df_numeric_i))
        push!(vtk_datasets, VtkDataset("design_variables" * string(isvalue), "CellData", s_i))
        push!(vtk_datasets, VtkDataset("topology" * string(isvalue), "CellData", t_i))
    end
    
    #
    output_vtu(physics, vtk_datasets, vtk_name)

end