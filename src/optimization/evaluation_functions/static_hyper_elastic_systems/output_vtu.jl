
#----------------------------------------------------------------------------
# 解析結果の出力(For element base TO)
#----------------------------------------------------------------------------
function output_vtu(
    eval::AbstractHyperElasticEvalFunc,
    physics::StaticSolidMechanics, 
    element_states::Vector{ElementState}, 
    opt_settings::OptimizationSettings,
    opt_filter::AbstractOptFilter,
    step::Int64,
    ::ElementBase
    )::Nothing

    # Set vtu name
    name = physics.output_file_name * "_res_step_" * string(step)

    # Set
    U = eval.U
    x = eval.x
    x_tilde = filtering(opt_filter, x)
    
    # 解析結果の情報作成
    mises = Vector{Float64}(undef, length(physics.elements))
    topology = similar(mises)
    
    # Element loop
    @floop for (e, element) in enumerate(physics.elements)
        
        # 要素変位を作成
        Ue = make_element_sol(element.dof_list, U)

        # ミーゼス応力の計算
        local mises_local = 0.0
        for (ip, material) in enumerate(element_states[e].materials)
            mises_local += compute_mises_stress(element, material, ip, Ue, element.type)
        end
        mises[e] = mises_local / Float64(element.shape.num_eval_points)

        # トポロジーの計算
        x_tilde_e = make_local_design_variables(x_tilde, opt_settings.num_svalue_types, e)
        topology[e] = compute(opt_settings.interpolate, x_tilde_e)
    end
    
    # 初期化    
    vtk_datasets = Vector{VtkDataset}()
    # データを追加
    push!(vtk_datasets, VtkDataset("displacement", "PointData", U))
    push!(vtk_datasets, VtkDataset("mises", "CellData", mises))
    push!(vtk_datasets, VtkDataset("topology", "CellData", topology))
    # 結果の出力
    output_vtu(physics, vtk_datasets, name)

    return nothing
    
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(output_vtu), 
    eval::AbstractHyperElasticEvalFunc,
    physics::StaticSolidMechanics, 
    element_states::Vector{ElementState},
    opt_settings::OptimizationSettings,
    opt_filter::AbstractOptFilter,
    step::Int64,
    ::ElementBase
    )
    return nothing
end
