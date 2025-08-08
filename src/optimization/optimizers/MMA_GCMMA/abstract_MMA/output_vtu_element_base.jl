
#------------------------------------------------------------------------
#
#------------------------------------------------------------------------
function add_design_variable_data!(
    vtk_datasets::Vector{VtkDataset},
    physics::P,
    opt_settings::OS,
    xval::Vector{Float64},
    x_tilde::Vector{Float64},
    ::ElementBase
    ) where {P<:AbstractPhysics, OS<:AbstractOptSettings}

    # データの作成
    svalue = Vector{Float64}(undef, length(physics.elements))
    topology = similar(svalue)
    Threads.@threads for e in eachindex(physics.elements)
        xe1 = make_local_design_variables(xval, opt_settings.num_svalue_types, e)
        xe2 = make_local_design_variables(x_tilde, opt_settings.num_svalue_types, e)
        svalue[e] = compute(opt_settings.interpolate, xe1)
        topology[e] = compute(opt_settings.interpolate, xe2)
    end

    # 出力データに追加
    push!(vtk_datasets, VtkDataset("design_variables", "CellData", svalue))
    push!(vtk_datasets, VtkDataset("topology", "CellData", topology))

end
#------------------------------------------------------------------------
#
#------------------------------------------------------------------------
function add_each_design_variable_data!(
    vtk_datasets::Vector{VtkDataset},
    physics::P,
    opt_settings::OS,
    xval::Vector{Float64},
    x_tilde::Vector{Float64},
    df0dx::Vector{Float64},
    ::ElementBase
    ) where {P<:AbstractPhysics, OS<:AbstractOptSettings}

    # 個別の設計変数に対する結果を出力
    for isvalue in 1 : opt_settings.num_svalue_types
        # 初期化
        obj_sens = Vector{Float64}(undef, length(physics.elements))
        s_i = similar(obj_sens)
        t_i = similar(s_i)
        Threads.@threads for e in eachindex(physics.elements)
            #
            index = make_index(opt_settings.num_svalue_types, isvalue, e)
            # 目的関数の感度
            obj_sens[e] = df0dx[index]
            # 設計変数とトポロジー
            s_i[e] = xval[index]
            t_i[e] = x_tilde[index]
        end
        
        # 出力データに追加
        push!(vtk_datasets, VtkDataset("obj_sens" * string(isvalue), "CellData", obj_sens))
        push!(vtk_datasets, VtkDataset("design_variables" * string(isvalue), "CellData", s_i))
        push!(vtk_datasets, VtkDataset("topology" * string(isvalue), "CellData", t_i))
    end

end