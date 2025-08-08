

#------------------------------------------------------------------------
# 
#------------------------------------------------------------------------
function output_vtu_for_mma_gcmma(
    physics::P, 
    opt_settings::OS, 
    opt_filter::F, 
    xval::Vector{Float64}, 
    df0dx::Vector{Float64}, 
    vtk_name::String,
    design_space::DS
    ) where {P<:AbstractPhysics, OS<:AbstractOptSettings, F<:AbstractOptFilter, DS<:AbstractDesignSpace}

    # vtkファイルに結果を書き出し
    vtk_datasets = Vector{VtkDataset}()

    # フィルタリング
    x_tilde = filtering(opt_filter, xval)

    # 全体のトポロジーのデータを出力に追加
    add_design_variable_data!(vtk_datasets, physics, opt_settings, xval, x_tilde, design_space)

    # 個別の設計変数に対する結果を出力
    add_each_design_variable_data!(vtk_datasets, physics, opt_settings, xval, x_tilde, df0dx, design_space)

    #
    output_vtu(physics, vtk_datasets, vtk_name)

end
#------------------------------------------------------------------------
# ベクトル版
#------------------------------------------------------------------------
function output_vtu_for_mma_gcmma(
    physics::Vector{P}, 
    opt_settings::OS, 
    opt_filter::F, 
    xval::Vector{Float64}, 
    df0dx::Vector{Float64}, 
    vtk_name::String,
    design_space::DS
    ) where {P<:AbstractPhysics, OS<:AbstractOptSettings, F<:AbstractOptFilter, DS<:AbstractDesignSpace}

    # vtkファイルに結果を書き出し
    vtk_datasets = Vector{VtkDataset}()

    # フィルタリング
    x_tilde = filtering(opt_filter, xval)

    # 全体のトポロジーのデータを出力に追加
    # 1つ目の物理モデルのデータを使う
    add_design_variable_data!(vtk_datasets, physics[1], opt_settings, xval, x_tilde, design_space)

    # 個別の設計変数に対する結果を出力
    # 1つ目の物理モデルのデータを使う
    add_each_design_variable_data!(vtk_datasets, physics[1], opt_settings, xval, x_tilde, df0dx, design_space)

    #
    # 1つ目の物理モデルのデータを使う
    output_vtu(physics[1], vtk_datasets, vtk_name)

end