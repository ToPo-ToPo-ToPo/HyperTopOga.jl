
#----------------------------------------------------------------
# VTKファイル出力に関するデータをまとめた構造体
#----------------------------------------------------------------
struct VtkDataset
    name::String             # データのラベル
    type                     # データの出力形式 point or cell
    field                    # 出力する変数
end
#----------------------------------------------------------------
# Paraviewのvtkファイルに結果を書き出す関数
#----------------------------------------------------------------
function output_vtu(physics::P, vtk_datasets::Vector{VtkDataset}, filename::String)::Nothing where {P<:SolidMechanics}
    
    # Set
    nodes = physics.nodes
    elements = physics.elements

    # 節点情報の設定
    num_point = length(nodes) 
    points = Matrix{Float64}(undef, 3, num_point)
    for i in 1 : 3
        for node in nodes
            id = node.id
            points[i, id] = node.coordinate[i]
        end
    end

    # 要素タイプの設定
    num_cell = length(elements)
    cells = Vector{MeshCell}(undef, num_cell)
    for element in elements
        id = element.id
        cells[id] = get_shape_data(element.shape)
    end

    # 結果の出力
    vtk_grid(filename, points, cells) do vtk
        for vtk_dataset in vtk_datasets
            if vtk_dataset.type == "PointData"
                vtk[vtk_dataset.name, VTKPointData()] = convert_data_to_output_format(nodes, vtk_dataset.field)
            elseif vtk_dataset.type == "ScalarPointData"
                    vtk[vtk_dataset.name, VTKPointData()] = convert_data_to_output_format_scalar(nodes, vtk_dataset.field)
            elseif vtk_dataset.type == "CellData"
                vtk[vtk_dataset.name, VTKCellData()] = convert_data_to_output_format(elements, vtk_dataset.field)
            else
                throw(ErrorException("No applicable dataset type, please set PointData or CellData."))
            end
        end
    end

    return nothing
end
#----------------------------------------------------------------
# VTKファイル出力用に該当する節点情報に一致するように中身を変更
#----------------------------------------------------------------
function convert_data_to_output_format(nodes::Vector{Node}, values::Vector{Float64})::Matrix{Float64}
    # 出力用のデータを作成
    output_values = zeros(Float64, 3, length(nodes))
    for (inode, node) in enumerate(nodes)
        for (i, idof) in enumerate(node.dof_ids)
            output_values[i, inode] = values[idof]
        end
    end
    return output_values
end
#----------------------------------------------------------------
# VTKファイル出力用に該当する節点情報に一致するように中身を変更
#----------------------------------------------------------------
function convert_data_to_output_format_scalar(::Vector{Node}, values::Vector{Float64})::Vector{Float64}
    # 出力用のデータを作成
    return values
end
#----------------------------------------------------------------
# VTKファイル出力用に該当する節点情報に一致するように中身を変更
#----------------------------------------------------------------
function convert_data_to_output_format(elements::EV, values::Vector{Float64})::Vector{Float64} where {EV<:Vector{StructuralElement}}
    return values
end
