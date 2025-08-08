
#----------------------------------------------------------------
# 
#----------------------------------------------------------------
mutable struct Quad4 <: AbstractShape 
    #
    dimension::Int64
    connects::Vector{Int64}                     # 構成する節点番号
    coord_matrix::Matrix{Float64}               # 要素内節点の座標を集めたマトリクス
    
    #
    num_eval_points::Int64                      # 評価点の数
    eval_coordinates::Matrix{Float64}           # 積分点の座標
    eval_weights::Vector{Float64}               # 積分点の重み
    thickness::Float64                          # 要素の厚さ
    
    # 境界に関する情報
    num_boundarys::Int64
    boundary_connects::Matrix{Int64}            # 境界を構成する節点番号
    
    #
    num_boundary_eval_points::Int64             # 境界の積分点の数
    boundary_eval_coordinaites::Array{Float64}  # 境界の評価点の座標
    boundary_eval_weights::Matrix{Float64}
    #----------------------------------------------------------------
    #
    #----------------------------------------------------------------
    function Quad4(nodes::Vector{Node}, thickness::Float64)
        
        # コネクティビティ
        connects = [nodes[1].id, nodes[2].id, nodes[3].id, nodes[4].id]

        # 座標をまとめたマトリクスを作成
        coord_matrix = Matrix{Float64}(undef, 4, 2)
        for (i, node) in enumerate(nodes)
            for j in 1 : 2
                coord_matrix[i, j] = node.coordinate[j]
            end
        end
        
        # 積分点に関する情報
        num_eval_points = 4
        eval_coordinates = Matrix{Float64}(undef, num_eval_points, 2)
        eval_coordinates[1, :] = [-1.0 / sqrt(3.0), -1.0 / sqrt(3.0)]
        eval_coordinates[2, :] = [ 1.0 / sqrt(3.0), -1.0 / sqrt(3.0)]
        eval_coordinates[3, :] = [ 1.0 / sqrt(3.0),  1.0 / sqrt(3.0)]
        eval_coordinates[4, :] = [-1.0 / sqrt(3.0),  1.0 / sqrt(3.0)]
        #
        eval_weights = fill(1.0, num_eval_points)

        # 境界に関する情報
        num_boundarys = 4
        boundary_connects = Matrix{Int64}(undef, num_boundarys, 2)
        boundary_connects[1, :] = [nodes[1].id, nodes[2].id]
        boundary_connects[2, :] = [nodes[2].id, nodes[3].id]
        boundary_connects[3, :] = [nodes[3].id, nodes[4].id]
        boundary_connects[4, :] = [nodes[4].id, nodes[1].id]
        
        # 境界積分に関する情報
        num_boundary_eval_points = 2
        #
        boundary_eval_coordinaites = Array{Float64}(undef, num_boundarys, num_boundary_eval_points, 2)
        boundary_eval_coordinaites[1, 1, :] = [-1.0 / sqrt(3.0), -1.0            ]
        boundary_eval_coordinaites[1, 2, :] = [ 1.0 / sqrt(3.0), -1.0            ]
        boundary_eval_coordinaites[2, 1, :] = [ 1.0,             -1.0 / sqrt(3.0)]
        boundary_eval_coordinaites[2, 2, :] = [ 1.0,              1.0 / sqrt(3.0)]
        boundary_eval_coordinaites[3, 1, :] = [-1.0 / sqrt(3.0),  1.0            ]
        boundary_eval_coordinaites[3, 2, :] = [ 1.0 / sqrt(3.0),  1.0            ]
        boundary_eval_coordinaites[4, 1, :] = [-1.0,             -1.0 / sqrt(3.0)]
        boundary_eval_coordinaites[4, 2, :] = [-1.0,              1.0 / sqrt(3.0)]
        #
        boundary_eval_weights = fill(1.0, num_boundarys, num_boundary_eval_points)
        
        # 構造体の作成
        return new(
            2, connects, coord_matrix, 
            num_eval_points, eval_coordinates, eval_weights, 
            thickness, 
            num_boundarys, boundary_connects, 
            num_boundary_eval_points, boundary_eval_coordinaites, boundary_eval_weights
        )
    end
end
#----------------------------------------------------------------
# 形状関数の定義
#----------------------------------------------------------------
function make_N_components(::Quad4, coordinate::Vector{Float64})::Vector{Float64}
    # 積分点座標の設定
    xi, et = coordinate
    # 成分の作成
    N = Vector{Float64}(undef, 4)
    N[1] = 0.25 * (1.0 - xi) * (1.0 - et)
    N[2] = 0.25 * (1.0 + xi) * (1.0 - et)
    N[3] = 0.25 * (1.0 + xi) * (1.0 + et)
    N[4] = 0.25 * (1.0 - xi) * (1.0 + et)
    return N
end
#----------------------------------------------------------------
# 形状関数の導関数の計算
#----------------------------------------------------------------
function make_dNdxi(::Quad4, coordinate::Vector{Float64})::Matrix{Float64}
    # 積分点座標の設定
    xi, et = coordinate
    # 形状関数の微分を計算
    dNdxi = Matrix{Float64}(undef, 2, 4)
    dNdxi[1, :] = [-0.25 * (1.0 - et) 0.25 * (1.0 - et) 0.25 * (1.0 + et) -0.25 * (1.0 + et)]
    dNdxi[2, :] = [-0.25 * (1.0 - xi) -0.25 * (1.0 + xi) 0.25 * (1.0 + xi) 0.25 * (1.0 - xi)]
    return dNdxi
end
#----------------------------------------------------------------
# Jacobi行列を計算
#----------------------------------------------------------------
function make_J(shape::Quad4, dNdxi::Matrix{Float64})::Matrix{Float64}
    return dNdxi * shape.coord_matrix
end
#----------------------------------------------------------------
# dNdxの計算: dNdxi = J * dNdx -> dNdx = J_invers * dNdxi
#----------------------------------------------------------------
function make_dNdx(::Quad4, dNdxi::Matrix{Float64}, J::Matrix{Float64})::Matrix{Float64}
    return my_inv22(J) * dNdxi
end
#----------------------------------------------------------------
# 積分の重み係数
#----------------------------------------------------------------
function integral_factor(shape::Quad4, ip::Int64)::Float64
    return shape.thickness * shape.eval_weights[ip]
end
#----------------------------------------------------------------
# 積分の重み係数
#----------------------------------------------------------------
function boundary_integral_factor(shape::Quad4, boundary_id::Int64, ip::Int64)::Float64
    return shape.thickness * shape.boundary_eval_weights[boundary_id, ip]
end
#----------------------------------------------------------------
# 外接円の半径を計算
#----------------------------------------------------------------
function circumscribed_circle_radius(shape::Quad4)::Float64
    #
    x1, y1 = shape.coord_matrix[1, :]
    x2, y2 = shape.coord_matrix[2, :]
    x3, y3 = shape.coord_matrix[3, :]
    x4, y4 = shape.coord_matrix[4, :]
    # 四角形の対角線の長さを求める
    diagonal1 = sqrt((x3 - x1)^2.0 + (y3 - y1)^2.0)
    diagonal2 = sqrt((x4 - x2)^2.0 + (y4 - y2)^2.0)
    # 対角線の長さの最大値を外接円の直径とする
    diameter = max(diagonal1, diagonal2)
    # 外接円の半径は直径の半分
    radius = diameter / 2.0
    #
    return radius
end
#----------------------------------------------------------------
#　VTKファイル出力時の要素形状の設定
#----------------------------------------------------------------
function get_shape_data(shape::Quad4)
    return MeshCell(VTKCellTypes.VTK_QUAD, shape.connects)
end