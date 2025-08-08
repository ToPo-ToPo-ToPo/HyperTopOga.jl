
#----------------------------------------------------------------
# 
#----------------------------------------------------------------
abstract type AbstractShape end
#----------------------------------------------------------------
# 要素中心座標を計算する関数
#----------------------------------------------------------------
function compute_center_coordinate(nodes::Vector{Node}, shape::S)::Vector{Float64} where {S<:AbstractShape}
    # 初期化
    center_coord = zeros(Float64, 3)
    # 主な処理
    for node_id in shape.connects
        center_coord += nodes[node_id].coordinate
    end
    # 平均化した値を返す
    return center_coord / Float64(length(shape.connects))
end
#----------------------------------------------------------------
# 要素中心座標を計算する関数
#----------------------------------------------------------------
function compute_boundary_center_coordinate(nodes::Vector{Node}, shape::S, iboundary::Int64)::Vector{Float64} where {S<:AbstractShape}
    # 初期化
    center_coord = zeros(Float64, 3)
    num_node = size(shape.boundary_connects)[2]
    # 主な処理
    for i in 1 : num_node
        node_id = shape.boundary_connects[iboundary, i]
        center_coord += nodes[node_id].coordinate
    end
    # 平均化した値を返す
    return center_coord / Float64(num_node)
end
#----------------------------------------------------------------
# 体積を計算する
#----------------------------------------------------------------
function compute_volume(shape::S)::Float64 where {S<:AbstractShape}
    #
    volume = 0.0
    # 積分点での計算
    for ip in 1 : shape.num_eval_points
        # 積分点座標を取得
        eval_coordinate = shape.eval_coordinates[ip, :]
        # 形状関数の微分
        dNdxi = make_dNdxi(shape, eval_coordinate)
        # Jacobi行列の計算
        J = make_J(shape, dNdxi)
        #
        factor = integral_factor(shape, ip) * my_det(J)
        #
        volume += factor
    end
    return volume
end