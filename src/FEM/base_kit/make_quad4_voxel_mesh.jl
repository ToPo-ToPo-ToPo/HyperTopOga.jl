
#----------------------------------------------------------------
# モデルの長さ、分割数からボクセルメッシュを作成
# 2D Quad4 version
#----------------------------------------------------------------
function make_quad4_voxel_mesh(length_x::Float64, length_y::Float64, division_x::Int64, division_y::Int64)
    # 節点情報の作成
    nodes,num_elements = make_nodes(length_x, length_y, division_x, division_y)
    # 要素の定義
    connectivities = make_quad4_connectivities(nodes, division_x, division_y)
    # 作成した構造体を返す
    return nodes, connectivities
end
#----------------------------------------------------------------
# モデルの長さ、分割数から節点情報を作成
#----------------------------------------------------------------
function make_nodes(length_x::Float64, length_y::Float64, division_x::Int64, division_y::Int64)
    # 初期化
    num_node = (division_x + 1) * (division_y + 1)
    nodes = Vector{Node}(undef, num_node)
    num_elements=division_x  * division_y 
    # 節点間距離の計算
    delta_x = length_x / Float64(division_x)
    delta_y = length_y / Float64(division_y)
    # 分割数のループ
    for j = 1 : division_y+1
        # y座標
        coord_y = Float64(j - 1) * delta_y
        for i = 1 : division_x+1
            # x座標
            coord_x = Float64(i - 1) * delta_x
            # 節点番号
            id = i + (j - 1) + division_x * (j - 1)
            # 情報の更新
            nodes[id] = Node(id, [coord_x, coord_y, 0.0]) 
        end
    end
    return nodes, num_elements
end
#----------------------------------------------------------------
# モデルの長さ、分割数から四角形要素の情報を作成
#----------------------------------------------------------------
function make_quad4_connectivities(nodes::Vector{Node}, division_x::Int64, division_y::Int64)
    # 初期化
    num_connectivity = division_x * division_y
    connectivities = Matrix{Node}(undef, num_connectivity, 4)
    for j = 1 : division_y
        for i = 1 : division_x
            id = i + (j - 1) * division_x
            n1 = i + (j - 1) * (division_x + 1)
            n2 = n1 + 1
            n3 = n1 + division_x + 1 + 1
            n4 = n1 + division_x + 1
            connectivities[id, :] = [nodes[n1] nodes[n2] nodes[n3] nodes[n4]]
        end
    end
    return connectivities
end