

#-----------------------------------------------------------------------------------------------------------
# フィルターの1要素あたりのデータをまとめた構造体
#-----------------------------------------------------------------------------------------------------------
struct OptFilterDataset
    radius::Float64                    # フィルター半径
    num_ids::Vector{Int64}             # 1つのフィルタリングに含まれる数
    id_list::Matrix{Int64}             # フィルタリングに関係するデータのリスト
    weights::Matrix{Float64}           # 重み係数
    volumes::Matrix{Float64}       
    volume_self::Vector{Float64}
end
#-----------------------------------------------------------------------------------------------------------
# データを作成
#-----------------------------------------------------------------------------------------------------------
function setup_opt_filter(
    nodes::Vector{Node}, 
    elements::Vector{E}, 
    radius_coeff::Float64, 
    ::ElementBase
    )::OptFilterDataset where {E<:AbstractElement}

    # フィルタリング半径を更新
    radius = set_filtering_radius(elements, radius_coeff)

    # 要素数と次元数を取得
    num_elements = length(elements)

    # 要素中心座標の D×N 行列を構築
    centers = zeros(Float64, 3, num_elements)
    for (i, element) in enumerate(elements)
        centers[:, i] = compute_center_coordinate(nodes, element.shape)
    end

    # KDTree の構築
    tree = KDTree(centers)

    # 全ての近傍数を事前に計算し、最大数を取得
    neighbor_lists = [inrange(tree, centers[:, i], radius) for i in 1 : num_elements]
    max_comp = maximum(length.(neighbor_lists))

    # 結果格納用配列（固定長）
    num_ids = zeros(Int64, num_elements)
    id_list = zeros(Int64, num_elements, max_comp)
    weights = zeros(Float64, num_elements, max_comp)
    volumes = zeros(Float64, num_elements, max_comp)
    volume_self = zeros(Float64, num_elements)

    # 要素ループ
    @floop for i in 1 : num_elements
        # Get element i
        center_coord_i = centers[:, i]
        volume_self[i] = compute_volume(elements[i].shape)
        # ソートして順序を安定化
        idxs = sort(neighbor_lists[i])
        # 数を記録
        num_ids[i] = length(idxs)
        # それ以外（j > length(idxs)）は初期化済みの0が残る（padding）
        for j in 1 : num_ids[i]
            #
            idx = idxs[j]
            id_list[i, j] = idx
            # Get element j
            element_j = elements[idx]
            volume_j = compute_volume(element_j.shape)
            center_coord_j = centers[:, idx]
            dist = norm(center_coord_i .- center_coord_j)
            #
            weights[i, j] = volume_j * (radius - dist) / radius
            volumes[i, j] = volume_j
        end
    end

    return OptFilterDataset(radius, num_ids, id_list, weights, volumes, volume_self)
end
#-----------------------------------------------------------------------------------------------------------
# フィルタリング半径を作成
#-----------------------------------------------------------------------------------------------------------
function set_filtering_radius(elements::Vector{E}, radius_coeff::Float64)::Float64 where {E<:AbstractElement}
    # フィルター半径を作成
    radius_min = 1.0e+10
    for element in elements
        # メッシュサイズを取得
        mesh_size = circumscribed_circle_radius(element.shape)
        # update
        if radius_min > mesh_size 
            radius_min = mesh_size
        end
    end
    return radius_coeff * radius_min
end
#-----------------------------------------------------------------------------------------------------------
# 自分自身を含めたフィルタリンググループの作成
#-----------------------------------------------------------------------------------------------------------
function find_neighbors_kdtree_with_self(centers::Matrix{Float64}, radius::Float64)
    D, N = size(centers)
    tree = KDTree(centers)
    neighbors = Vector{Vector{Int}}(undef, N)

    for i in 1 : N
        # 自分自身も含めて取得
        neighbors[i] = inrange(tree, centers[:, i], radius)
    end

    return neighbors
end