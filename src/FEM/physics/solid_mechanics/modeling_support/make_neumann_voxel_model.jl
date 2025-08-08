

#----------------------------------------------------------------
# ボクセル要素向けの境界条件を設定するプログラム
# 入力した範囲から該当する節点を取得し、境界条件のデータを返す
#----------------------------------------------------------------
function make_neumann_voxel_model(
    nodes::Vector{Node}, 
    elements::Vector{E}, 
    range_min::Vector{Float64}, 
    range_max::Vector{Float64},
    num_step::Int64,   
    values::Vector{Function}
    )::NeumannBc where {E<:StructuralElement}

    # 初期化
    target_element_ids = Vector{Int64}()
    target_boundary_ids = Vector{Int64}()
    
    # 要素ループを行い、該当する面を探索
    for element in elements
        # 要素境界ループ
        for iboundary in 1 : element.shape.num_boundarys
            # 情報を取得
            boundary_id = iboundary
            # 要素中心を取得
            center_coord = compute_boundary_center_coordinate(nodes, element.shape, iboundary)
            # 範囲に含まれるかを確認
            if range_min[1] - 1.0e-05 < center_coord[1] < range_max[1] + 1.0e-05 && 
               range_min[2] - 1.0e-05 < center_coord[2] < range_max[2] + 1.0e-05 &&
               range_min[3] - 1.0e-05 < center_coord[3] < range_max[3] + 1.0e-05
                   # 該当する要素を追加
                   push!(target_element_ids, element.id)
                   # 回答する要素内境界番号を追加
                   push!(target_boundary_ids, boundary_id)
            end
        end
    end
    # 条件値の設定
    value_set = zeros(Float64, length(values), num_step)
    for istep in 1 : num_step
        for j in 1 : length(values)
            value_set[j, istep] = values[j](istep)
        end
    end
    #
    return NeumannBc(target_element_ids, target_boundary_ids, value_set)
end
#----------------------------------------------------------------
# ボクセル要素向けの境界条件を設定するプログラム
# 入力した範囲から該当する節点を取得し、境界条件のデータを返す
#----------------------------------------------------------------
function make_neumann_voxel_model(
    nodes::Vector{Node}, 
    range_min::Vector{Float64}, 
    range_max::Vector{Float64},
    num_step::Int64,  
    values::Vector{Function}
    )::NeumannPointBc

    # 初期化
    target_node_ids = Vector{Int64}()
    # 要素ループを行い、該当する面を探索
    for node in nodes
        # 座標を取得
        coord = node.coordinate
        # 範囲に含まれるかを確認
        if range_min[1] - 1.0e-05 < coord[1] < range_max[1] + 1.0e-05 && 
            range_min[2] - 1.0e-05 < coord[2] < range_max[2] + 1.0e-05 &&
            range_min[3] - 1.0e-05 < coord[3] < range_max[3] + 1.0e-05
            # 該当する要素を追加
            push!(target_node_ids, node.id)
        end
    end
    # 条件値の設定
    value_set = zeros(Float64, length(values), num_step)
    for istep in 1 : num_step
        for j in 1 : length(values)
            value_set[j, istep] = values[j](istep)
        end
    end
    #
    return NeumannPointBc(target_node_ids, value_set)
end
#----------------------------------------------------------------------------
# 条件を更新（点荷重用）
#----------------------------------------------------------------------------
function update_neumann_voxel_model!(
    bc::AbstractNeumannBc, 
    nodes::Vector{Node}, 
    range_min::Vector{Float64}, 
    range_max::Vector{Float64},
    num_step::Int64,  
    values::Vector{Function}
    )::Nothing

    # 初期化
    target_node_ids = Vector{Int64}()
    # 要素ループを行い、該当する面を探索
    for node in nodes
        # 座標を取得
        coord = node.coordinate
        # 範囲に含まれるかを確認
        if range_min[1] - 1.0e-05 < coord[1] < range_max[1] + 1.0e-05 && 
            range_min[2] - 1.0e-05 < coord[2] < range_max[2] + 1.0e-05 &&
            range_min[3] - 1.0e-05 < coord[3] < range_max[3] + 1.0e-05
            # 該当する要素を追加
            push!(target_node_ids, node.id)
        end
    end
    # 条件値の設定
    value_set = zeros(Float64, length(values), num_step)
    for istep in 1 : num_step
        for j in 1 : length(values)
            value_set[j, istep] = values[j](istep)
        end
    end
    # Update
    bc.node_ids = target_node_ids
    bc.values = value_set
    #
    return nothing
end