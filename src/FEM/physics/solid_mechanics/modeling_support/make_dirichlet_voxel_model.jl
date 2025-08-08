

#----------------------------------------------------------------
# ボクセル要素向けの境界条件を設定するプログラム
# 入力した範囲から該当する節点を取得し、境界条件のデータを返す
#----------------------------------------------------------------
function make_dirichlet_voxel_model(
    nodes::Vector{Node}, 
    range_min::Vector{Float64}, 
    range_max::Vector{Float64},
    num_step::Int64,  
    flag::Vector{Bool}, 
    values::Vector{Function}
    )::DirichletBc

    # 該当する節点を探索
    node_ids = Vector{Int64}()
    for node in nodes
        if range_min[1] - 1.0e-05 < node.coordinate[1] < range_max[1] + 1.0e-05 && 
           range_min[2] - 1.0e-05 < node.coordinate[2] < range_max[2] + 1.0e-05 &&
           range_min[3] - 1.0e-05 < node.coordinate[3] < range_max[3] + 1.0e-05
            push!(node_ids, node.id)
        end
    end
    # 条件値の設定
    value_set = zeros(Float64, length(values), num_step)
    for istep in 1 : num_step
        for j in 1 : length(values)
            value_set[j, istep] = values[j](istep)
        end
    end
    # 境界条件を作成
    return DirichletBc(node_ids, flag, value_set)
end