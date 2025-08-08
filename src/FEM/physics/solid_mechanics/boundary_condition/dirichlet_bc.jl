
#----------------------------------------------------------------------------
#　境界条件
#----------------------------------------------------------------------------
struct DirichletBc <: AbstractDirichletBc
    node_ids::Vector{Int64}                # 条件を与える節点番号
    flag::Vector{Bool}                     # 条件を与えるかどうかのフラグ
    values::Matrix{Float64}                # 条件値
end
#----------------------------------------------------------------
# 境界条件から自由度番号を付け直す
#----------------------------------------------------------------
function count_equation(conditions::Vector{DirichletBc}, nodes::Vector{Node})::Float64

    # 総自由度数の計算
    num_total_eq = 0
    for node in nodes
        num_total_eq += length(node.dof_ids)
    end
    
    # ディレクレ条件に含まれる節点の状態場は未知数でない
    # 該当する節点に印をつける
    for condition in conditions
        # ディレクレ境界条件を設定された節点ループ 
        for node_id in condition.node_ids
            # 自由度ループを行い、自由度番号を借り設定
            for i in 1 : length(nodes[node_id].dof_ids)
                # フラグが設定されているか確認し、されていない場合は、処置を行わない
                if condition.flag[i] == false
                    continue
                end
                nodes[node_id].dof_ids[i] = - 1
            end
        end
    end

    # 自由度番号の再定義
    eq_counter = 0
    for node_id in 1 : length(nodes)
        for idof in 1 : length(nodes[node_id].dof_ids)
            # ディレクレ境界に含まれる場合は、除外
            if nodes[node_id].dof_ids[idof] == - 1
                continue
            end
            # カウンターの更新
            eq_counter += 1
            # 自由度番号をつける
            nodes[node_id].dof_ids[idof] = eq_counter
        end
    end

    # 未知数の数
    num_free_eq = eq_counter

    # ディレクレ条件に設定に設定された節点に自由度番号を割り当てる
    spc_counter = eq_counter
    for condition in conditions
        # ディレクレ境界条件を設定された節点ループ 
        for node_id in condition.node_ids
            # 自由度ループを行い、自由度番号を設定
            for idof in 1 : length(nodes[node_id].dof_ids)
                if nodes[node_id].dof_ids[idof] == -1
                    spc_counter += 1
                    nodes[node_id].dof_ids[idof] = spc_counter
                end 
            end
        end
    end

    # 結果の確認
    if num_total_eq != spc_counter
        println("num_total_eq: ", num_total_eq)
        println("spc_counter: ", spc_counter)
        error("The number of total degrees of freedom is incorrect.")
    end

    return num_total_eq
end
#----------------------------------------------------------------
# 境界条件の考慮
#----------------------------------------------------------------
function compute_dirichlet_bc(
    bcs::Vector{DirichletBc}, 
    nodes::Vector{Node}, 
    num_total_eq::Int64, 
    lhs_bc::SparseMatrixCSC{Float64, Int64}, 
    rhs_bc::Vector{Float64}, 
    istep::Int64,
    U::Vector{Float64}=zeros(Float64, num_total_eq)
    )::Tuple{SparseMatrixCSC{Float64, Int64}, Vector{Float64}}

    # Early check
    if length(bcs) == 0
        return lhs_bc, rhs_bc
    end

    # ディレクレ境界条件の適用
    for bc in bcs
        # 条件が与えられる節点ループ
        for node_id in bc.node_ids
            # 自由度番号のループ
            for (i, idof) in enumerate(nodes[node_id].dof_ids)
                # フラグが立っているか確認 立っていない場合は処理なし
                if bc.flag[i] == false
                    continue
                end
                # 条件値を取得
                value_dof = bc.values[i, istep] - U[idof]
                if value_dof != 0.0
                    for jdof = 1 : num_total_eq
                        rhs_bc[jdof] -= value_dof * lhs_bc[jdof, idof]
                    end
                end
                # 左辺係数行列の修正
                lhs_bc[:,    idof] = zeros(Float64, num_total_eq)
                lhs_bc[idof,    :] = zeros(Float64, num_total_eq)
                lhs_bc[idof, idof] = 1.0
                # 右辺ベクトルの修正
                rhs_bc[idof] = value_dof
            end
        end
    end
    #
    return lhs_bc, rhs_bc
end