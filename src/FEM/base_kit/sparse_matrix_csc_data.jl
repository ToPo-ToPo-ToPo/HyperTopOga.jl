#---------------------------------------------------------------------------------------
# 疎行列データをまとめた構造体
#---------------------------------------------------------------------------------------
struct SparseMatrixCSCData
    nrows::Int64
    ncols::Int64
    colptr::Vector{Int64}                        # 列インデックス
    rowval::Vector{Int64}                        # 行ポインタ
    index_cache::Dict{Int64, Dict{Int64, Int64}} # インデックスの探索情報
    nnz::Int64
end
#----------------------------------------------------------------------------
# 疎行列を初期化
#----------------------------------------------------------------------------
function my_spzeros(sp_data::SparseMatrixCSCData)::SparseMatrixCSC{Float64, Int64}
    return SparseMatrixCSC(
        sp_data.nrows, 
        sp_data.ncols, 
        sp_data.colptr, 
        sp_data.rowval, 
        zeros(Float64, sp_data.nnz)
        )
end
#----------------------------------------------------------------------------
# 疎行列を初期化
#----------------------------------------------------------------------------
function my_sparse(sp_data::SparseMatrixCSCData, nzval::Vector{Float64})::SparseMatrixCSC{Float64, Int64}
    return SparseMatrixCSC(
        sp_data.nrows, 
        sp_data.ncols, 
        sp_data.colptr, 
        sp_data.rowval, 
        nzval
        )
end
#---------------------------------------------------------------------------------------
# FEMにおける疎行列のデータ構造をあらかじめ作成
#---------------------------------------------------------------------------------------
function create_sparse_matrix_csc_data(nodes::Vector, elements::Vector, num_total_eq::Int64)::SparseMatrixCSCData

    # create node-node table
    node2nodes = Vector{Data2Data}(undef, length(nodes))
    for i in 1 : length(nodes)
        node2nodes[i] = Data2Data()
    end
    for element in elements
        for node_id_i in element.shape.connects
            for node_id_j in element.shape.connects
                # 節点ID_iに関連する節点番号ID_jを追加する
                push!(node2nodes[node_id_i].data, node_id_j)
            end
        end
    end

    # 重複する節点番号の削除とソートを行う
    node2nodes = delete_duplicate_and_sort(node2nodes)

    # create dof-dof table
    # 節点に関連する自由度のリストを作成
    dof2dofs = Vector{Data2Data}(undef, num_total_eq)
    for i in 1 : num_total_eq
        dof2dofs[i] = Data2Data()
    end

    for (i, node2node) in enumerate(node2nodes)
        for idof in nodes[i].dof_ids
            for ID in node2node.data
                for jdof in nodes[ID].dof_ids
                    push!(dof2dofs[idof].data, jdof)
                end
            end
        end
    end

    # 重複する自由度の削除と自由度番号のソートを行う
    dof2dofs = delete_duplicate_and_sort(dof2dofs)

    # 非ゼロ要素の数を取得
    nnz = 0
    for dof2dof in dof2dofs
        nnz += length(dof2dof.data)
    end

    # indexベクトルの作成
    colptr = Vector{Int64}(undef, num_total_eq+1)
    colptr[1] = 1
    for i in 2 : length(colptr)
        colptr[i] = colptr[i-1] + length(dof2dofs[i-1].data)
    end

    # columnベクトルの作成
    rowval = Vector{Int64}(undef, nnz)
    counter = 0
    for dof2dof in dof2dofs
        for idof in dof2dof.data
            counter += 1
            rowval[counter] = idof
        end
    end

    # インデックスのキャッシュを作成
    index_cache = Dict{Int64, Dict{Int64, Int64}}()
    for col in 1 : num_total_eq
        
        # 非ゼロの位置に関する情報
        row_start = colptr[col]
        row_end = colptr[col+1] - 1
        row_indices = rowval[row_start:row_end]

        # この列の辞書を作成
        col_cache = Dict{Int64, Int64}()
        for (idx, row) in enumerate(row_indices)
            col_cache[row] = row_start + idx - 1
        end
        index_cache[col] = col_cache
    end

    # 疎行列データを作成して返す
    return SparseMatrixCSCData(num_total_eq, num_total_eq, colptr, rowval, index_cache, nnz)
end