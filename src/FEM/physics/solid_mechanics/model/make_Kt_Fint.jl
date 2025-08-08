
#-------------------------------------------------------------------------------
# 接線と内力を同時に作成
#-------------------------------------------------------------------------------
function make_Kt_Fint(
    physics::StaticSolidMechanics, 
    element_states::Vector{ElementState}, 
    U::Vector{Float64}
    )::Tuple{SparseMatrixCSC{Float64, Int64}, Vector{Float64}}

    # Zero clear
    K_nzval = zeros(Float64, physics.sparse_data.nnz)
    F = zeros(Float64, physics.num_total_eq)
    
    # 要素ループ
    for color in physics.multi_color
        Threads.@threads for id in color.element_ids
            # Get element
            element = physics.elements[id]
            # 要素内自由度リストを作成
            dof_list = element.dof_list
            # 要素変位を作成
            Ue = make_element_sol(dof_list, U)
            # 要素剛性マトリクスを計算する
            Ke, Fe = make_element_Kt_Fint(element, element_states[id], Ue, element.type)
            # アセンブリング
            for (j, jdof) in enumerate(dof_list)
                #
                F[jdof] += Fe[j]
                # キャッシュ取得
                row_cache = physics.sparse_data.index_cache[jdof]  
                #
                for (i, idof) in enumerate(dof_list)
                    if haskey(row_cache, idof)
                        idx = row_cache[idof]
                        K_nzval[idx] += Ke[i, j]
                    end
                end
            end
        end
    end
    # 
    return my_sparse(physics.sparse_data, K_nzval), F
end