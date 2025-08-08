
#-------------------------------------------------------------------------------
# 内部仮想仕事ベクトルの作成
#-------------------------------------------------------------------------------
function make_Fint(physics::P, element_states::Vector{ElementState}, U::Vector{Float64}
    )::Vector{Float64} where {P<:SolidMechanics}

    # zero clear
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
            Fe = make_element_Fint(element, element_states[id], Ue, element.type)
            # アセンブリング
            F[dof_list] += Fe
        end
    end
    # 
    return F
end