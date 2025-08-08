
#--------------------------------------------------------------------------------------------------------------
# 要素の状態の生成(ベースプログラム)
#--------------------------------------------------------------------------------------------------------------
function make_element_states(
    x::Vector{Float64}, 
    elements::Vector{E}, 
    num_svalue_types::Int64, 
    design_space::D, 
    filter::F
    )::Vector{ElementState} where {E<:AbstractElement, D<:AbstractDesignSpace, F<:AbstractOptFilter}
    
    return make_element_states_sub(x, elements, num_svalue_types, design_space, filter)
end
#--------------------------------------------------------------------------------------------------------------
# 要素の状態の生成
#--------------------------------------------------------------------------------------------------------------
function make_element_states_sub(x::Vector{Float64}, elements::Vector{E}, num_svalue_types::Int64, ::ElementBase, ::F
    )::Vector{ElementState} where {E<:AbstractElement, F<:AbstractOptFilter}

    # 要素の状態情報の作成
    element_states = Vector{ElementState}(undef, length(elements))
    
    # 要素ループ
    for (e, element) in enumerate(elements)
        # zerosで初期化した方が、自動微分の計算が安定化するため使用
        xe = zeros(Float64, num_svalue_types, element.shape.num_eval_points)
        # 位置eにおける全ての設計変数のベクトルを取得
        for ip in 1 : element.shape.num_eval_points
            local_x = make_local_design_variables(x, num_svalue_types, e)
            xe[:, ip] = local_x
        end
        # 要素の状態の更新
        element_states[e] = setup(element.material, xe)
    end
    
    #
    return element_states
end