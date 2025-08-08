
#----------------------------------------------------------------
# 設計変数の情報をベクトルとして表現するため必要
#----------------------------------------------------------------
function make_index(num_types::Int64, type::Int64, e::Int64)::Int64
    return type + (e - 1) * num_types
end
#----------------------------------------------------------------
# 位置eにおける材料データのベクトルを取得
#----------------------------------------------------------------
function make_local_design_variables(x::Vector{Float64}, num_types::Int64, e::Int64)::Vector{Float64}
    # 初期化
    xe = zeros(Float64, num_types)
    # ある位置eにおける設計変数の組を設定
    for type in eachindex(xe)
        index = make_index(num_types, type, e)
        xe[type] = x[index]
    end
    #
    return xe
end