
#----------------------------------------------------------------
# 要素変位ベクトルの作成
#----------------------------------------------------------------
function make_element_sol(dof_list::Vector{Int64}, sol::Vector{Float64})::Vector{Float64}
    return sol[dof_list]
end