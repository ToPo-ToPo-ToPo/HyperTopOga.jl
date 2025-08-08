

#----------------------------------------------------------------------------
# 出力用関数 
#----------------------------------------------------------------------------
function newton_rapson_printer(iter::Int64, Fext_norm::Float64, R_norm::Float64, dU::Vector{Float64}, U::Vector{Float64})::Nothing
    dU_norm = my_norm(dU)
    U_norm = my_norm(U)
    @printf("Step: %-2d |R0| = %-.4e, |R| = %-.4e, |dU| = %-.4e, |U| = %-.4e\n", iter, Fext_norm, R_norm, dU_norm, U_norm)
    return nothing
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(
    ::typeof(newton_rapson_printer), 
    iter::Int64, 
    Fext_norm::Float64, 
    R_norm::Float64, 
    dU::Vector{Float64}, 
    U::Vector{Float64}
    )
    return nothing
end