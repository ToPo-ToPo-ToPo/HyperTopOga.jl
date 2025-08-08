

#----------------------------------------------------------------------------
# 荷重境界を作成する
#----------------------------------------------------------------------------
function make_Ft(bcs::TV, physics::P, istep::Int64)::Vector{Float64} where {TV<:Vector{AbstractNeumannBc}, P<:SolidMechanics}
    # 初期化
    F = zeros(Float64, physics.num_total_eq)
    # 境界条件のループ
    for bc in bcs
        # 条件を設定
        assemble_Ft!(bc, physics, istep, F)
    end
    #
    return F
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(make_Ft), bcs::TV, physics::P, istep::Int64) where {TV<:Vector{AbstractNeumannBc}, P<:SolidMechanics}
    return nothing
end