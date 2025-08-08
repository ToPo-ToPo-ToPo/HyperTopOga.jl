

#----------------------------------------------------------------------------
#　境界条件
#----------------------------------------------------------------------------
mutable struct NeumannBc <: AbstractNeumannBc
    element_ids::Vector{Int64}             # 条件を与える要素番号の配列
    boundary_ids::Vector{Int64}            # 条件を与える要素内面番号の配列
    values::Matrix{Float64}                # 条件値
end
#----------------------------------------------------------------------------
# 荷重境界を作成する
#----------------------------------------------------------------------------
function assemble_Ft!(bc::NeumannBc, physics::P, step::Int64, F::Vector{Float64})::Nothing where {P<:SolidMechanics}
    # 荷重値を取得
    value = bc.values[:, step]
    # 要素ループ
    for (e, element_id) in enumerate(bc.element_ids)
        # 条件を与える要素
        element = physics.elements[element_id]
        # 条件を与える要素内面番号
        boundary_id = bc.boundary_ids[e]
        # 要素ベクトルの作成
        Fe = make_element_Fext(element, boundary_id, value)
        # アセンブリング
        F[element.dof_list] += Fe
    end
    #
    return nothing
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(
    ::typeof(assemble_Ft!), 
    bc::NeumannBc, 
    physics::P,
    step::Int64, 
    F::Vector{Float64}
    ) where {P<:SolidMechanics}
    return nothing
end