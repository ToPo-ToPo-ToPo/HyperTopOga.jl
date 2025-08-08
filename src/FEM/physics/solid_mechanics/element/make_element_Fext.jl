
#----------------------------------------------------------------
# 要素ベクトルを求めるための線積分
#----------------------------------------------------------------
function make_element_Fext(element::E, boundary_id::Int64, value::Vector{Float64}
    )::Vector{Float64} where {E<:StructuralElement}
    
    # zero clear
    Fe = zeros(Float64, element.num_element_dof)
    # 線積分の計算 積分点のループ
    for ip in 1 : element.shape.num_boundary_eval_points
        # 境界におけるipの積分点座標を取得
        eval_coordinate = element.shape.boundary_eval_coordinaites[boundary_id, ip, :]
        # 形状関数の計算
        N = make_N(element, eval_coordinate)
        # 形状関数の微分
        dNdxi = make_dNdxi(element.shape, eval_coordinate)
        # Jacobi行列の計算
        J = make_J(element.shape, dNdxi)
        # Jacobianの計算
        detJ = make_boundary_detJ(element.shape, boundary_id, J)
        # 要素ベクトルの計算
        Fe .+= (boundary_integral_factor(element.shape, boundary_id, ip) * detJ) .* (N' * value)
    end
    return Fe
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(
    ::typeof(make_element_Fext), 
    element::E, 
    boundary_id::Int64, 
    value::Vector{Float64}
    ) where {E<:StructuralElement}
    return nothing
end