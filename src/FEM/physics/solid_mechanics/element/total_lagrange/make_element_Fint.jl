
#---------------------------------------------------------------------------------------------
# トータルラグランジュに基づく有限変形用
#---------------------------------------------------------------------------------------------
function make_element_Fint(
    element::E, 
    element_state::ElementState, 
    Ue::Vector{Float64}, 
    problem_type::TotalLagrange
)::Vector{Float64} where {E<:StructuralElement}

    # 初期化
    Fe = zeros(Float64, element.num_element_dof)
    # 積分点ループ
    for (ip, material) in enumerate(element_state.materials)
        # 積分点座標を取得
        eval_coordinate = element.shape.eval_coordinates[ip, :]
        # 形状関数の微分
        dNdxi = make_dNdxi(element.shape, eval_coordinate)
        # Jacobi行列の計算
        J = make_J(element.shape, dNdxi)
        # 
        dNdx = make_dNdx(element.shape, dNdxi, J)
        # 変形勾配テンソルF
        F = make_F(element, dNdx, Ue, problem_type)
        # BマトリクスとB_NLマトリクスの計算
        B_l, _ = make_BL_and_BNL(element, dNdx, F)
        # 固体の構成則の計算
        S, _, _ = compute_constitutive_law(material, F)
        #
        factor = integral_factor(element.shape, ip) * my_det(J)
        # 要素内力ベクトルの計算
        Fe .+= factor .* (B_l'*S)
    end
    #
    return Fe
end