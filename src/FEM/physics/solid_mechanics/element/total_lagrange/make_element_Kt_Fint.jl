
#---------------------------------------------------------------------------------------------
# 要素接線剛性マトリクスの作成
# トータルラグランジュに基づく有限変形用
#---------------------------------------------------------------------------------------------
function make_element_Kt_Fint(
    element::E, 
    element_state::ElementState, 
    Ue::Vector{Float64}, 
    problem_type::TotalLagrange
)::Tuple{Matrix{Float64}, Vector{Float64}} where {E<:StructuralElement}

    # 要素自由度
    num_dof = element.num_element_dof
    # 初期化
    Ke = zeros(Float64, num_dof, num_dof)
    Fe = zeros(Float64, num_dof)
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
        B_l, B_nl = make_BL_and_BNL(element, dNdx, F)
        # 固体の構成則の計算
        S, S_bar, Ct = compute_constitutive_law(material, F)
        #
        factor = integral_factor(element.shape, ip) * my_det(J)
        # 要素内力ベクトルの計算
        Fe .+= factor .* (B_l' * S)
        # 要素剛性行列の計算
        Ke .+= factor .* (B_l' * Ct * B_l .+ B_nl' * S_bar * B_nl)
    end
    #
    return Ke, Fe
end