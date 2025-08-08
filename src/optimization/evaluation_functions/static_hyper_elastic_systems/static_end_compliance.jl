
#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
mutable struct StaticEndComplianceHyperElastic <: AbstractHyperElasticEvalFunc
    x::Vector{Float64}
    U::Vector{Float64}
    value::Float64
    num_step::Int64
    method::AbstractMethod
    #--------------------------------------------------------
    # 内部コンストラクタ
    #--------------------------------------------------------
    function StaticEndComplianceHyperElastic(
        x0::Vector{Float64}, 
        physics::StaticSolidMechanics,
        num_step::Int64,
        method::AbstractMethod=NewtonRapson()
        )
        return new(
            x0, 
            zeros(Float64, physics.num_total_eq), 
            0.0, 
            num_step, 
            method
        )
    end
end
#-------------------------------------------------------------------------------------------------
# 評価関数
#-------------------------------------------------------------------------------------------------
function compute_eval_func(
    eval::StaticEndComplianceHyperElastic, 
    physics::StaticSolidMechanics, 
    opt_settings::OptimizationSettings,
    opt_filter::F
    )::Float64 where{F<:AbstractOptFilter}

    # 設計変数のフィルタリング
    x_tilde = filtering(opt_filter, eval.x)

    # 材料物性のデータ作成
    element_states = make_element_states(
        x_tilde, 
        physics.elements, 
        opt_settings.num_svalue_types, 
        opt_settings.design_space_type, 
        opt_filter
    )

    # 非線形FEMにより、最終ステップの変位を求める
    compute_nonlinear_fem!(eval, physics, element_states, opt_settings, opt_filter)

    # 表面力ベクトルによる外力ベクトルを作成
    Ft = make_Ft(physics.condition.neumann, physics, eval.num_step)

    # 最終的に得られた変位でエンドコンプライアンスを計算
    eval.value = dot(Ft, eval.U) 
    return eval.value

end
#-------------------------------------------------------------------------------------------------
# 評価関数
#-------------------------------------------------------------------------------------------------
function compute_grad_by_automatic_adjoint(
    eval::StaticEndComplianceHyperElastic, 
    physics::StaticSolidMechanics, 
    opt_settings::OptimizationSettings,
    opt_filter::AbstractOptFilter
    )::Vector{Float64}

    # 設計変数のフィルタリング
    x_tilde = filtering(opt_filter, eval.x)

    # 材料物性のデータ作成
    element_states = make_element_states(
        x_tilde, 
        physics.elements, 
        opt_settings.num_svalue_types, 
        opt_settings.design_space_type, 
        opt_filter
    )

    # 随伴方程式の作成
    f_u(U) = compute_compliance(U, physics, element_states, eval.num_step)
    dRdU, _ = make_Kt_Fint(physics, element_states, eval.U)
    dfdU = - Enzyme.gradient(Enzyme.Reverse, Const(f_u), eval.U)[1]

    # ディレクレ境界の考慮
    lhs, rhs = compute_dirichlet_bc(
        physics.condition.dirichlet, 
        physics.nodes, 
        physics.num_total_eq, 
        dRdU, 
        dfdU, 
        eval.num_step, 
        zero(eval.U)
    )

    # 随伴方程式を解く
    lambda = linear_solve(physics.linear_solver, lhs, rhs)

    # 関数の定義
    L(x) = compute_lagrangian_compliance(x, physics, opt_settings, eval.U, lambda, opt_filter, eval.num_step)
    dfdx_tilde = Enzyme.gradient(Enzyme.set_runtime_activity(Enzyme.Reverse), Const(L), x_tilde)[1]

    # フィルタリングの感度を考慮
    df = filtering_df(opt_filter, dfdx_tilde, eval.x)

    #
    return df
end
#-------------------------------------------------------------------------------------------------
# 評価関数
#-------------------------------------------------------------------------------------------------
function compute_grad_by_adjoint(
    eval::StaticEndComplianceHyperElastic, 
    physics::StaticSolidMechanics, 
    opt_settings::OptimizationSettings,
    opt_filter::AbstractOptFilter
    )::Vector{Float64}

    # 設計変数のフィルタリング
    x_tilde = filtering(opt_filter, eval.x)

    # 材料物性のデータ作成
    element_states = make_element_states(
        x_tilde, 
        physics.elements, 
        opt_settings.num_svalue_types, 
        opt_settings.design_space_type, 
        opt_filter
    )

    # 随伴方程式の左辺
    dRdU, _ = make_Kt_Fint(physics, element_states, eval.U)

    # 随伴方程式の右辺
    dfdU = - make_Ft(physics.condition.neumann, physics, eval.num_step)

    # ディレクレ境界の考慮
    lhs, rhs = compute_dirichlet_bc(
        physics.condition.dirichlet, 
        physics.nodes, 
        physics.num_total_eq, 
        dRdU, 
        dfdU, 
        eval.num_step, 
        zero(eval.U)
    )

    # 随伴方程式を解く
    lambda = linear_solve(physics.linear_solver, lhs, rhs)

    # 感度解析
    dfdx_tilde = zero(x_tilde)
    # 要素ループ
    for (e, element) in enumerate(physics.elements)
        # 要素内自由度リストを作成
        dof_list = element.dof_list
        # 要素変位を作成
        Ue = make_element_sol(dof_list, eval.U)
        lam_e = make_element_sol(dof_list, lambda)
        # 材料ループ
        for itype in 1 : opt_settings.num_svalue_types
            # 位置eにおける全ての設計変数のベクトルを取得
            xe = make_local_design_variables(x_tilde, opt_settings.num_svalue_types, e)
            # 積分点ループ
            temp = compute_element_explicit_sensitivity(xe, itype, element, element_states[e], Ue, lam_e, eval, element.type)
            #
            index = make_index(opt_settings.num_svalue_types, itype, e)
            dfdx_tilde[index] += temp
        end
    end

    # フィルタリングの感度を考慮
    df = filtering_df(opt_filter, dfdx_tilde, eval.x)

    #
    return df
end
#-------------------------------------------------------------------------------------------------
# 評価関数
#-------------------------------------------------------------------------------------------------
function compute_element_explicit_sensitivity(
    xe::Vector{Float64},
    itype::Int64,
    element::StructuralElement,
    element_state::ElementState,
    Ue::Vector{Float64},
    lam_e::Vector{Float64},
    ::StaticEndComplianceHyperElastic,
    ::TotalLagrange 
    )::Float64
    
    # 積分点ループ
    temp = 0.0
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
        F = make_F(element, dNdx, Ue, element.type)
        # Bマトリクスの計算
        B_l, _ = make_BL_and_BNL(element, dNdx, F)
        # 空間勾配の計算
        lam_grad = B_l * lam_e
        # 応力感度の計算
        dS = compute_stress_sensitivity(element.material, material, F, xe, itype)
        # 陽的項の計算
        factor = integral_factor(element.shape, ip) * det(J)
        temp += dot(lam_grad, dS) * factor
    end
    return temp
end
#-------------------------------------------------------------------------------------------------
# 評価関数
#-------------------------------------------------------------------------------------------------
function compute_compliance(
    U::Vector{Float64}, 
    physics::StaticSolidMechanics, 
    ::Vector{ElementState}, 
    istep::Int64
    )::Float64

    # 表面力ベクトルによる外力ベクトルを作成
    Ft = make_Ft(physics.condition.neumann, physics, istep)
    # 関数を計算する
    return dot(Ft, U)
end
#-------------------------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------------------------
function compute_lagrangian_compliance(
    x::Vector{Float64}, 
    physics::StaticSolidMechanics, 
    opt_settings::OptimizationSettings,
    U::Vector{Float64}, 
    lambda::Vector{Float64},
    opt_filter::F,
    istep::Int64
    )::Float64 where{F<:AbstractOptFilter}

    # 材料物性のデータ作成
    element_states = make_element_states(x, physics.elements, opt_settings.num_svalue_types, opt_settings.design_space_type, opt_filter)
    # 表面力ベクトルによる外力ベクトルを作成
    Ft = make_Ft(physics.condition.neumann, physics, istep)
    #
    Fext = Ft
    # 関数の計算
    f = dot(Fext, U)
    R = make_Fint(physics, element_states, U) .- Fext
    lamR = dot(lambda, R)
    #
    return f + lamR
end