

#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
mutable struct StaticMisesCauchyStressPnormHyperElastic <: AbstractHyperElasticEvalFunc
    x::Vector{Float64}
    U::Vector{Float64}
    value::Float64
    num_step::Int64
    interpolate_sigma_limit::AbstractInterpolation
    p::Float64
    method::AbstractMethod
    #--------------------------------------------------------
    # 内部コンストラクタ
    #--------------------------------------------------------
    function StaticMisesCauchyStressPnormHyperElastic(
        x0::Vector{Float64}, 
        physics::StaticSolidMechanics, 
        interpolate_sigma_limit::AbstractInterpolation,
        p::Float64, 
        num_step::Int64,
        method::AbstractMethod=NewtonRapson()
        )
        return new(
            x0, 
            zeros(Float64, physics.num_total_eq), 
            0.0,
            num_step, 
            interpolate_sigma_limit, 
            p,
            method
        )     
    end
end
#----------------------------------------------------------------------------
# 評価関数
#----------------------------------------------------------------------------
function compute_eval_func(
    eval::StaticMisesCauchyStressPnormHyperElastic, 
    physics::StaticSolidMechanics, 
    opt_settings::OptimizationSettings,
    opt_filter::AbstractOptFilter
    )::Float64

    # 設計変数のフィルタリング
    x_tilde = filtering(opt_filter, eval.x)

    # 材料物性のデータ作成
    materials = make_material_model(
        x_tilde, 
        physics.elements, 
        opt_settings.num_svalue_types, 
        opt_settings.design_space_type, 
        opt_filter
    )

    # 非線形FEMにより、最終ステップの変位を求める
    compute_nonlinear_fem!(eval, physics, materials, opt_settings, opt_filter)

    # 許容値の計算
    sigma_limit = make_sigma_limit(
        x_tilde, 
        eval.interpolate_sigma_limit, 
        opt_settings.num_svalue_types, 
        length(physics.elements)
    )

    #
    eval.value = static_mises_cauchy_stress_pnorm(physics, materials, sigma_limit, eval.p, eval.U)
    return eval.value
end
#----------------------------------------------------------------------------
# 随伴法に基づく感度解析
#----------------------------------------------------------------------------
function compute_grad_by_automatic_adjoint(
    eval::StaticMisesCauchyStressPnormHyperElastic, 
    physics::StaticSolidMechanics, 
    opt_settings::OptimizationSettings,
    opt_filter::AbstractOptFilter
    )

    # 設計変数のフィルタリング
    x_tilde = filtering(opt_filter, eval.x)

    # 材料物性のデータ作成
    materials = make_material_model(
        x_tilde, 
        physics.elements, 
        opt_settings.num_svalue_types, 
        opt_settings.design_space_type, 
        opt_filter
    )

    # 許容値の計算
    sigma_limit = make_sigma_limit(
        x_tilde, 
        eval.interpolate_sigma_limit, 
        opt_settings.num_svalue_types,
        length(physics.elements)
    )

    # 随伴方程式の作成
    f_u(U) = compute_mises_cauchy_stress_pnorm(U, physics, materials, sigma_limit, eval.p)
    dRdU, _ = make_Kt_Fint(physics, materials, eval.U)
    dfdU = - gradient(Reverse, Const(f_u), eval.U)

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
    lambda = linear_solve(physics.linear_solver_type, lhs, rhs)

    # 関数の定義
    L(x) = compute_lagrangian_mises_cauchy_stress_pnorm(
        x, 
        physics, 
        eval.interpolate_sigma_limit, 
        eval.p, 
        eval.num_step, 
        opt_settings.num_svalue_types, 
        eval.U, 
        lambda,
        opt_filter, 
        opt_settings.design_space_type
    )
    dfdx_tilde = gradient(Enzyme.Reverse, Const(L), x_tilde)

    # フィルタリングの感度を考慮
    df = filtering_df(opt_filter, dfdx_tilde, eval.x)

    #
    return df
end
#----------------------------------------------------------------------------
# 評価関数
#----------------------------------------------------------------------------
function compute_mises_cauchy_stress_pnorm(
    U::Vector{Float64}, 
    physics::StaticSolidMechanics, 
    materials::Vector{AbstractMaterial}, 
    sigma_limit::Vector{Float64}, 
    p::Float64
    )

    return static_mises_cauchy_stress_pnorm(physics, materials, sigma_limit, p, U)
end
#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
function compute_lagrangian_mises_cauchy_stress_pnorm(
    x_tilde::Vector{Float64}, 
    physics::SolidMechanics, 
    interpolate_sigma_limit::AbstractInterpolation, 
    p::Float64, 
    num_step::Int64, 
    num_svalue_types::Int64, 
    U::Vector{Float64}, 
    lambda::Vector{Float64},
    opt_filter::AbstractOptFilter,
    design_space_type::ElementBase
    )

    # 材料物性のデータ作成
    materials = make_material_model(
        x_tilde, 
        physics.elements, 
        num_svalue_types, 
        design_space_type, 
        opt_filter
    )

    # 許容応力のデータ作成
    sigma_limit = make_sigma_limit(
        x_tilde, 
        interpolate_sigma_limit, 
        num_svalue_types, 
        length(physics.elements)
    )

    # 表面力ベクトルによる外力ベクトルを作成
    Ft = make_Ft(physics.condition.neumann, physics, num_step)
    
    #
    Fext = Ft
    R = make_Fint(physics, materials, U) .- Fext

    # 関数の計算
    f = static_mises_cauchy_stress_pnorm(physics, materials, sigma_limit, p, U)
    lamR = dot(lambda, R)

    #
    return f + lamR
    
end
#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
function static_mises_cauchy_stress_pnorm(
    physics::StaticSolidMechanics, 
    materials::Vector{AbstractMaterial}, 
    sigma_limit::Vector{Float64}, 
    p::Float64, 
    U::Vector{Float64} 
    )

    # 要素ループ
    value = 0.0
    for (e, element) in enumerate(physics.elements)
        # 要素変位を作成
        Ue = make_element_sol(element.dof_list, U)
        # 積分点ループ
        for ip in 1 : element.shape.num_eval_points
            mises = compute_element_mises_cauchy_stress(
                element, 
                materials[e], 
                ip, 
                Ue, 
                physics.problem_type
            )
            value += (mises / sigma_limit[e])^p
        end
    end

    #
    return value^(1.0 / p)
end
#----------------------------------------------------------------------------
# 評価関数
#----------------------------------------------------------------------------
function compute_element_mises_cauchy_stress(
    element::AbstractElement, 
    material::AbstractMaterial, 
    ip::Int64, 
    Ue::Vector{Float64}, 
    problem_type::TotalLagrange
    )
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
    # 固体の構成則の計算
    S, F = compute_stress(material, F, ip)
    # Push forward
    sigma = (F * S * F') ./ my_det(F)
    # ミーゼス応力の計算
    mises = compute_mises_stress(sigma)
    #
    return mises
end