
#----------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------
mutable struct StaticMisesStressPnormHyperElastic <: AbstractHyperElasticEvalFunc
    x::Vector{Float64}
    U::Vector{Float64}
    value::Float64
    num_step::Int64
    sigma_limit::AbstractInterpolation
    p::Float64
    method::AbstractMethod
    #--------------------------------------------------------
    # 内部コンストラクタ
    #--------------------------------------------------------
    function StaticMisesStressPnormHyperElastic(
        x0::Vector{Float64}, 
        physics::StaticSolidMechanics, 
        sigma_limit::AbstractInterpolation, 
        p::Float64, 
        num_step::Int64,
        method::AbstractMethod=NewtonRapson()
        )
        return new(
            x0, zeros(Float64, physics.num_total_eq), 
            0.0, num_step, 
            sigma_limit, p,
            method
        )   
    end
end
#----------------------------------------------------------------------------
# 評価関数
#----------------------------------------------------------------------------
function compute_eval_func(
    eval::StaticMisesStressPnormHyperElastic, 
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

    #
    eval.value = compute_mises_stress_pnorm(
        x_tilde,
        opt_settings.num_svalue_types,
        physics, 
        element_states, 
        eval.sigma_limit, 
        eval.p, 
        eval.U
    )

    #
    return eval.value
end
#----------------------------------------------------------------------------
# 随伴法に基づく感度解析
#----------------------------------------------------------------------------
function compute_grad_by_automatic_adjoint(
    eval::StaticMisesStressPnormHyperElastic, 
    physics::StaticSolidMechanics, 
    opt_settings::OptimizationSettings,
    opt_filter::F
)::Vector{Float64} where{F<:AbstractOptFilter}

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
    f_u(U) = compute_mises_stress_pnorm(
        U,
        x_tilde, opt_settings.num_svalue_types,
        physics, 
        element_states, 
        eval.sigma_limit, 
        eval.p
    )
    dRdU, _ = make_Kt_Fint(physics, element_states, eval.U)
    dfdU = - gradient(Enzyme.set_runtime_activity(Reverse), Const(f_u), eval.U)[1]

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
    L(x) = compute_lagrangian_mises_stress_pnorm(
        x, 
        physics, 
        eval.sigma_limit, eval.p, 
        eval.num_step, 
        opt_settings.num_svalue_types, 
        eval.U, 
        lambda,
        opt_filter,
        opt_settings.design_space_type
    )
    dfdx_tilde = gradient(Enzyme.set_runtime_activity(Enzyme.Reverse), Const(L), x_tilde)[1]

    # フィルタリングの感度を考慮
    df = filtering_df(opt_filter, dfdx_tilde, eval.x)

    #
    return df
end
#----------------------------------------------------------------------------
# 評価関数
#----------------------------------------------------------------------------
function compute_mises_stress(
    element::E, 
    material::M, 
    ip::Int64, 
    Ue::Vector{Float64}, 
    problem_type::TotalLagrange
)::Float64 where {E<:StructuralElement, M<:AbstractMaterial}

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
    S, _ = compute_stress(material, F)
    # ミーゼス応力の計算
    mises = compute_mises_stress(S)
    #
    return mises
end
#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
function compute_mises_stress_pnorm(
    x_tilde::Vector{Float64},
    num_svalue_types::Int64,
    physics::StaticSolidMechanics, 
    element_states::Vector{ElementState},
    sigma_limit::AbstractInterpolation,
    p::Float64, 
    U::Vector{Float64} 
)::Float64

    # 要素ループ
    value = 0.0
    for (e, element) in enumerate(physics.elements)
        #
        xe = make_local_design_variables(x_tilde, num_svalue_types, e)
        # 要素変位を作成
        Ue = make_element_sol(element.dof_list, U)
        # 積分点ループ
        for (ip, material) in enumerate(element_states[e].materials)
            # ミーゼス応力の計算
            mises = compute_mises_stress(element, material, ip, Ue, element.type)
            # 許容応力の計算
            sigma_limit_local = compute(sigma_limit, xe)
            #
            value += (mises / sigma_limit_local)^p
        end
    end

    #
    return value^(1.0 / p)
end
#----------------------------------------------------------------------------
# 評価関数
#----------------------------------------------------------------------------
function compute_mises_stress_pnorm(
    U::Vector{Float64},
    x_tilde::Vector{Float64},
    num_svalue_types::Int64, 
    physics::SolidMechanics, 
    element_states::Vector{ElementState}, 
    sigma_limit::AbstractInterpolation,
    p::Float64
    )::Float64
    return compute_mises_stress_pnorm(
        x_tilde, 
        num_svalue_types, 
        physics, 
        element_states,
        sigma_limit, 
        p, 
        U
    )
end
#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
function compute_lagrangian_mises_stress_pnorm(
    x_tilde::Vector{Float64}, 
    physics::SolidMechanics, 
    sigma_limit::AbstractInterpolation, 
    p::Float64, 
    num_step::Int64, 
    num_svalue_types::Int64, 
    U::Vector{Float64}, 
    lambda::Vector{Float64},
    opt_filter::AbstractOptFilter, 
    design_space_type::ElementBase
    )::Float64

    # 材料物性のデータ作成
    element_states = make_element_states(
        x_tilde, 
        physics.elements, 
        num_svalue_types, 
        design_space_type, 
        opt_filter
    )

    # 表面力ベクトルによる外力ベクトルを作成
    Ft = make_Ft(physics.condition.neumann, physics, num_step)
    
    #
    Fext = Ft
    R = make_Fint(physics, element_states, U) .- Fext

    # 関数の計算
    f = compute_mises_stress_pnorm(
        x_tilde, num_svalue_types, physics, element_states, sigma_limit, p, U
    )
    lamR = dot(lambda, R)

    #
    return f + lamR
    
end
