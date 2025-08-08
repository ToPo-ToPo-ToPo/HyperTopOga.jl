
#------------------------------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------------------------------
struct NewtonRapson <: AbstractMethod end
#------------------------------------------------------------------------------------------------------
# Newton-Rapson法による収束計算を行う
#------------------------------------------------------------------------------------------------------
function newton_rapson(
    ::NewtonRapson, 
    physics::StaticSolidMechanics, 
    element_states::Vector{ElementState}, 
    istep::Int64, 
    U0::Vector{Float64}, 
    output_flag::Bool=false
    )::Vector{Float64}

    # print
    if output_flag == true
        println("--------------------------------------------------------------------------------------")
        println("Incremental step: $(istep)")
        println("--------------------------------------------------------------------------------------")
    end

    # 解ベクトルの初期化
    U = copy(U0)
    dU = similar(U)

    # 現ステップの外部仮想仕事ベクトルを作成
    Fext = make_Fext(physics, element_states, istep)

    # 全体接線マトリクスと内部仮想仕事ベクトルを作成
    K, Fint = make_Kt_Fint(physics, element_states, U)

    # 初期残差ベクトルとノルムを計算
    R = Fint - Fext

    # 支配方程式を作成
    lhs = K
    rhs = - R

    # ディレクレ境界の考慮
    lhs, rhs = compute_dirichlet_bc(
        physics.condition.dirichlet, 
        physics.nodes, 
        physics.num_total_eq, 
        lhs, 
        rhs, 
        istep, 
        U
    )
    
    # 収束判定値の計算
    eps1 = 1.0e-07
    eps2 = 1.0e-05
    Fext_norm = my_norm(Fext)
    tol = eps1 * Fext_norm + eps2
    R_norm = 1.0e+10
    #-----------------------------------------------
    # 収束計算
    #-----------------------------------------------
    iter_max = 20
    for iter in 1 : iter_max
        
        # 方程式を解き、変位の修正量を求める
        dU .= linear_solve(physics.linear_solver, lhs, rhs)

        # 変位を更新する
        U .+= dU

        # 全体接線マトリクスと内部仮想仕事ベクトルの更新
        K, Fint = make_Kt_Fint(physics, element_states, U)

        # 残差の更新とノルムの計算
        R = Fint - Fext

        # 支配方程式を作成
        lhs .= K
        rhs .= - R

        # ディレクレ境界の考慮
        lhs, rhs = compute_dirichlet_bc(
            physics.condition.dirichlet, 
            physics.nodes, 
            physics.num_total_eq, 
            lhs, 
            rhs, 
            istep, 
            U
        )

        # 残差の更新とノルムの計算
        R_norm = my_norm(rhs)

        # 計算途中の情報を出力
        if output_flag == true
            newton_rapson_printer(iter, Fext_norm, R_norm, dU, U)
        end

        # 収束判定
        if R_norm < tol
            return U
        end
    end

    # 収束したか確認
    println("R_norm = $(R_norm)")
    println("tol = $(tol)")
    error("newton_rapson(): Not converge!")
end
#----------------------------------------------------------------------------
# 外力ベクトルの作成
#----------------------------------------------------------------------------
function make_Fext(
    physics::P, 
    element_states::Vector{ElementState}, 
    istep::Int64
    )::Vector{Float64} where {P<:SolidMechanics}
    #
    Ft = make_Ft(physics.condition.neumann, physics, istep)
    #
    return Ft
end