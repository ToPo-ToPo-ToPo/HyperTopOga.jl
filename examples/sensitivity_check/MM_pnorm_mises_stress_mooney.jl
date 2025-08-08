
using Revise
using Random
using HyperTopOga
#-------------------------------------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------------------------------------
function main()
    
    #--------------------------------------------------------------------------------------------
    # 解析モデルの設定
    length_x = 200.0
    length_y = 100.0
    division_x = 10
    division_y = 5
    thickness = 1.0

    # 材料モデルの定義 固体領域
    young1 = 5.0
    poisson1 = 0.3
    # lame定数
    K1 = (poisson1 * young1) / ((1.0+poisson1) * (1.0-2.0*poisson1))
    mu1 = 0.5 * young1 / (1.0 + poisson1)
    #
    young2 = 10.0
    poisson2 = 0.3
    # lame定数
    K2 = (poisson2 * young1) / ((1.0+poisson2) * (1.0-2.0*poisson2))
    mu2 = 0.5 * young2 / (1.0 + poisson2)
    #
    C1 = HyperTopOga.DMOSIMP([1.0e-03, 0.5*mu1, 0.5*mu2], 3.0)
    C2 = HyperTopOga.DMOSIMP([1.0e-03, 0.2*mu1, 0.2*mu2], 3.0)
    K = HyperTopOga.DMOSIMP([1.0e-03, 2.0*K1, 2.0*K2], 3.0)
    density = HyperTopOga.Constant(1.0e-05)
    gamma_x = HyperTopOga.Constant(1.0)
    sigma_limit = HyperTopOga.DMOSIMP([1.0e-03, 2.0, 10.0], 2.5)
    material_data = HyperTopOga.MooneyRivlinPlaneStrainData(
        C1, C2, HyperTopOga.Constant(0.0), HyperTopOga.Constant(0.0), K, density, gamma_x
    )
    interpolate = HyperTopOga.DMOSIMP([0.0, 1.0, 2.0], 1.0)

    # 検証内容
    num_step = 10
    linear_solver = HyperTopOga.DirectSolverUMFPACK()
    output_file_name = "output/grad_check"
    output_flag = true
    #
    problem_type = HyperTopOga.TotalLagrange()
    #problem_type = HyperTopOga.StabilisedTotalLagrange()

    #--------------------------------------------------------------------------------------------
    # 最適化の設定
    num_svalue_types = 2
    p = 16.0
    #
    opt_filter = HyperTopOga.HeavisideProjectionFilter()
    radius = 3.0

    # 解析モデルの作成
    nodes, connects = HyperTopOga.make_quad4_voxel_mesh(
        length_x, length_y, division_x, division_y
    )

    # 要素情報の作成
    elements = HyperTopOga.make_structural_ele_quad4(
        connects, material_data, problem_type, thickness
    )

    # 境界条件・ソース項の情報を設定
    bc = HyperTopOga.ConditionPack()

    # ディレクレ条件
    range_min = [0.0, 0.0, 0.0]
    range_max = [0.0, 100.0, 0.0]
    dirichlet_flag1 = [true, true]
    u_d1_x(x) = 0.0
    u_d1_y(x) = 0.0
    dirichlet_values1 = [u_d1_x, u_d1_y]
    HyperTopOga.add_condition(
        bc, 
        HyperTopOga.make_dirichlet_voxel_model(
            nodes, 
            range_min, 
            range_max,
            num_step, 
            dirichlet_flag1, 
            dirichlet_values1
        )
    )

    # ノイマン条件
    range_min = [200.0, 0.0, 0.0]
    range_max = [200.0, 100.0, 0.0]
    traction_x(x) = 0.0
    traction_y(x) = -0.5e-01 * (x / Float64(num_step))
    neumann_values = [traction_x, traction_y]
    HyperTopOga.add_condition(
        bc, 
        HyperTopOga.make_neumann_voxel_model(
            nodes, 
            elements,
            range_min, 
            range_max,
            num_step,  
            neumann_values
        )
    )

    # モデルの作成
    physics = HyperTopOga.StaticSolidMechanics(
        nodes, 
        elements, 
        bc, 
        linear_solver, 
        output_file_name
    )

    #--------------------------------------------------------------------------------------------
    # Setting of topology optimization
    opt_settings = HyperTopOga.OptimizationSettings(num_svalue_types, interpolate, output_flag)

    # 最適化用のフィルターのオブジェクトの作成と設定
    HyperTopOga.filter_setup!(opt_filter, num_svalue_types, length(elements), physics, radius)

    # 微係数を求めたい値
    Random.seed!(10) # 同じシード値で設計変数を生成したい場合は必要
    x = rand(length(elements)*num_svalue_types)
    #x = fill(0.4, length(elements))

    # 評価関数の定義を更新
    eval = HyperTopOga.StaticMisesStressPnormHyperElastic(x, physics, sigma_limit, p, num_step)
    
    #--------------------------------------------------------------------------------------------
    HyperTopOga.check_gradient(
        x, eval, physics, opt_settings, opt_filter, 
        HyperTopOga.AutomaticDiff(), 
        HyperTopOga.AutomaticAdjointDiff(), 
        HyperTopOga.NumericalDiff()
    )

end
#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
main()
