
using Revise
using Random
using HyperTopOga
#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
function main()

    #----------------------------------------------------------
    # 解析モデルの設定
    length_x = 200.0
    length_y = 100.0
    division_x = 80
    division_y = 40 
    thickness = 1.0

    # 材料モデルの定義 固体領域
    interpolate = HyperTopOga.DMOSIMP([0.0, 1.0, 2.0], 1.0)
    young = HyperTopOga.DMOSIMP([1.0e-03, 5.0, 7.0], 3.0)
    poisson = HyperTopOga.Constant(0.3)
    density = HyperTopOga.Constant(1.0e-05)
    gamma_x = HyperTopOga.Constant(1.0)
    material_data = HyperTopOga.NeoHookeanPlaneStrainData(young, poisson, density, gamma_x)

    # 検証内容
    num_step = 10
    linear_solver = HyperTopOga.DirectSolverUMFPACK()
    output_file_name = "output/grad_check"
    output_flag = false
    #
    problem_type = HyperTopOga.TotalLagrange()

    #----------------------------------------------------------
    # 最適化の設定
    num_svalue_types = 2
    design_space_type = HyperTopOga.ElementBase()
    opt_filter = HyperTopOga.HeavisideProjectionFilter()
    radius = 3.0

    #----------------------------------------------------------
    # 解析モデルの作成
    nodes, connects = HyperTopOga.make_quad4_voxel_mesh(length_x, length_y, division_x, division_y)

    # 要素情報の作成
    elements = HyperTopOga.make_structural_ele_quad4(connects, material_data, problem_type, thickness)

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

    #----------------------------------------------------------
    # Setting of topology optimization
    opt_settings = HyperTopOga.OptimizationSettings(
        num_svalue_types, 
        interpolate, 
        design_space_type, 
        output_flag
    )

    # 最適化用のフィルターのオブジェクトの作成と設定
    HyperTopOga.filter_setup!(
        opt_filter, 
        num_svalue_types, 
        length(elements), 
        physics, 
        radius, 
        design_space_type
    )

    # 微係数を求めたい値
    x = fill(0.5, length(elements) * num_svalue_types)

    # 評価関数の定義
    eval = HyperTopOga.StaticEndComplianceHyperElastic(
        x, 
        physics, 
        num_step
    )

    #----------------------------------------------------------
    # 
    #df_method=HyperTopOga.AutomaticDiff()
    #df_method=HyperTopOga.AutomaticAdjointDiff()
    df_method=HyperTopOga.AdjointDiff()
    f0, df1 = HyperTopOga.compute_f_and_grad(
        x, eval, physics, opt_settings, opt_filter, df_method
    )

end
#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
main()
