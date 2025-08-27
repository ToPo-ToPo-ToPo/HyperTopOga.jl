
using Revise
using HyperTopOga
#-------------------------------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------------------------------
function main()
    #--------------------------------------------------------------------------------------------
    # 解析モデルの設定
    length_x = 200.0
    length_y = 100.0
    division_x = 100
    division_y = 50
    thickness = 1.0

    # 材料モデルの定義
    interpolate = HyperTopOga.DMOSIMP([0.0, 1.0, 2.0], 1.0)
    young = HyperTopOga.DMOSIMP([0.0, 5.0, 7.0], 3.0)
    poisson = HyperTopOga.Constant(0.3)
    density = HyperTopOga.Constant(1.0e-05)
    gamma_x = HyperTopOga.Constant(1.0)
    material_data = HyperTopOga.NeoHookeanPlaneStrainData(young, poisson, density, gamma_x)

    # 解析条件の設定
    num_step = 5
    linear_solver = HyperTopOga.DirectSolverUMFPACK()
    problem_type = HyperTopOga.TotalLagrange()
    method = HyperTopOga.NewtonRapson()
    output_flag = true
    output_file_name = "output/MMTO_hyper_elastic_comp"

    #--------------------------------------------------------------------------------------------
    # 最適化の設定
    num_svalue_types = 2
    xmin = 1.0e-03
    max_eval = 300
    optimizer = HyperTopOga.MyGCMMA(0.25, 1.08, 0.65, 0.25, 0.2)

    # For objective function
    df_method = HyperTopOga.AutomaticAdjointDiff()

    # For constraint
    V_limit = 0.2

    # For opt filter
    radius = 3.0
    opt_filter = HyperTopOga.HeavisideProjectionFilter()
    
    #--------------------------------------------------------------------------------------------
    # 解析モデルの作成
    nodes, connects = HyperTopOga.make_quad4_voxel_mesh(length_x, length_y, division_x, division_y)

    # 要素情報の作成
    elements = HyperTopOga.make_structural_ele_quad4(connects, material_data, problem_type, thickness)

    # 境界条件・ソース項の情報を設定
    bc = HyperTopOga.ConditionPack()

    # ディレクレ境界を設定
    range_min = [0.0, 0.0, 0.0]
    range_max = [0.0, 100.0, 0.0]
    dirichlet_flag = [true, true]
    u_x(x) = 0.0
    u_y(x) = 0.0
    dirichlet_values = [u_x, u_y]
    HyperTopOga.add_condition(
        bc, 
        HyperTopOga.make_dirichlet_voxel_model(
            nodes, 
            range_min, 
            range_max,
            num_step, 
            dirichlet_flag, 
            dirichlet_values
        )
    )

    # ディレクレ境界を設定
    range_min = [200.0, 0.0, 0.0]
    range_max = [200.0, 100.0, 0.0]
    dirichlet_flag = [true, false]
    HyperTopOga.add_condition(
        bc, 
        HyperTopOga.make_dirichlet_voxel_model(
            nodes, 
            range_min, 
            range_max,
            num_step, 
            dirichlet_flag, 
            dirichlet_values
        )
    )

    # ノイマン境界を設定
    range_min = [190.0, 100.0, 0.0]
    range_max = [200.0, 100.0, 0.0]
    # 荷重条件を指定
    traction_x(x) = 0.0
    traction_y(x) = -0.04e-00 * (x / Float64(num_step))
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

    # 設計変数の初期値
    x0 = fill(0.5, length(elements) * num_svalue_types)

    # 最適化アルゴリズムの定義と設計変数の数の設定
    HyperTopOga.set_optimizer!(optimizer, length(x0))

    # 設計変数に直接与えられる制約
    lower_bounds = fill(xmin, length(x0))
    upper_bounds = fill(1.0, length(x0))
    # 固体領域を設定
    range_min = [190.0, 95.0, 0.0]
    range_max = [200.0, 100.0, 0.0]
    for isvalue in 1 : num_svalue_types
        for (e, element) in enumerate(elements)
            # 設計変数位置を取得
            index = HyperTopOga.make_index(num_svalue_types, isvalue, e)
            # 要素中心を取得
            center_coord = HyperTopOga.compute_center_coordinate(physics.nodes, element.shape)
            # 範囲に含まれるかを確認
            if range_min[1] - 1.0e-05 < center_coord[1] < range_max[1] + 1.0e-05 && 
                range_min[2] - 1.0e-05 < center_coord[2] < range_max[2] + 1.0e-05 &&
                range_min[3] - 1.0e-05 < center_coord[3] < range_max[3] + 1.0e-05
                if isvalue == num_svalue_types
                    x0[index] = 1.0
                    lower_bounds[index] = 1.0 - 1.0e-04
                    upper_bounds[index] = 1.0
                else
                    x0[index] = xmin
                    lower_bounds[index] = xmin - 1.0e-04
                    upper_bounds[index] = xmin + 1.0e-04

                end
            end
        end
    end
    optimizer.model.lower_bounds = lower_bounds
    optimizer.model.upper_bounds = upper_bounds

    # 最適化アルゴリズムのオプションを設定
    HyperTopOga.set_option!(optimizer, max_eval)

    # 目的関数の定義と設定
    objective = HyperTopOga.StaticEndComplianceHyperElastic(x0, physics, num_step, method)
    HyperTopOga.set_objective_function!(optimizer.model, "min", objective, 1.0, df_method)

    # 体積制約関数を定義
    constraint = HyperTopOga.Volume(x0, 1)
    HyperTopOga.add_inequality_constraint!(optimizer.model, constraint, V_limit, HyperTopOga.AutomaticAdjointDiff())

    # 体積制約関数を定義
    constraint = HyperTopOga.Volume(x0, 2)
    HyperTopOga.add_inequality_constraint!(optimizer.model, constraint, V_limit, HyperTopOga.AutomaticAdjointDiff())

    #--------------------------------------------------------------------------------------------
    # Running Optimization
    history, minx = HyperTopOga.optimize!(optimizer, physics, opt_settings, opt_filter, x0)

    # 最終的な最適化結果の出力
    HyperTopOga.output_opt_result(minx, physics, opt_settings, opt_filter)

    # 評価関数の履歴を出力
    HyperTopOga.output(history, physics.output_file_name, length(optimizer.model.objective), length(optimizer.model.constraint))
end
#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
@time main()
