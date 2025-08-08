
#-------------------------------------------------------------------------------
# 固体力学のモデル：型パラメータ化
#-------------------------------------------------------------------------------
struct StaticSolidMechanics{EV<:Vector{StructuralElement}, L<:AbstractLinearSolver} <: SolidMechanics
    #
    nodes::Vector{Node}                                # 節点情報
    elements::EV                                       # 要素情報（型パラメータ E）
    condition::ConditionPack                           # 境界条件情報
    #
    num_total_eq::Int64
    sparse_data::SparseMatrixCSCData                   # 疎行列に関するデータをまとめた構造体
    multi_color::Vector{MultiColorizedElementDataset}  # 並列計算用のデータ
    linear_solver::L                                   # 線形ソルバー（型パラメータ L）
    #
    output_file_name::String
    #---------------------------------------------------------------
    # 内部コンストラクタ
    #---------------------------------------------------------------
    function StaticSolidMechanics(
        nodes::Vector{Node}, 
        elements::EV, 
        condition::ConditionPack,
        linear_solver::L,
        output_file_name::String
    ) where {EV<:Vector{StructuralElement}, L<:AbstractLinearSolver}

        # 自由度番号の設定と方程式の数を設定
        num_total_eq = 0
        for node in nodes
            num_total_eq += length(node.dof_ids)
        end

        # 並列計算用のカラーリングデータ
        multi_color = create_multi_colored_elem_ids(nodes, elements)

        # 疎行列向けのデータを作成
        sparse_data = create_sparse_matrix_csc_data(nodes, elements, num_total_eq)

        #
        return new{EV, L}(
            nodes, 
            elements, 
            condition, 
            num_total_eq, 
            sparse_data, 
            multi_color, 
            linear_solver, 
            output_file_name
        )
    end
end
#----------------------------------------------------------------
# 方程式を解く
#----------------------------------------------------------------
function solve_static_linear_problem(
    physics::StaticSolidMechanics, 
    element_states::Vector{ElementState}, 
    istep::Int64
    )::Vector{Float64}

    # 支配方程式を作成
    lhs, rhs = make_governing_equation(physics, element_states, istep)

    # 境界条件の処理
    lhs_bc, rhs_bc = compute_dirichlet_bc(
        physics.condition.dirichlet, 
        physics.nodes, 
        physics.num_total_eq, 
        lhs, 
        rhs, 
        istep
    )

    # 方程式の解を導出
    U = linear_solve(physics.linear_solver, lhs_bc, rhs_bc)
    
    #
    return U
end