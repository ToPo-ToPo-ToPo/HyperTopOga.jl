
#----------------------------------------------------------------
# 構造用Quad4要素の定義
#----------------------------------------------------------------
struct StructuralQuad4{ET<:AbstractElementType, MD<:AbstractMaterialData} <: StructuralElement
    id::Int64
    num_element_dof::Int64
    dof_list::Vector{Int64}
    shape::Quad4
    type::ET
    material::MD        # 材料モデルのデータ
    #----------------------------------------------------------------
    #
    #----------------------------------------------------------------
    function StructuralQuad4(
        id::Int64, 
        nodes::Vector{Node}, 
        type::ET, 
        material::MD, 
        thickness::Float64
        ) where {ET<:AbstractElementType, MD<:AbstractMaterialData}

        # 形状データ作成
        shape = Quad4(nodes, thickness)
        
        # 自由度リスト作成
        dof_list = Vector{Int64}(undef, 8)
        for (i, id) in enumerate(shape.connects)
            dof_x = 2 * id - 1
            dof_y = 2 * id
            # 節点の自由度番号を作成(nodesは要素節点のみなのでここでは4つしかない)
            if length(nodes[i].dof_ids) < shape.dimension
                push!(nodes[i].dof_ids, dof_x)
                push!(nodes[i].dof_ids, dof_y)
            end
            # 要素内自由度番号の作成
            dof_list[2*i-1] = dof_x
            dof_list[2*i  ] = dof_y
        end
        
        # データを作成
        return new{ET, MD}(id, 8, dof_list, shape, type, material)
    end
end
#----------------------------------------------------------------
# 形状関数の定義
#----------------------------------------------------------------
function make_N(element::StructuralQuad4, coordinate::Vector{Float64})::Matrix{Float64}
    # 形状関数の要素を取得
    N_comp = make_N_components(element.shape, coordinate)
    # 形状関数の成形
    N = zeros(Float64, 2, 8)
    for i in 1 : 4
        N[1, 2*i-1] = N_comp[i]
        N[2, 2*i]   = N_comp[i]
    end
    return N
end
#----------------------------------------------------------------
# Bマトリクスの作成
#----------------------------------------------------------------
function make_B(::StructuralQuad4, dNdx::Matrix{Float64})::Matrix{Float64}
    B = zeros(Float64, 3, 8)
    for i in 1 : 4
        B[1, 2*i-1] = dNdx[1, i]
        B[2, 2*i  ] = dNdx[2, i]
        B[3, 2*i-1] = dNdx[2, i]
        B[3, 2*i  ] = dNdx[1, i]
    end
    return B
end