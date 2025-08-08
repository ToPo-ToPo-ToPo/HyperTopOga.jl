

#----------------------------------------------------------------------------
#　境界条件をまとめた構造体
#----------------------------------------------------------------------------
struct ConditionPack
    dirichlet::Vector{DirichletBc}
    neumann::Vector{AbstractNeumannBc}
    #------------------------------------------------------------
    # 内部コンストラクタ
    #------------------------------------------------------------
    function ConditionPack()
        return new(
            Vector{DirichletBc}(), 
            Vector{AbstractNeumannBc}()
        )
    end
end
#----------------------------------------------------------------
# 条件を追加する関数
#----------------------------------------------------------------
function add_condition(pack::ConditionPack, bc::DirichletBc)
    push!(pack.dirichlet, bc)
end
#----------------------------------------------------------------
# 条件を追加する関数
#----------------------------------------------------------------
function add_condition(pack::ConditionPack, bc::NeumannBc)
    push!(pack.neumann, bc)
end