
#----------------------------------------------------------------------------
# Solver
#----------------------------------------------------------------------------
function my_linear_solve_UMFPACK(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64})::Vector{Float64}
    prob = LinearProblem(A, b)
    sol = solve(prob, UMFPACKFactorization())
    return sol.u
end
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
function EnzymeRules.augmented_primal(
    config::EnzymeRules.RevConfig, 
    func::Const{typeof(my_linear_solve_UMFPACK)}, 
    ::Type{RT}, 
    A::Annotation{SparseMatrixCSC{Float64, Int64}}, 
    b::Annotation{Vector{Float64}}
    ) where RT

    # 関数の戻り値
    primal = if EnzymeRules.needs_primal(config)
        func.val(A.val, b.val)
    else
        func.val(A.val, b.val)
    end

    # 引数A
    cache_A = if EnzymeRules.overwritten(config)[2]
        copy(A.val)
    else
        nothing
    end

    # 引数B
    cache_b = if EnzymeRules.overwritten(config)[3]
        copy(b.val)
    else
        nothing
    end

    # dy_bar
    dc = if EnzymeRules.width(config) == 1
        zero(primal)
    else
        ntuple(Val(EnzymeRules.width(config))) do i
            Base.@_inline_meta
            zero(primal)
        end
    end

    # 微分値を計算する際に必要なデータをキャッシュとして保存
    cache = NamedTuple{(
        Symbol("1"), Symbol("2"), Symbol("3"), Symbol("4")), 
        Tuple{typeof(primal), typeof(dc), typeof(cache_A), typeof(cache_b)}}(
        (primal, dc, cache_A, cache_b)
    )

    # 関数の計算結果と戻り値に対する微分値を入力するデータ構造、微分値計算に必要なデータを返す
    return EnzymeRules.AugmentedReturn{typeof(primal), typeof(dc), typeof(cache)}(primal, dc, cache)
end
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
function EnzymeRules.reverse(
    config::EnzymeRules.RevConfig, 
    func::Const{typeof(my_linear_solve_UMFPACK)}, 
    ::Type{RT}, 
    cache,
    A::Annotation{<:SparseMatrixCSC{Float64, Int64}}, 
    b::Annotation{<:Vector{Float64}}
    ) where RT

    # 関数の計算にて使用したキャッシュを取得
    y, dys, cache_A, cache_b = cache
    
    if !EnzymeRules.overwritten(config)[2]
        cache_A = A.val
    end

    if !EnzymeRules.overwritten(config)[3]
        cache_b = b.val
    end

    if EnzymeRules.width(config) == 1
        dys = (dys,)
    end

    dAs = if EnzymeRules.width(config) == 1
        if typeof(A) <: Const
            (nothing,)
        else
            (A.dval,)
        end
    else
        if typeof(A) <: Const
            ntuple(Val(EnzymeRules.width(config))) do i
                Base.@_inline_meta
                nothing
            end
        else
            A.dval
        end
    end

    dbs = if EnzymeRules.width(config) == 1
        if typeof(b) <: Const
            (nothing,)
        else
            (b.dval,)
        end
    else
        if typeof(b) <: Const
            ntuple(Val(EnzymeRules.width(config))) do i
                Base.@_inline_meta
                nothing
            end
        else
            b.dval
        end
    end

    # 関数に入力される変数ごとの微分値を計算
    # y = inv(A) * b
    # dA −= z * y^T  
    # db += z, where z = inv(A^T) * dy
    for (dA, db, dy) in zip(dAs, dbs, dys)
        #
        z = func.val(sparse(transpose(cache_A)), dy)
        #
        if !(typeof(A) <: Const)
            # 疎行列はスパース性が変更されないようにする
            n, pp, ii = cache_A.n, cache_A.colptr, cache_A.rowval
            for j = 1 : n
                for p = pp[j] : pp[j+1] - 1
                    # 元のdAの情報を直接書き換える
                    dA.nzval[p] -= z[ii[p]] * y[j]
                end
            end
        end
        if !(typeof(b) <: Const)
            # .=を使用して値を更新
            db .+= z
        end
        # .=を使用して値を更新
        dy .= eltype(dy)(0)
    end

    return (nothing, nothing)
end