
import Base.+
import Base.-
#----------------------------------------------------------------------------
# すでに値が保存されている要素の値を更新
#----------------------------------------------------------------------------
function +(A::SparseMatrixCSC{Float64, Int64}, B::SparseMatrixCSC{Float64, Int64})::SparseMatrixCSC{Float64, Int64}
    if getproperty.(Ref(A), (:m,:n,:colptr,:rowval)) == getproperty.(Ref(B), (:m,:n,:colptr,:rowval))
        nzval = A.nzval + B.nzval
        return SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, nzval)
    else
        if A.nzval > B.nzval
            m, n, colptr, rowval = A.m, A.n, A.colptr, A.rowval
            sp = SparseMatrixCSC(m, n, colptr, rowval, zero(A.nzval))
            for j in 1 : n
                for i in 1 : m
                    sp[i, j] = A[i, j] + B[i, j]
                end
            end
            return sp
        else
            m, n, colptr, rowval = B.m, B.n, B.colptr, B.rowval
            sp = SparseMatrixCSC(m, n, colptr, rowval, zero(B.nzval))
            for j in 1 : n
                for i in 1 : m
                    sp[i, j] = A[i, j] + B[i, j]
                end
            end
            return sp
        end
    end
end
#----------------------------------------------------------------------------
# すでに値が保存されている要素の値を更新
#----------------------------------------------------------------------------
function -(A::SparseMatrixCSC{Float64, Int64}, B::SparseMatrixCSC{Float64, Int64})::SparseMatrixCSC{Float64, Int64}
    if getproperty.(Ref(A), (:m,:n,:colptr,:rowval)) == getproperty.(Ref(B), (:m,:n,:colptr,:rowval))
        nzval = A.nzval .- B.nzval
        return SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, nzval)
    else
        if A.nzval > B.nzval
            m, n, colptr, rowval = A.m, A.n, A.colptr, A.rowval
            sp = SparseMatrixCSC(m, n, colptr, rowval, zero(A.nzval))
            for j in 1 : n
                for i in 1 : m
                    sp[i, j] = A[i, j] - B[i, j]
                end
            end
            return sp
        else
            m, n, colptr, rowval = B.m, B.n, B.colptr, B.rowval
            sp = SparseMatrixCSC(m, n, colptr, rowval, zero(B.nzval))
            for j in 1 : n
                for i in 1 : m
                    sp[i, j] = A[i, j] - B[i, j]
                end
            end
            return sp
        end
    end
end
#------------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------------
function make_diag(A::SparseMatrixCSC{Float64, Int64})::Vector{Float64}
    diag = Vector{Int64}(undef, A.n)
    @inbounds for i in eachindex(diag)
        for k in A.colptr[i] : A.colptr[i+1] - 1
            if A.rowval[k] == i
                diag[i] = k
            end
        end
    end
    return diag
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(make_diag), A::SparseMatrixCSC{Float64, Int64})
    return nothing
end