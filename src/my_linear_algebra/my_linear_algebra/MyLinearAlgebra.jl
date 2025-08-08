
#----------------------------------------------------------------------------
# 関数を定義して行列のdeterminantを計算する
# 行列式の計算 (2x2 または 3x3 の行列に対応)
#----------------------------------------------------------------------------
function my_det(matrix::Matrix{Float64})::Float64
    det = 0.0
    # Check 22
    if size(matrix) == (2, 2)
        det = matrix[1, 1] * matrix[2, 2] - matrix[1, 2] * matrix[2, 1]
    # Check 33
    elseif size(matrix) == (3, 3)
        det = matrix[1, 1] * (matrix[2, 2] * matrix[3, 3] - matrix[2, 3] * matrix[3, 2]) -
              matrix[1, 2] * (matrix[2, 1] * matrix[3, 3] - matrix[2, 3] * matrix[3, 1]) +
              matrix[1, 3] * (matrix[2, 1] * matrix[3, 2] - matrix[2, 2] * matrix[3, 1])
    else
        error("Unsupported matrix size for determinant calculation.")
    end
    #
    return det
end
#----------------------------------------------------------------------------
# ノルムの計算 (norm()が自動微分ライブラリEnzyme.jlに対応していないため)
#----------------------------------------------------------------------------
function my_norm(A::Vector{Float64})::Float64
    return sqrt(dot(A, A))
end
#----------------------------------------------------------------------------
# 逆行列の計算
#----------------------------------------------------------------------------
function my_inv(A0::Matrix{Float64})::Matrix{Float64}
    #
    A = copy(A0)
    #
    n = size(A)[1]
    invA = zero(A)
    for i in 1 : n
        invA[i, i] = 1.0
    end
    # ガウス・ジョルダン法を用いて逆行列を計算
    for i in 1 : n
        # 対角成分を1に正規化
        if A[i, i] == 0.0
            error("Error! A[i, i] = $(A[i, i]). Matrix is singular")
        end
        factor = A[i, i]
        A[i, :] = A[i, :] / factor
        invA[i, :] = invA[i, :] / factor

        # 他の行を消去
        for j in 1 : n
            if j == i
                continue
            end
            factor = A[j, i]
            A[j, :] = A[j, :] - factor * A[i, :]
            invA[j, :] = invA[j, :] - factor * invA[i, :]
        end
    end
    #
    return invA
    #return inv(A0)
end
#----------------------------------------------------------------------------
# 逆行列
#----------------------------------------------------------------------------
function my_inv22(A::Matrix{Float64})::Matrix{Float64}

    if size(A) != (2, 2)
        error("Error: Input matrix must be a 2x2 matrix.")
    end

    A_inv = similar(A)
    detA = A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
    
    if detA == 0.0
        error("Error: The determinant of the matrix is zero. It's not invertible.")
    end

    A_inv[1, 1] =   A[2, 2]
    A_inv[1, 2] = - A[1, 2]
    A_inv[2, 1] = - A[2, 1]
    A_inv[2, 2] =   A[1, 1]

    return A_inv / detA
end
#----------------------------------------------------------------------------
# 3x3行列の逆行列を計算する関数
#----------------------------------------------------------------------------
function my_inv33(matrix::Matrix{Float64})::Matrix{Float64}
    if size(matrix) != (3, 3)
        error("Error: Input matrix must be a 3x3 matrix.")
    end
    
    det = matrix[1,1]*(matrix[2,2]*matrix[3,3] - matrix[2,3]*matrix[3,2]) -
          matrix[1,2]*(matrix[2,1]*matrix[3,3] - matrix[2,3]*matrix[3,1]) +
          matrix[1,3]*(matrix[2,1]*matrix[3,2] - matrix[2,2]*matrix[3,1])

    if det == 0.0
        error("Error: The determinant of the matrix is zero. It's not invertible.")
    end
    
    inv_matrix = similar(matrix)
    inv_matrix[1,1] = (matrix[2,2]*matrix[3,3] - matrix[2,3]*matrix[3,2])
    inv_matrix[1,2] = (matrix[1,3]*matrix[3,2] - matrix[1,2]*matrix[3,3])
    inv_matrix[1,3] = (matrix[1,2]*matrix[2,3] - matrix[1,3]*matrix[2,2])
    inv_matrix[2,1] = (matrix[2,3]*matrix[3,1] - matrix[2,1]*matrix[3,3])
    inv_matrix[2,2] = (matrix[1,1]*matrix[3,3] - matrix[1,3]*matrix[3,1])
    inv_matrix[2,3] = (matrix[1,3]*matrix[2,1] - matrix[1,1]*matrix[2,3])
    inv_matrix[3,1] = (matrix[2,1]*matrix[3,2] - matrix[2,2]*matrix[3,1])
    inv_matrix[3,2] = (matrix[1,2]*matrix[3,1] - matrix[1,1]*matrix[3,2])
    inv_matrix[3,3] = (matrix[1,1]*matrix[2,2] - matrix[1,2]*matrix[2,1])
    
    return inv_matrix / det
end