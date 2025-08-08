
#---------------------------------------------------------------------------------------------
# 4階のテンソルから2階のテンソルへ変換する
#---------------------------------------------------------------------------------------------
function tensor_4th_to_matrix_2th(tensor_4th::Array{Float64})::Matrix{Float64}
    # make IJ data
    IJdata = Matrix{Int64}(undef, 3, 3)
    IJdata[1, 1] = 1
    IJdata[1, 2] = 4
    IJdata[1, 3] = 6
 
    IJdata[2, 1] = 4
    IJdata[2, 2] = 2
    IJdata[2, 3] = 5
 
    IJdata[3, 1] = 6
    IJdata[3, 2] = 5
    IJdata[3, 3] = 3
 
    # make 2th tensor
    matrix_2th = Matrix{Float64}(undef, 6, 6)
    for i in 1 : 3
        for j in 1 : 3
            I = IJdata[i, j]
            for k in 1 : 3
                for l in 1 : 3
                    J = IJdata[k, l]
                    matrix_2th[I, J] = 0.5 * (tensor_4th[i, j, k, l] + tensor_4th[i, j, l, k])
                end
            end
        end
    end
    return matrix_2th
end
#---------------------------------------------------------------------------------------------
# 2階のテンソルから4階のテンソルへ変換する
#---------------------------------------------------------------------------------------------
function matrix_2th_to_tensor_4th(matrix_2th::Matrix{Float64})::Matrix{Float64}

    # make IJ data
    IJdata = Matrix{Int64}(undef, 3, 3)
    IJdata[1, 1] = 1
    IJdata[1, 2] = 4
    IJdata[1, 3] = 6
 
    IJdata[2, 1] = 4
    IJdata[2, 2] = 2
    IJdata[2, 3] = 5
 
    IJdata[3, 1] = 6
    IJdata[3, 2] = 5
    IJdata[3, 3] = 3

    # Tensor to matrix
    matrix_2th_tmp = copy(matrix_2th)
    for i in 4 : 6
        for j in 4 : 6
            matrix_2th_tmp[i, j] = 0.5 * matrix_2th_tmp[i, j]
        end
    end
    #
    tensor_4th = Array{Float64}(undef, 3, 3, 3, 3)
    for i in 1 : 3
        for j in 1 : 3
            for k in 1 : 3
                for l in 1 : 3
                    tensor_4th[i, j, k, l] = matrix_2th_tmp[IJdata[i, j], IJdata[k, l]]
                end
            end
        end
    end
    #
    return tensor_4th
end