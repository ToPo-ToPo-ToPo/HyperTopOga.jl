

#---------------------------------------------------------------------------------------------
#
#---------------------------------------------------------------------------------------------
function voigt_to_tensor_2th(vector::Vector{Float64})::Matrix{Float64}
    
    tensor = Matrix{Float64}(undef, 3, 3)
    
    tensor[1, 1] = vector[1]
    tensor[2, 2] = vector[2]
    tensor[3, 3] = vector[3]

    tensor[1, 2] = vector[4]
    tensor[1, 3] = vector[6]
    tensor[2, 1] = vector[4]
    tensor[2, 3] = vector[5]
    tensor[3, 1] = vector[6]
    tensor[3, 2] = vector[5]

    return tensor
end
#---------------------------------------------------------------------------------------------
#
#---------------------------------------------------------------------------------------------
function voigt_to_tensor_2th_2d(vector::Vector{Float64})::Matrix{Float64}
    
    tensor = Matrix{Float64}(undef, 3, 3)
    
    tensor[1, 1] = vector[1]
    tensor[2, 2] = vector[2]
    tensor[3, 3] = vector[4]

    tensor[1, 2] = vector[3]
    tensor[1, 3] = 0.0
    tensor[2, 1] = vector[3]
    tensor[2, 3] = 0.0
    tensor[3, 1] = 0.0
    tensor[3, 2] = 0.0

    return tensor
end
#---------------------------------------------------------------------------------------------
# テンソルからvoigt表記に変更する
# 応力の変換に使用
#---------------------------------------------------------------------------------------------
function tensor_2th_to_voigt(tensor::Matrix{Float64})::Vector{Float64}

    vector = Vector{Float64}(undef, 6)

    vector[1] = tensor[1, 1]
    vector[2] = tensor[2, 2]
    vector[3] = tensor[3, 3]
    vector[4] = tensor[1, 2]
    vector[5] = tensor[2, 3]
    vector[6] = tensor[3, 1]

    return vector
end
#---------------------------------------------------------------------------------------------
# テンソルからvoigt表記に変更する
# 応力の変換に使用
#---------------------------------------------------------------------------------------------
function tensor_2th_to_voigt_plane_strain(tensor::Matrix{Float64})::Vector{Float64}

    vector = Vector{Float64}(undef, 3)

    vector[1] = tensor[1, 1]
    vector[2] = tensor[2, 2]
    vector[3] = tensor[1, 2]

    return vector
end
#---------------------------------------------------------------------------------------------
# テンソルからvoigt表記に変更する
# zz成分を考慮
#---------------------------------------------------------------------------------------------
function tensor_2th_to_voigt_plane_strain_with_zz(tensor::Matrix{Float64})::Vector{Float64}

    vector = Vector{Float64}(undef, 4)

    vector[1] = tensor[1, 1]
    vector[2] = tensor[2, 2]
    vector[3] = tensor[1, 2]
    vector[4] = tensor[3, 3]

    return vector
end
#---------------------------------------------------------------------------------------------
# テンソルからvoigt表記に変更する
# 応力用（2D問題の内力計算に使用する）
#---------------------------------------------------------------------------------------------
function tensor_2th_to_voigt_2d_stress(tensor::Matrix{Float64})::Vector{Float64}

    vector = Vector{Float64}(undef, 3)
    
    vector[1] = tensor[1, 1]
    vector[2] = tensor[2, 2]
    vector[3] = tensor[1, 2]

    return vector
end