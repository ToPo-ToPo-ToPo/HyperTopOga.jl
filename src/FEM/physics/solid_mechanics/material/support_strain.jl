
#---------------------------------------------------------------------------------------------
#
#---------------------------------------------------------------------------------------------
function voigt_to_tensor_2th_strain(vector::Vector{Float64})::Matrix{Float64}

    tensor = Matrix{Float64}(undef, 3, 3)
    tensor[1, 1] = vector[1]
    tensor[2, 2] = vector[2]
    tensor[3, 3] = vector[3]
    
    tensor[1, 2] = 0.5 * vector[4]
    tensor[1, 3] = 0.5 * vector[6]
    tensor[2, 1] = 0.5 * vector[4]
    tensor[2, 3] = 0.5 * vector[5]
    tensor[3, 1] = 0.5 * vector[6]
    tensor[3, 2] = 0.5 * vector[5]
    
    return tensor
end
#---------------------------------------------------------------------------------------------
# 平面ひずみ用にvoigt表記からテンソルに変換する
# 工学ひずみからの変換を行う
#---------------------------------------------------------------------------------------------
function voigt_to_tensor_2th_plane_strain(vector::Vector{Float64})::Matrix{Float64}
    
    tensor = Matrix{Float64}(undef, 3, 3)
    tensor[1, 1] = vector[1]
    tensor[2, 2] = vector[2]
    tensor[3, 3] = 0.0

    # 以下はvectorが工学ひずみで保存されていることが前提
    tensor[1, 2] = 0.5 * vector[3]
    tensor[1, 3] = 0.0
    tensor[2, 1] = 0.5 * vector[3]
    tensor[2, 3] = 0.0
    tensor[3, 1] = 0.0
    tensor[3, 2] = 0.0
    
    return tensor
end
#---------------------------------------------------------------------------------------------
# テンソルからvoigt表記に変更する
#---------------------------------------------------------------------------------------------
function tensor_2th_to_voigt_strain(tensor::Matrix{Float64})::Vector{Float64}

    vector = Vector{Float64}(undef, 6)

    vector[1] = tensor[1, 1]
    vector[2] = tensor[2, 2]
    vector[3] = tensor[3, 3]
    vector[4] = 2.0 * tensor[1, 2]
    vector[5] = 2.0 * tensor[2, 3]
    vector[6] = 2.0 * tensor[3, 1]

    return vector
end
#---------------------------------------------------------------------------------------------
# テンソルからvoigt表記に変更する
# ひずみ用
# 2D問題を対象
#---------------------------------------------------------------------------------------------
function tensor_2th_to_voigt_2d_strain(tensor::Matrix{Float64})::Vector{Float64}

    vector = Vector{Float64}(undef, 3)
    vector[1] = tensor[1, 1]
    vector[2] = tensor[2, 2]
    vector[3] = 2.0 * tensor[1, 2]

    return vector
end