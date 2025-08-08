

#--------------------------------------------------------------------
# 相当ひずみを計算
#--------------------------------------------------------------------
function compute_equivalent_strain(strain::Matrix{Float64})::Float64
    val1 = sqrt(2.0) / 3.0
    val2 = (strain[1, 1] - strain[2, 2])^2.0 + (strain[2, 2] - strain[3, 3])^2.0 + (strain[3, 3] - strain[1, 1])^2.0
    val3 = 6.0 * (strain[1, 2]^2.0 + strain[2, 3]^2.0 + strain[3, 1]^2.0)
    return val1 * sqrt(val2 + val3)  
end
#--------------------------------------------------------------------
# 相当ひずみを計算 2D
#--------------------------------------------------------------------
function compute_equivalent_strain_plane_strain(strain::Vector{Float64})::Float64
    val1 = sqrt(2.0) / 3.0
    val2 = (strain[1] - strain[2])^2.0 + strain[2]^2.0 + strain[1]^2.0
    val3 = 6.0 * ((0.5*strain[3])^2.0)
    return val1 * sqrt(val2 + val3)  
end
#--------------------------------------------------------------------
# 相当ひずみを計算 3D
#--------------------------------------------------------------------
function compute_equivalent_strain(strain::Vector{Float64})::Float64
    val1 = sqrt(2.0) / 3.0
    val2 = (strain[1] - strain[2])^2.0 + (strain[2] - strain[3])^2.0 + (strain[3] - strain[1])^2.0
    val3 = 6.0 * (strain[4]^2.0 + strain[5]^2.0 + strain[6]^2.0)
    return val1 * sqrt(val2 + val3)  
end