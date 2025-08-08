

#---------------------------------------------------------------------------------------------
# モデルのデータ
#---------------------------------------------------------------------------------------------
struct NeoHookeanPlaneStrainData <: AbstractMaterialData
    young::AbstractInterpolation
    poisson::AbstractInterpolation
    density::AbstractInterpolation
    gamma_x::AbstractInterpolation
end
#---------------------------------------------------------------------------------------------
# データセットから実際の材料モデルを作成する
#---------------------------------------------------------------------------------------------
function setup(data::NeoHookeanPlaneStrainData, xe::Matrix{Float64})::ElementState
    #
    num_p = size(xe, 2)
    materials = Vector{NeoHookeanPlaneStrain}(undef, num_p)
    # Make materials
    for ip in 1 : num_p
        #
        xip = xe[:, ip]
        #
        materials[ip] = NeoHookeanPlaneStrain(
            compute(data.young, xip), 
            compute(data.poisson, xip), 
            compute(data.density, xip), 
            compute(data.gamma_x, xip)
        )
    end
    #
    return ElementState(materials)
end
#---------------------------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------------------------
struct NeoHookeanPlaneStrain <: AbstractMaterial
    young::Float64
    poisson::Float64
    density::Float64
    gamma_x::Float64
end
#---------------------------------------------------------------------------------------------
# 構成則の計算 : 応力と接線剛性を計算する
#---------------------------------------------------------------------------------------------
function compute_constitutive_law(
    material::NeoHookeanPlaneStrain, 
    F::Matrix{Float64}
    )::Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Float64}}

    # 材料物性の取得
    young = material.young
    poisson = material.poisson
    # lame定数
    lambda  = (poisson * young) / ((1.0+poisson) * (1.0-2.0*poisson))
    mu = 0.5 * young / (1.0 + poisson)
    # 単位行列
    Imat = Matrix{Float64}(I, 3, 3)
    # cal deformation gradient F
    F = update_F_for_plane_strain(F)
    # cal determinant of F
    detF = my_det(F)
    # calc C
    C = F' * F
    # calc C inverse
    C_inv = my_inv33(C)
    # calc 2-PK stress
    S = mu * (Imat - C_inv) + lambda * detF * (detF - 1.0) * C_inv
    # create S bar matrix
    S_bar = zeros(Float64, 4, 4)
    for i in 1 : 2
        for j in 1 : 2
            S_bar[i, j] = S[i, j]
            S_bar[i+2, j+2] = S[i, j]
        end
    end
    # make stress vector
    S_vec = Vector{Float64}(undef, 3)
    S_vec[1] = S[1, 1]
    S_vec[2] = S[2, 2]
    S_vec[3] = S[1, 2]
    # make tangent modulus
    A1 = lambda * detF * (2.0 * detF - 1.0)
    A2 = 2.0 * (mu - lambda * detF * (detF - 1.0))
    C_tensor = Array{Float64}(undef, 3, 3, 3, 3)
    for i in 1 : 3
        for j in 1 : 3
            for k in 1 : 3
                for l in 1 : 3
                    C_tensor[i, j, k, l] = A1 * C_inv[i, j] * C_inv[k, l] + A2 * 0.5 * (C_inv[i, k] * C_inv[j, l] + C_inv[i, l] * C_inv[j, k])
                end
            end
        end
    end
    # 4th tensor -> matrix
    C_matrix = tensor_4th_to_matrix_2th(C_tensor)
    # data set
    mat_tmp = Vector{Float64}(undef, 3)
    mat_tmp[1] = C_matrix[1, 4]
    mat_tmp[2] = C_matrix[2, 4]
    mat_tmp[3] = C_matrix[4, 4]

    # update
    C_matrix_33 = Matrix{Float64}(undef, 3, 3)
    for i in 1 : 3
        for j in 1 : 3
            C_matrix_33[i, j] = C_matrix[i, j]
        end
    end
    for i in 1 : 3
        C_matrix_33[i, 3] = mat_tmp[i]
        C_matrix_33[3, i] = mat_tmp[i]
    end

    #
    return S_vec, S_bar, C_matrix_33
end
#---------------------------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------------------------
function compute_stress(
    material::NeoHookeanPlaneStrain, 
    F::Matrix{Float64}
    )::Tuple{Matrix{Float64}, Matrix{Float64}}

    # 材料物性の取得
    young = material.young
    poisson = material.poisson
    # lame定数
    lambda  = (poisson * young) / ((1.0+poisson) * (1.0-2.0*poisson))
    mu = 0.5 * young / (1.0 + poisson)
    # 単位行列
    Imat = Matrix{Float64}(I, 3, 3)
    # cal deformation gradient F
    F = update_F_for_plane_strain(F)
    # cal determinant of F
    detF = my_det(F)
    # calc C
    C = F' * F
    # calc C inverse
    C_inv = my_inv33(C)
    # calc 2-PK stress
    S = mu * (Imat - C_inv) + lambda * detF * (detF - 1.0) * C_inv
    #
    return S, F
end
#---------------------------------------------------------------------------------------------
# ミーゼス応力を計算
#---------------------------------------------------------------------------------------------
function compute_mises_stress(::NeoHookeanPlaneStrain, stress::Matrix{Float64})::Float64
    return sqrt((stress[1, 1]-stress[2, 2])^2.0 + (stress[2, 2]-stress[3, 3])^2.0 + (stress[3, 3]-stress[1, 1])^2.0 + 6.0*stress[1, 2]^2.0)
end
#---------------------------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------------------------
function make_C(material::NeoHookeanPlaneStrain)::Matrix{Float64}
    return make_C_plane_strain(material.young, material.poisson)
end
#---------------------------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------------------------
function make_C0(material::NeoHookeanPlaneStrain)::Matrix{Float64}
    return make_C0_plane_strain(material.poisson)
end
#---------------------------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------------------------
function compute_linear_stress(material::NeoHookeanPlaneStrain, strain::Vector{Float64})::Matrix{Float64}
    # 材料物性の取得
    young = material.young
    poisson = material.poisson
    # lame定数
    lambda  = (poisson * young) / ((1.0+poisson) * (1.0-2.0*poisson))
    mu = 0.5 * young / (1.0 + poisson)
    # 単位行列
    Imat = Matrix{Float64}(I, 3, 3)
    # strain
    E = zeros(Float64, 3, 3)
    E[1, 1] = strain[1]
    E[1, 2] = strain[3]
    E[2, 1] = strain[3]
    E[2, 2] = strain[2]
    
    #
    S = lambda * tr(E) * Imat + 2.0 * mu * E
    #
    return S
end
#---------------------------------------------------------------------------------------------
# 次のステップのために値を更新しておく
#---------------------------------------------------------------------------------------------
function update!(::NeoHookeanPlaneStrain)::Nothing
    return nothing
end
#---------------------------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------------------------
function compute_stress_sensitivity(
    mat_data::NeoHookeanPlaneStrainData,
    material::NeoHookeanPlaneStrain, 
    F::Matrix{Float64},
    xe::Vector{Float64},
    type::Int64
    )::Vector{Float64}

    # 材料物性の取得
    dyoung = compute_derivative(mat_data.young, xe, type)
    poisson = material.poisson
    # lame定数
    dlambda  = (poisson * dyoung) / ((1.0+poisson) * (1.0-2.0*poisson))
    dmu = 0.5 * dyoung / (1.0 + poisson)
    # 単位行列
    Imat = Matrix{Float64}(I, 3, 3)
    # cal deformation gradient F
    F = update_F_for_plane_strain(F)
    # cal determinant of F
    detF = my_det(F)
    # calc C
    C = F' * F
    # calc C inverse
    C_inv = my_inv33(C)
    # calc 2-PK stress
    dS = dmu * (Imat - C_inv) + dlambda * detF * (detF - 1.0) * C_inv
    #
    dS_voigt = Vector{Float64}(undef, 3)
    dS_voigt[1] = dS[1, 1]
    dS_voigt[2] = dS[2, 2]
    dS_voigt[3] = dS[1, 2]
    #
    return dS_voigt
end