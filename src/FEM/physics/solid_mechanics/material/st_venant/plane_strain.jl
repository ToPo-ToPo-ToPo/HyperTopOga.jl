

#---------------------------------------------------------------------------------------------
# モデルのデータ
#---------------------------------------------------------------------------------------------
struct StVenantPlaneStrainData <: AbstractMaterialData
    young::AbstractInterpolation
    poisson::AbstractInterpolation
    density::AbstractInterpolation
    gamma_x::AbstractInterpolation
end
#---------------------------------------------------------------------------------------------
# データセットから実際の材料モデルを作成する
#---------------------------------------------------------------------------------------------
function setup(dataset::StVenantPlaneStrainData, xe::Vector{Float64})::ElementState
    #
    num_p = size(xe, 2)
    containers = Vector{StVenantPlaneStrain}(undef, num_p)
    # Make containers
    for ip in 1 : num_p
        #
        xip = xe[:, ip]
        #
        containers[ip] = StVenantPlaneStrain(
            compute(dataset.young, xip), 
            compute(dataset.poisson, xip), 
            compute(dataset.density, xip), 
            compute(dataset.gamma_x, xip)
        )
    end
    #
    return ElementState(materials)
end
#---------------------------------------------------------------------------------------------
# 2D平面応力問題の線形弾性体モデル
#---------------------------------------------------------------------------------------------
struct StVenantPlaneStrain <: AbstractMaterial
    young::Float64
    poisson::Float64
    density::Float64
    gamma_x::Float64
end
#---------------------------------------------------------------------------------------------
# 構成則の計算 : 応力と接線剛性を計算する
#---------------------------------------------------------------------------------------------
function compute_constitutive_law(material::StVenantPlaneStrain, F::Matrix{Float64}
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
    # calc C
    C = F' * F
    # calc E
    E = 0.5 * (C - Imat)
    trE = E[1, 1] + E[2, 2] + E[3, 3]
    # calc 2-PK stress
    S = lambda * trE * Imat + 2.0 * mu * E
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
    C_mat = make_C(material)

    #
    return S_vec, S_bar, C_mat
end
#---------------------------------------------------------------------------------------------
# 構成則の計算 : 応力と接線剛性を計算する
#---------------------------------------------------------------------------------------------
function compute_stress(material::StVenantPlaneStrain, F::Matrix{Float64}
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
    # calc C
    C = F' * F
    # calc E
    E = 0.5 * (C - Imat)
    trE = E[1, 1] + E[2, 2] + E[3, 3]
    # calc 2-PK stress
    S = lambda * trE * Imat + 2.0 * mu * E
    #
    return S, F
end
#---------------------------------------------------------------------------------------------
# ミーゼス応力を計算
#---------------------------------------------------------------------------------------------
function compute_mises_stress(::StVenantPlaneStrain, stress::Matrix{Float64})::Float64
    return sqrt((stress[1, 1]-stress[2, 2])^2.0 + (stress[2, 2]-stress[3, 3])^2.0 + (stress[3, 3]-stress[1, 1])^2.0 + 6.0*stress[1, 2]^2.0)
end
#---------------------------------------------------------------------------------------------
# 弾性剛性マトリクスを作成
#---------------------------------------------------------------------------------------------
function make_C(material::StVenantPlaneStrain)::Matrix{Float64}
    return make_C_plane_strain(material.containers[ip].young, material.containers[ip].poisson)
end
#---------------------------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------------------------
function compute_linear_stress(material::StVenantPlaneStrain, strain::Vector{Float64})::Matrix{Float64}
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
function update!(::StVenantPlaneStrain, ::Int64)::Nothing
    return nothing
end
#---------------------------------------------------------------------------------------------
# 構成則の計算 : 応力と接線剛性を計算する
#---------------------------------------------------------------------------------------------
function compute_stress_sensitivity(
    mat_data::StVenantPlaneStrainData,
    material::StVenantPlaneStrain, 
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
    # calc C
    C = F' * F
    # calc E
    E = 0.5 * (C - Imat)
    trE = E[1, 1] + E[2, 2] + E[3, 3]
    # calc 2-PK stress
    dS = dlambda * trE * Imat + 2.0 * dmu * E
    #
    dS_voigt = Vector{Float64}(undef, 3)
    dS_voigt[1] = dS[1, 1]
    dS_voigt[2] = dS[2, 2]
    dS_voigt[3] = dS[1, 2]
    #
    return dS_voigt
end