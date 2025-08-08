

#---------------------------------------------------------------------------------------------
# モデルのデータ
#---------------------------------------------------------------------------------------------
struct MooneyRivlinPlaneStrainData <: AbstractMaterialData
    a1::AbstractInterpolation
    a2::AbstractInterpolation
    b1::AbstractInterpolation
    b2::AbstractInterpolation
    K::AbstractInterpolation
    density::AbstractInterpolation
    gamma_x::AbstractInterpolation
end
#---------------------------------------------------------------------------------------------
# データセットから実際の材料モデルを作成する
#---------------------------------------------------------------------------------------------
function setup(data::MooneyRivlinPlaneStrainData, xe::Matrix{Float64})::ElementState
    #
    num_p = size(xe, 2)
    materials = Vector{MooneyRivlinPlaneStrain}(undef, num_p)
    for ip in 1 : num_p
        #
        xip = xe[:, ip]
        #
        materials[ip] = MooneyRivlinPlaneStrain(
            compute(data.a1, xip), 
            compute(data.b1, xip),
            compute(data.a2, xip), 
            compute(data.b2, xip),
            compute(data.K, xip), 
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
struct MooneyRivlinPlaneStrain <: AbstractMaterial
    a1::Float64
    b1::Float64
    a2::Float64
    b2::Float64
    K::Float64
    density::Float64
    gamma_x::Float64
end
#---------------------------------------------------------------------------------------------
# ミーゼス応力を計算
#---------------------------------------------------------------------------------------------
function compute_mises_stress(::MooneyRivlinPlaneStrain, stress::Matrix{Float64})::Float64
    return sqrt((stress[1, 1]-stress[2, 2])^2.0 + (stress[2, 2]-stress[3, 3])^2.0 + (stress[3, 3]-stress[1, 1])^2.0 + 6.0*stress[1, 2]^2.0)
end
#---------------------------------------------------------------------------------------------
# 次のステップのために値を更新しておく
#---------------------------------------------------------------------------------------------
function update!(::MooneyRivlinPlaneStrain)::Nothing
    return nothing
end