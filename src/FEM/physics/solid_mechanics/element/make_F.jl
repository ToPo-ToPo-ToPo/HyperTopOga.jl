
#------------------------------------------------------------------------
# 平面ひずみモデルの変形勾配テンソルFを作成
#------------------------------------------------------------------------
function update_F_for_plane_strain(F::Matrix{Float64})::Matrix{Float64}
    # 初期化
    F_plane_strain = zeros(Float64, 3, 3)
    for i in 1 : 2
        for j in 1 : 2
            F_plane_strain[i, j] = F[i, j]
        end
    end
    # Z方向を考慮
    F_plane_strain[3, 3] = 1.0
    return F_plane_strain
end
