
#----------------------------------------------------------------
# 線積分におけるjacobianの計算
#----------------------------------------------------------------
function make_boundary_detJ(::Quad4, boundary_id::Int64, J::Matrix{Float64})::Float64
    # 面番号によって分岐
    detJ = 0.0
    if boundary_id == 1 || boundary_id == 3
        detJ = sqrt(J[1, 2]^2.0 + J[2, 2]^2.0)
    elseif boundary_id == 2 || boundary_id == 4
        detJ = sqrt(J[1, 1]^2.0 + J[2, 1]^2.0)
    end
    # Jacobianの正負チェック
    if detJ <= 0.0 
        error("Jacobian is negative. Please change the nodal setting.")
    end
    return detJ
end