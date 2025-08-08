
#--------------------------------------------------------------------
# Mooney-Rivlin Hyperelastic Material Model (Compressible)
#--------------------------------------------------------------------
function compute_stress(
    material::MooneyRivlinPlaneStrain, F::Matrix{Float64}
)::Tuple{Matrix{Float64}, Matrix{Float64}}

    # Constants
    onethd = 1.0 / 3.0

    # Material properties
    a1 = material.a1
    b1 = material.b1
    a2 = material.a2
    b2 = material.b2
    K = material.K

    # Initialize tensors
    FF = update_F_for_plane_strain(F)

    eye = Matrix{Float64}(I, 3, 3)
    Pis0 = zeros(Float64, 3, 3, 3, 3)

    # Determinant of F
    dj = my_det(FF)
    if dj <= 0.0
        error("Error: det(F)=0 at step $istp, element $inel")
    end

    # Compute Fbar
    Fbar = (dj)^(-onethd) .* FF

    # C and Cbar
    CC = FF' * FF
    Cbar = Fbar' * Fbar

    # Inverses
    Cinv = my_inv(CC)

    # Tensor constructions
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        Pis0[i,j,k,l] = eye[i,k]*eye[j,l] - onethd*Cinv[i,j]*CC[k,l]
    end

    # Invariants
    zbar01 = tr(Cbar)
    zbar02 = 0.5 * (zbar01^2 - tr(Cbar * Cbar))
    zinv03 = dj^2

    # Energy derivatives
    dw01 = a1 + 2 * a2 * (zbar01 - 3)
    dw02 = b1 + 2 * b2 * (zbar02 - 3)
    gam01 = 2 * (dw01 + zbar01 * dw02)
    gam02 = -2 * dw02
    gav01 = K * (dj - 1)

    # Sbar
    Sbar = gam01 .* eye .+ gam02 .* Cbar

    # Svol, Siso, SS
    Svol = (dj * gav01) .* Cinv
    Siso = zeros(Float64, 3, 3)
    fac = zinv03^(-onethd)
    for i in 1:3, j in 1:3
        Siso[i, j] = fac * sum(Pis0[i,j,ia,ib] * Sbar[ia,ib] for ia in 1:3, ib in 1:3)
    end
    SS = Siso .+ Svol

    return SS, FF
end

#---------------------------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------------------------
function compute_stress_sensitivity(
    mat_data::MooneyRivlinPlaneStrainData,
    ::MooneyRivlinPlaneStrain, 
    F::Matrix{Float64},
    xe::Vector{Float64},
    type::Int64
    )::Vector{Float64}

    # Constants
    onethd = 1.0 / 3.0

    # Material properties
    da1 = compute_derivative(mat_data.a1, xe, type)
    db1 = compute_derivative(mat_data.b1, xe, type)
    da2 = compute_derivative(mat_data.a2, xe, type)
    db2 = compute_derivative(mat_data.b2, xe, type)
    dK = compute_derivative(mat_data.K, xe, type)

    # Initialize tensors
    FF = update_F_for_plane_strain(F)
    #
    eye = Matrix{Float64}(I, 3, 3)
    Pis0 = zeros(Float64, 3, 3, 3, 3)

    # Determinant of F
    dj = my_det(FF)
    if dj <= 0.0
        error("Error: det(F)=0 at step $istp, element $inel")
    end

    # Compute Fbar
    Fbar = (dj)^(-onethd) .* FF

    # C and Cbar
    CC = FF' * FF
    Cbar = Fbar' * Fbar

    # Inverses
    Cinv = my_inv(CC)

    # Tensor constructions
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        Pis0[i,j,k,l] = eye[i,k]*eye[j,l] - onethd*Cinv[i,j]*CC[k,l]
    end

    # Invariants
    zbar01 = tr(Cbar)
    zbar02 = 0.5 * (zbar01^2 - tr(Cbar * Cbar))
    zinv03 = dj^2

    # Energy derivatives
    dw01 = da1 + 2 * da2 * (zbar01 - 3)
    dw02 = db1 + 2 * db2 * (zbar02 - 3)
    gam01 = 2 * (dw01 + zbar01 * dw02)
    gam02 = -2 * dw02
    gav01 = dK * (dj - 1)

    # Sbar
    Sbar = gam01 .* eye .+ gam02 .* Cbar

    # Svol, Siso, SS
    Svol = (dj * gav01) .* Cinv
    Siso = zeros(Float64, 3, 3)
    fac = zinv03^(-onethd)
    for i in 1:3, j in 1:3
        Siso[i, j] = fac * sum(Pis0[i,j,ia,ib] * Sbar[ia,ib] for ia in 1:3, ib in 1:3)
    end
    SS = Siso .+ Svol
    #
    dS_voigt = Vector{Float64}(undef, 3)
    dS_voigt[1] = SS[1, 1]
    dS_voigt[2] = SS[2, 2]
    dS_voigt[3] = SS[1, 2]

    return dS_voigt

end
