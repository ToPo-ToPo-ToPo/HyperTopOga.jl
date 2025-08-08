
#--------------------------------------------------------------------
# Mooney-Rivlin Hyperelastic Material Model (Compressible)
#--------------------------------------------------------------------
function compute_constitutive_law(
    material::MooneyRivlinPlaneStrain, F::Matrix{Float64}
)::Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Float64}}

    # Constants
    onethd = 1.0 / 3.0
    twothd = 2.0 / 3.0

    # Material properties
    a1 = material.a1
    b1 = material.b1
    a2 = material.a2
    b2 = material.b2
    K = material.K

    # Initialize tensors
    FF = update_F_for_plane_strain(F)

    #
    eye = Matrix{Float64}(I, 3, 3)
    eye4 = zeros(Float64, 3, 3, 3, 3)
    CinCin = zeros(Float64, 3, 3, 3, 3)
    CimCim = zeros(Float64 ,3, 3, 3, 3)
    Pis0 = zeros(Float64, 3, 3, 3, 3)
    Pis1 = zeros(Float64, 3, 3, 3, 3)
    Pis2 = zeros(Float64, 3, 3, 3, 3)

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
        CinCin[i,j,k,l] = Cinv[i,j] * Cinv[k,l]
        CimCim[i,j,k,l] = Cinv[i,k] * Cinv[j,l]
        eye4[i,j,k,l] = eye[i,k] * eye[j,l]
        Pis0[i,j,k,l] = eye4[i,j,k,l] - onethd*Cinv[i,j]*CC[k,l]
        Pis1[i,j,k,l] = eye[i,l] * eye[j,k] - onethd*Cinv[k,l]*CC[i,j]
        Pis2[i,j,k,l] = CimCim[i,j,k,l] - onethd*CinCin[i,j,k,l]
    end

    # Invariants
    zbar01 = tr(Cbar)
    zbar02 = 0.5 * (zbar01^2 - tr(Cbar * Cbar))
    zinv03 = dj^2

    # Energy derivatives
    dw01 = a1 + 2 * a2 * (zbar01 - 3)
    dw02 = b1 + 2 * b2 * (zbar02 - 3)
    dw0101 = 2 * a2
    dw0202 = 2 * b2
    gam01 = 2 * (dw01 + zbar01 * dw02)
    gam02 = -2 * dw02
    gav01 = K * (dj - 1)
    gav02 = K

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

    # Create S bar matrix
    S_bar_bar = zeros(Float64, 4, 4)
    for i in 1 : 2
        for j in 1 : 2
            S_bar_bar[i, j] = SS[i, j]
            S_bar_bar[i+2, j+2] = SS[i, j]
        end
    end

    # Voigt S
    stress = zeros(Float64, 3)
    stress[1] = SS[1,1] 
    stress[2] = SS[2,2]
    stress[3] = SS[1,2]

    # Tangent stiffness tensor
    c2bar = zeros(Float64 ,3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        c2bar[i,j,k,l] = 4 * zinv03^(-twothd) * (dw02 * (eye[i,j]*eye[k,l] - eye[i,k]*eye[j,l]) +
            (dw0101 + (zbar02^2)*dw0202) * eye[i,j]*eye[k,l] + dw0202 * (Cbar[i,j]*Cbar[k,l] - zbar01 * (Cbar[i,j]*eye[k,l] + eye[i,j]*Cbar[k,l]))
        )
    end

    pppp = zeros(Float64, 3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        pppp[i,j,k,l] = sum(Pis0[i,j,ia,ib] * c2bar[ia,ib,ic,id] * Pis1[ic,id,k,l] for ia in 1:3, ib in 1:3, ic in 1:3, id in 1:3)
    end
    qqqq = ((2/3) * dot(Cbar, Sbar)) .* Pis2

    c2pk = zeros(Float64, 3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        avol = dj * ((gav01 + dj * gav02) * CinCin[i,j,k,l] - 2 * gav01 * CimCim[i,j,k,l])
        aiso = pppp[i,j,k,l] + qqqq[i,j,k,l] - twothd * (Cinv[i,j] * Siso[k,l] + Siso[i,j] * Cinv[k,l])
        c2pk[i,j,k,l] = avol + aiso
    end

    # Voigt map
    C_matrix = tensor_4th_to_matrix_2th(c2pk)

    for i in 1:6, j in i+1:6
        C_matrix[j,i] = C_matrix[i,j]
    end

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

    return stress, S_bar_bar, C_matrix_33
end
