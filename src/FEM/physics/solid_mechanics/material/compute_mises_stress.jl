
#--------------------------------------------------------------------
# ミーゼス応力を計算
#--------------------------------------------------------------------
function compute_mises_stress(stress::Matrix{Float64})::Float64
    stress_voigt = tensor_2th_to_voigt(stress)
    sig_xx, sig_yy, sig_zz, tau_xy, tau_yz, tau_zx = stress_voigt
    return sqrt(0.5 * ((sig_xx-sig_yy)^2.0 + (sig_yy-sig_zz)^2.0 + (sig_zz-sig_xx)^2.0 + 6.0*(tau_xy^2.0+tau_yz^2.0+tau_zx^2.0)))
end
#--------------------------------------------------------------------
# ミーゼス応力を計算 2D 平面ひずみ用
#--------------------------------------------------------------------
function compute_mises_stress_2d_plane_strain(stress::Vector{Float64})::Float64
    sig_xx, sig_yy, tau_xy, sig_zz = stress
    return sqrt(0.5 * ((sig_xx-sig_yy)^2.0 + (sig_yy-sig_zz)^2.0 + (sig_zz-sig_xx)^2.0 + 6.0*tau_xy^2.0))
end
#---------------------------------------------------------------------------------------------
# ミーゼス応力を計算 3D
#---------------------------------------------------------------------------------------------
function compute_mises_stress(stress::Vector{Float64})::Float64
    sig_xx, sig_yy, sig_zz, tau_xy, tau_yz, tau_zx = stress
    term1 = (sig_xx-sig_yy)^2.0 + (sig_yy-sig_zz)^2.0 + (sig_zz-sig_xx)^2.0
    term2 = 6.0*(tau_xy^2.0 + tau_yz^2.0 + tau_zx^2.0)
    return sqrt(term1 + term2)
end