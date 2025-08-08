
#----------------------------------------------------------------
# BLマトリクスの作成
#----------------------------------------------------------------
function make_BL_and_BNL(
    ::StructuralQuad4, 
    dNdx::Matrix{Float64}, 
    F::Matrix{Float64}
    )::Tuple{Matrix{Float64}, Matrix{Float64}}
    
    # B_lマトリクスの作成
    B_l = zeros(Float64, 3, 8)
    for i in 1 : 4
        B_l[1, 2*i-1] = F[1, 1] * dNdx[1, i]
        B_l[1, 2*i  ] = F[2, 1] * dNdx[1, i]
       
        B_l[2, 2*i-1] = F[1, 2] * dNdx[2, i]
        B_l[2, 2*i  ] = F[2, 2] * dNdx[2, i]
            
        B_l[3, 2*i-1] = F[1, 2] * dNdx[1, i] + F[1, 1] * dNdx[2, i]
        B_l[3, 2*i  ] = F[2, 2] * dNdx[1, i] + F[2, 1] * dNdx[2, i]
    end

    # B_nlマトリクスの作成
    B_nl = zeros(Float64, 4, 8)
    for i in 1 : 4
        B_nl[1, 2*i-1] = dNdx[1, i]
        B_nl[2, 2*i-1] = dNdx[2, i]
        
        B_nl[3, 2*i  ] = dNdx[1, i]
        B_nl[4, 2*i  ] = dNdx[2, i]
    end

    return B_l, B_nl
end
#-----------------------------------------------------------------------------------------------------
# 2Dの変形勾配テンソルFを作成
#-----------------------------------------------------------------------------------------------------
function make_F(element::StructuralQuad4, dNdx::Matrix{Float64}, Ue::Vector{Float64}, ::T
    )::Matrix{Float64} where {T<:AbstractTotalLagrange}
    # 初期化
    dudX = zeros(Float64, 2, 2)
    Imat = Matrix{Float64}(I, 2, 2)
    # make dudX
    for i in 1 : length(element.shape.connects)
        # for x
        dudX[1, 1] += dNdx[1, i] * Ue[2*i-1]
        dudX[1, 2] += dNdx[2, i] * Ue[2*i-1]
        # for y
        dudX[2, 1] += dNdx[1, i] * Ue[2*i]
        dudX[2, 2] += dNdx[2, i] * Ue[2*i]
    end
    #
    return Imat + dudX
end