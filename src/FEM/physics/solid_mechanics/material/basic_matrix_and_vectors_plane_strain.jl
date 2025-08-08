
#---------------------------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------------------------
function make_I(n::Int64)::Matrix{Float64}
    Imat = Matrix{Float64}(I, n, n)
    return Imat
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(make_I), n::Int64)
    return nothing
end
#---------------------------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------------------------
function make_Is_tensor(n::Int64)::Array{Float64, 4}
    #
    I = make_I(n)
    #
    Is = Array{Float64}(undef, n, n, n, n)
    for i in 1 : n
        for j in 1 : n
            for k in 1 : n
                for l in 1 : n
                    Is[i, j, k, l] = 0.5 * (I[i, k] * I[j, l] + I[i, l] * I[j, k])
                end
            end
        end
    end
    #
    return Is
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(make_Is_tensor), n::Int64)
    return nothing
end
#---------------------------------------------------------------------------------------------
# 2階の恒等テンソルのベクトル形式
#---------------------------------------------------------------------------------------------
function make_I_vector_plane_strain()::Vector{Float64}
    SOID = Vector{Float64}(undef, 4)
    SOID[1] = 1.0
    SOID[2] = 1.0
    SOID[3] = 0.0
    SOID[4] = 1.0
    return SOID
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(make_I_vector_plane_strain))
    return nothing
end
#---------------------------------------------------------------------------------------------
# 4階の恒等テンソルのベクトル形式
#---------------------------------------------------------------------------------------------
function make_Is_vector_plane_strain()::Matrix{Float64}
    FOID = zeros(Float64, 4, 4)
    FOID[1, 1] = 1.0
    FOID[2, 2] = 1.0
    FOID[3, 3] = 0.5
    FOID[4, 4] = 1.0
    return FOID
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(make_Is_vector_plane_strain))
    return nothing
end
#---------------------------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------------------------
function make_deviatoric_projection_tensor_plane_strain()::Matrix{Float64}
    # make base vector
    SOID = make_I_vector_plane_strain()
    # make base matrix
    FOID = make_Is_vector_plane_strain()
    # Set deviatoric projection tensor
    R1D3 = 1.0 / 3.0
    DEVPRJ = Matrix{Float64}(undef, 3, 3)
    for i in 1 : 3
        for j in 1 : 3
            DEVPRJ[i, j] = FOID[i, j] - SOID[i] * SOID[j] * R1D3
        end
    end
    #
    return DEVPRJ
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(make_deviatoric_projection_tensor_plane_strain))
    return nothing
end