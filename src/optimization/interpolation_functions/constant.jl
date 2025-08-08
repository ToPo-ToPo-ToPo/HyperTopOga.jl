
#----------------------------------------------------------------
# 定数
#----------------------------------------------------------------
mutable struct Constant <: AbstractInterpolation 
    param::Float64
end
#----------------------------------------------------------------
# 定数
#----------------------------------------------------------------
function compute(method::Constant, ::Vector{Float64})::Float64
    return method.param
end
#----------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------
function EnzymeRules.inactive(::typeof(compute), method::Constant, ::Vector{Float64})
    return nothing
end
#----------------------------------------------------------------
#
#----------------------------------------------------------------
function compute_derivative(::Constant, ::Vector{Float64}, ::Int64)::Float64
    return 0.0
end