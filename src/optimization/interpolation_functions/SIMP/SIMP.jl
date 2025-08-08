
#----------------------------------------------------------------
# SIMP法
#----------------------------------------------------------------
mutable struct SIMP <: AbstractInterpolation
    params::Vector{Float64}
    penalty::Float64
end
function compute(::SIMP, ::Vector{Float64})::Float64
    return 0.0
end
function compute_derivative(::SIMP, ::Vector{Float64}, ::Int64)::Float64
    return 0.0
end