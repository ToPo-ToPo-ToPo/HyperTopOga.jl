
#----------------------------------------------------------------
# extended SIMP
#----------------------------------------------------------------
mutable struct ExtendedSIMP <: AbstractInterpolation
    params::Vector{Float64}
    base::Float64
    penalty::Float64
end
function compute(::ExtendedSIMP, ::Vector{Float64})
    return 0.0
end