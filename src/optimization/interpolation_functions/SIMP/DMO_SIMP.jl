
#----------------------------------------------------------------
# DMO SIMP
#----------------------------------------------------------------
mutable struct DMOSIMP <: AbstractInterpolation
    params::Vector{Float64}
    penalty::Float64
end
#----------------------------------------------------------------
# DMO SIMP
#----------------------------------------------------------------
function compute(method::DMOSIMP, x::Vector{Float64})::Float64
    # check
    if length(x) + 1 != length(method.params)
        error("DMOSIMP: size of x + 1 != size of params")
    end
    # main
    base = method.params[1]
    temp = base
    for i in eachindex(x)
        param_i = method.params[i+1]
        dens_i = x[i]^method.penalty
        dens_j = 1.0
        for j in eachindex(x)
            dens_j *= inner_calc_dmo_simp(i, j, x[j], method.penalty)
        end
        temp += dens_i * dens_j * (param_i - base)
    end
    return temp
end
#----------------------------------------------------------------
# DMO SIMP sub
#----------------------------------------------------------------
function inner_calc_dmo_simp(
    i::Int64, 
    j::Int64, 
    x_j::Float64, 
    penalty::Float64
    )::Float64
    
    val = 1.0
    if i != j
        val = 1.0 - x_j^penalty
    end
    return val
end
#----------------------------------------------------------------
#
#----------------------------------------------------------------
function compute_derivative(method::DMOSIMP, x::Vector{Float64}, id::Int64)::Float64

    #
    dens_i = x[id]
    params = method.params
    p = method.penalty
    base = method.params[1]

    #
    temp = 0.0
    for j in eachindex(x)
        w_j = 1.0
        if id == j
            w_j *= p* dens_i^(p-1.0)
            for k in eachindex(x)
                if id == k
                    continue
                end
                w_j *= (1.0 - x[k]^p)
            end
        else
            w_j *= - p * dens_i^(p-1.0) * x[j]^p
            for k in eachindex(x)
                if id == k || j == k
                    continue
                end
                w_j *= (1.0 - x[k]^p)
            end
        end
        #
        temp += w_j * (params[j+1]-base)
    end
    #
    return temp
end