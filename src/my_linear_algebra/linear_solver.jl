
abstract type AbstractLinearSolver end
struct DirectSolverUMFPACK <: AbstractLinearSolver end
#------------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------------
function linear_solve(::DirectSolverUMFPACK, A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, x0::Vector{Float64}=zero(b))::Vector{Float64}
    return my_linear_solve_UMFPACK(A, b)
end