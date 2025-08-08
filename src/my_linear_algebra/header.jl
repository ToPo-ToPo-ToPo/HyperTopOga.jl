
#-----------------------------------------------------------------------------------
# Dense matirx library
#-----------------------------------------------------------------------------------
include("my_linear_algebra/MyLinearAlgebra.jl")
include("my_linear_algebra/tensor_4th_to_matrix_2th.jl")

#-----------------------------------------------------------------------------------
# Sparse matrix library for CSC format 
#-----------------------------------------------------------------------------------
include("sparse_matrix_csc/MySparseMatrixCSC.jl")
include("sparse_matrix_csc/linear_solvers/UMFPACK.jl")
include("linear_solver.jl")

