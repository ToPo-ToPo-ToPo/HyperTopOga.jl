
#----------------------------------------------------------------------
#
#----------------------------------------------------------------------
module HyperTopOga

#----------------------------------------------------
# Packages
#----------------------------------------------------
# 線形代数関連
using LinearAlgebra
using SparseArrays
using LinearSolve
# 自動微分関連   
using Enzyme
import .EnzymeRules: augmented_primal, reverse
using .EnzymeRules 
using Checkpointing
# 並列計算
using FLoops
# その他
using NearestNeighbors
using Printf
using CPUTime
using WriteVTK

#----------------------------------------------------
# 自作ライブラリ
#----------------------------------------------------
# Dense and Sparse matrix 
include("my_linear_algebra/header.jl")

# 他の物理モデルのモジュールと共有する部分
include("FEM/base_kit/header.jl")

# 物理モデル
include("FEM/physics/header.jl")

# 最適化関連
include("optimization/header.jl")

end # module HyperTopOga
