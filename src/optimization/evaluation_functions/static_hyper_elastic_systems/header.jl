

#----------------------------------------------------------------------------
# 抽象型の定義
#----------------------------------------------------------------------------
abstract type AbstractHyperElasticEvalFunc <: AbstractEvaluationFunction end

#----------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------
include("output_vtu.jl")
include("nonlinear_fem.jl")
#
include("static_end_compliance.jl")
include("static_pnorm_mises_stress.jl")