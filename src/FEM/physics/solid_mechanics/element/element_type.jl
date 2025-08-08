
#
abstract type AbstractElementType end

#----------------------------------------------------------------------------------
# トータルラグランジュに基づく有限変形
#----------------------------------------------------------------------------------
# 抽象型
abstract type AbstractTotalLagrange <: AbstractElementType end

# 通常の解析用
struct TotalLagrange <: AbstractTotalLagrange end