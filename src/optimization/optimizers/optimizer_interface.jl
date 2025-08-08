
#--------------------------------------------------------------------------
# 
#--------------------------------------------------------------------------
abstract type AbstractOptimizer end
#--------------------------------------------------------------------------
# 
#--------------------------------------------------------------------------
mutable struct MyGCMMA <: AbstractOptimizer
    model
    asyinit::Float64 
    asyinc::Float64 
    asydec::Float64
    move::Float64
    move_limit::Float64
    innerit_max::Int64 
    function MyGCMMA(asyinit::Float64=0.5, asyinc::Float64=1.2, asydec::Float64=0.7, move::Float64=0.5, move_limit::Float64=0.1, innerit_max::Int64=30)
        return new(nothing, asyinit, asyinc, asydec, move, move_limit, innerit_max)
    end
end

#--------------------------------------------------------------------------
# 最適化アルゴリズム全体のインターフェースを定義
#--------------------------------------------------------------------------
function set_optimizer!(optimizer::MyGCMMA, num_x::Int64)
    optimizer.model = GCMMA(num_x, optimizer.asyinit, optimizer.asyinc, optimizer.asydec, optimizer.move, optimizer.move_limit, optimizer.innerit_max)
end
#--------------------------------------------------------------------------
# 最適化アルゴリズム全体のインターフェースを定義
#--------------------------------------------------------------------------
function set_option!(optimizer::MyGCMMA, max_eval::Int64)
    optimizer.model.max_eval = max_eval
end
#--------------------------------------------------------------------------
# 
#--------------------------------------------------------------------------
function optimize!(optimizer::MyGCMMA, physics::AbstractPhysics, opt_settings::AbstractOptSettings, opt_filter::AbstractOptFilter, x0::Vector{Float64})
    return my_optimize(optimizer.model, physics, opt_settings, opt_filter, x0)
end
