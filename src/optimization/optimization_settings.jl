
#----------------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------------
abstract type AbstractOptSettings end
#----------------------------------------------------------------------------------------------
# 最適化に必要なパラメータをまとめたコンテナ
# output_flagはデフォルトはfalseとなっており、結果の出力はしない
#----------------------------------------------------------------------------------------------
mutable struct OptimizationSettings{IP<:AbstractInterpolation, DS<:AbstractDesignSpace} <: AbstractOptSettings
    num_svalue_types::Int64                   # 設計変数の種類
    interpolate::IP                           # 結果出力時の内挿関数
    opt_steps::Int64                          # 最適化のステップ番号
    design_space_type::DS                     # 要素ベースまたは節点ベースか
    output_flag::Bool                         # 最適化計算中にvtkファイルを出力するかどうかのフラグ
    #-----------------------------------------------------------------------
    # コンストラクタ
    #-----------------------------------------------------------------------
    function OptimizationSettings(num_svalue_types::Int64, interpolate::IP, design_space_type::DS, output_flag::Bool
        ) where {IP<:AbstractInterpolation, DS<:AbstractDesignSpace}
        return new{IP, DS}(num_svalue_types, interpolate, 0, design_space_type, output_flag)
    end
    #-----------------------------------------------------------------------
    # コンストラクタ
    #-----------------------------------------------------------------------
    function OptimizationSettings(num_svalue_types::Int64, interpolate::IP, design_space_type::DS
        ) where {IP<:AbstractInterpolation, DS<:AbstractDesignSpace}
        return new{IP, DS}(num_svalue_types, interpolate, 0, design_space_type, false)
    end
    #-----------------------------------------------------------------------
    # コンストラクタ
    #-----------------------------------------------------------------------
    function OptimizationSettings(num_svalue_types::Int64, interpolate::IP, output_flag::Bool
        ) where {IP<:AbstractInterpolation}
        return new{IP, ElementBase}(num_svalue_types, interpolate, 0, ElementBase(), output_flag)
    end
    #-----------------------------------------------------------------------
    # コンストラクタ
    #-----------------------------------------------------------------------
    function OptimizationSettings(num_svalue_types::Int64, interpolate::IP
        ) where {IP<:AbstractInterpolation}
        return new{IP, ElementBase}(num_svalue_types, interpolate, 0, ElementBase(), false)
    end
end