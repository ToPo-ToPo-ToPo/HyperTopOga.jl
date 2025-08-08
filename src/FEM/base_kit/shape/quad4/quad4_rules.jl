
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(make_N_components), ::Quad4, coordinate::Vector{Float64})
    return nothing
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(make_dNdxi), ::Quad4, coordinate::Vector{Float64})
    return nothing
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(make_J), shape::Quad4, dNdxi::Matrix{Float64})
    return nothing
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(make_dNdx), ::Quad4, dNdxi::Matrix{Float64}, J::Matrix{Float64})
    return nothing
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(integral_factor), shape::Quad4, ip::Int64)
    return nothing
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(boundary_integral_factor), shape::Quad4, boundary_id::Int64, ip::Int64)
    return nothing
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(circumscribed_circle_radius), shape::Quad4)
    return nothing
end