
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(make_N), ::StructuralQuad4, coordinate::Vector{Float64})
    return nothing
end
#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(make_B), ::StructuralQuad4, dNdx::Matrix{Float64})
    return nothing
end