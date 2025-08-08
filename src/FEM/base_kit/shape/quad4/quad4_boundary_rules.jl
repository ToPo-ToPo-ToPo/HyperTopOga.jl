

#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(make_boundary_detJ), ::Quad4, boundary_id::Int64, J::Matrix{Float64})
    return nothing
end