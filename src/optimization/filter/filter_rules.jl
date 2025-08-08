

#----------------------------------------------------------------------------
# カスタムルールの定義
#----------------------------------------------------------------------------
function EnzymeRules.inactive(::typeof(set_filtering_radius), elements::Vector{E}, radius_coeff::Float64) where {E<:AbstractElement}
    return nothing
end