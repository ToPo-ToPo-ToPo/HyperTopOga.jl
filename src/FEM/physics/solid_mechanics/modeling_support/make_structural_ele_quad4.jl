
#------------------------------------------------------------------------------------------------------------------------
# ボクセル要素向け（Quad4）
#------------------------------------------------------------------------------------------------------------------------
function make_structural_ele_quad4(
    connectivities::Matrix{Node}, 
    material::MD, 
    element_type::ET, 
    thickness::Float64
    )::Vector{StructuralElement} where {MD<:AbstractMaterialData, ET<:AbstractElementType}

    #
    elements = Vector{StructuralElement}(undef, size(connectivities)[1])
    #
    for e in eachindex(elements)
        elements[e] = StructuralQuad4(e, connectivities[e, :], element_type, material, thickness)
    end
    #
    return elements
end