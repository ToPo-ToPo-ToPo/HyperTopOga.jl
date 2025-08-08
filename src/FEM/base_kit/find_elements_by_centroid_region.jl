#要素の中心座標を元に, 領域に入る要素番号を出力

function find_elements_by_centroid_region(elements, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64)
    selected = Int[]
    for (i, elem) in enumerate(elements)
        coords = elem.shape.coord_matrix  # [nnode × 2] matrix
        centroid = sum(eachrow(coords)) ./ size(coords, 1)
        x, y = centroid
        if xmin ≤ x ≤ xmax && ymin ≤ y ≤ ymax
            push!(selected, i)  # 要素番号（1始まり）
        end
    end
    return selected
end