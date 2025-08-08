
#---------------------------------------------------------------------------------------
# 要素のマルチカラーの情報を保存する構造体
#---------------------------------------------------------------------------------------
struct MultiColorizedElementDataset
    element_ids::Vector{Int64}                   # 1つの色に含まれる要素番号のリスト
    function MultiColorizedElementDataset()
        return new(Vector{Int64}())
    end
end
#---------------------------------------------------------------------------------------
# データを取得
#---------------------------------------------------------------------------------------
function create_multi_colored_elem_ids(nodes::Vector{Node}, elements::Vector)::Vector{MultiColorizedElementDataset}

    # create node-element table
    node2elem = Vector{Data2Data}(undef, length(nodes))
    for i in 1 : length(node2elem)
        node2elem[i] = Data2Data()
    end
    for (e, element) in enumerate(elements)
        for node_id in element.shape.connects
            push!(node2elem[node_id].data, e)
        end
    end

    # create elemet-node table
    elem2node = Vector{Data2Data}(undef, length(elements))
    for i in 1 : length(elem2node)
        elem2node[i] = Data2Data()
    end
    for i in 1 : length(nodes)
        for elem_id in node2elem[i].data
            push!(elem2node[elem_id].data, i)
        end
    end

    # create element-element table (graph)
    elem2elem = Vector{Data2Data}(undef, length(elements))
    for i in 1 : length(elem2elem)
        elem2elem[i] = Data2Data()
    end
    for e in 1 : length(elements)
        is_registered = fill(false, length(elements))
        for node_id in elem2node[e].data
            for elem_id in node2elem[node_id].data
                if e == elem_id
                    continue
                end
                if is_registered[elem_id] == true
                    continue
                end
                push!(elem2elem[e].data, elem_id)
                is_registered[elem_id] = true
            end
        end
    end

    # sort by number of nodes belong to the element
    deg_size = [(length(elem2elem[i].data), i) for i in 1 : length(elements)]
    sort!(deg_size, rev=true)

    # Welsh-Powell Algorithm
    max_color = 0
    color_elem = fill(-1, length(elements))
    for i in 1 : size(deg_size)[1]
        elem_id1 = deg_size[i][2]
        color = 0
        is_used = true

        while is_used
            is_used = false
            for elem_id2 in elem2elem[elem_id1].data
                if elem_id1 == elem_id2
                    continue
                end
                if color_elem[elem_id2] == color
                    is_used = true
                    color += 1
                    break
                end
            end
        end
        color_elem[elem_id1] = color
        max_color = max(max_color, color)
    end

    # set element group by color
    colored_elem_IDs = Vector{MultiColorizedElementDataset}(undef, max_color + 1)
    for i in 1 : max_color + 1
        colored_elem_IDs[i] = MultiColorizedElementDataset()
    end
    for (e, color) in enumerate(color_elem)
        push!(colored_elem_IDs[color+1].element_ids, e)
    end
    
    #
    return colored_elem_IDs
end