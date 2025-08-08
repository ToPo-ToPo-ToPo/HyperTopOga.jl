
#---------------------------------------------------------------------------------------
# 重複データの削除とソートを行う
#---------------------------------------------------------------------------------------
function delete_duplicate_and_sort(datasets::Vector{Data2Data}) 
    for i in 1 : length(datasets)
        # ソート
        sort!(datasets[i].data)
        # 重複する要素を削除
        datasets[i].data = unique(datasets[i].data)
    end
    return datasets
end