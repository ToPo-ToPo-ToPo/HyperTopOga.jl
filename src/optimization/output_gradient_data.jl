
#--------------------------------------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------------------------------------
function output_gradient_data(file_name::String, method_name::String, num_svalue_types::Int64, num_domain::Int64, df::Vector{Float64})

    for isvalue in 1 : num_svalue_types
        # ファイルを開く
        out_file = open(file_name * "_" * method_name * "_gradient_for_s" * string(isvalue) * ".dat", "w")
        # データを書き込む
        for e in 1 : num_domain
            # 設計変数位置を取得
            index = make_index(num_svalue_types, isvalue, e)
            # 要素番号と数値微分の結果の出力
            println(out_file, "$(e) $(df[index])")
        end
        close(out_file)
    end
end