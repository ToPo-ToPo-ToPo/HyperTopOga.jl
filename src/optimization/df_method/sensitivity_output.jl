
function output_sensitivity_result(num_svalue_types::Int64, num_domain::Int64, output_file_name::String, df_numeric::Vector{Float64}, df_auto::Vector{Float64}, df_adjoint::Vector{Float64})
    # テキストファイルにデータを出力する
    for isvalue in 1 : num_svalue_types
        # ファイルを開く
        out_file = open(output_file_name * "_gradient_for_s" * string(isvalue) * ".dat", "w")
        # データを書き込む
        println(out_file, "df_numeric df_auto df_adjoint")
        for e in 1 : num_domain
            # 設計変数位置を取得
            index = Problem.make_index(num_svalue_types, isvalue, e)
            println(out_file, "$(df_numeric[index]) $(df_auto[index]) $(df_adjoint[index])")
        end
        close(out_file)
    end
end