
#
mutable struct MmaHistory
    opt_step::Int64
    obj_func::Vector{Float64}
    con_func::Vector{Float64}
    function MmaHistory(num_obj::Int64, num_con::Int64)
        return new(0, Vector{Float64}(undef, num_obj), Vector{Float64}(undef, num_con)) 
    end
end

# 途中の関数値の保存
function add_history_data!(history::Vector{MmaHistory}, step::Int64, f0vals::Vector{Float64}, fval::Vector{Float64})
    history[step+1].opt_step = step
    for iobj in eachindex(f0vals)
        history[step+1].obj_func[iobj] = f0vals[iobj]
    end
    for icon in eachindex(fval)
        history[step+1].con_func[icon] = fval[icon]
    end
end

# 最適化計算結果を出力
function output(history::Vector{MmaHistory}, output_file_name::String, num_obj::Int64, num_con::Int64)

    # 目的関数履歴を出力
    for iobj in 1 : num_obj
        out_file = open(output_file_name * "_obj_" * string(iobj) * ".dat", "w")
        for (step, h) in enumerate(history)
            println(out_file, "$(step) $(h.obj_func[iobj])")
        end
        close(out_file)
    end

    # 制約関数の履歴を出力
    for icon in 1 : num_con
        out_file = open(output_file_name * "_const_" * string(icon) * ".dat", "w")
        for (step, h) in enumerate(history)
            println(out_file, "$(step) $(h.con_func[icon])")
        end
        close(out_file)
    end

end