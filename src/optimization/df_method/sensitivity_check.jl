
#------------------------------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------------------------------
function sensitivity_check(s::Vector{Float64}, df1::Vector{Float64}, df2::Vector{Float64}, df3::Vector{Float64})

    #
    max_abs_value = maximum(abs.(df1))
    
    #
    println("数値微分の感度の最大値: ", max_abs_value)
    for i in eachindex(s)
        error1 = 100.0 * (df1[i]-df2[i]) / max_abs_value
        error2 = 100.0 * (df1[i]-df3[i]) / max_abs_value
        @printf("x = %.5e, df1 = %.5e, df2 = %.5e, df3 = %.5e, df1 vs df2 = %.5e, df1 vs df3 = %.5e\n", s[i], df1[i], df2[i], df3[i], error1, error2)
    end

    #
    println("norm of error (df1 vs df2): $(100.0 * norm(df1-df2) / max_abs_value) [%]")
    println("norm of error (df1 vs df3): $(100.0 * norm(df1-df3) / max_abs_value) [%]")

    #
    return nothing

end
#------------------------------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------------------------------
function sensitivity_check(s::Vector{Float64}, df1::Vector{Float64}, df2::Vector{Float64})
    #
    max_abs_value = maximum(abs.(df1))
    
    #
    println("数値微分の感度の最大値: ", max_abs_value)
    for i in eachindex(s)
        error = 100.0 * (df1[i]-df2[i]) / max_abs_value
        @printf("x = %.5e, df1 = %.5e, df2 = %.5e, df1 vs df2 = %.5e\n", s[i], df1[i], df2[i], error)
    end

    #
    println("norm of error (df1 vs df2): $(100.0 * norm(df1-df2) / max_abs_value) [%]")

    #
    return nothing

end