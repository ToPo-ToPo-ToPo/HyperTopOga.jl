
#--------------------------------------------------------------------------------------------------------
# MMAに関するデータをまとめた構造体
#--------------------------------------------------------------------------------------------------------
mutable struct GCMMA <: AbstractMMA
    n::Int64                                     # 設計変数の総数
    m::Int64                                     # 制約条件の数
    max_eval::Int64                              # 最大ステップ
    innerit_max::Int64                           # 内部ループの最大ステップ
    lower_bounds::Vector{Float64}                # 設計変数の下限制約値
    upper_bounds::Vector{Float64}                # 設計変数の上限制約値
    asyinit::Float64
    asyinc::Float64
    asydec::Float64
    move::Float64                            
    move_limit::Float64                             # ムーブリミット
    kkttol::Float64                                 # 収束判定値
    xtol_rel::Float64                               # 設計変数の変化量に対する収束判定値
    ftol_rel::Float64                               # 目的関数の変化量に対する収束判定値
    objective::Vector{AbstractEvaluationFunction}     # 目的関数
    objective0::Vector{Float64}                     # 目的関数の初期値
    objective_weight::Vector{Function}              # 目的関数の重み
    objective_df_methods::Vector{AbstractDfMethod}  # 目的関数の感度解析手法
    constraint::Vector{AbstractEvaluationFunction}    # 制約関数
    constraint_lim::Vector{Function}                # 制約条件の制約値
    constraint_df_methods::Vector{AbstractDfMethod} # 制約関数の感度解析手法
    #---------------------------------------------------------------------
    # MMAに関するデータをまとめた構造体
    # 具体的な関数は別で定義
    #---------------------------------------------------------------------
    function GCMMA(
        n::Int64, asyinit::Float64, asyinc::Float64, asydec::Float64, move::Float64, move_limit::Float64, innerit_max::Int64, 
        max_eval::Int64=200, kkttol=1.0e-06, xtol_rel=0.0, ftol_rel=0.0
        )
        return new(
            n, 0, max_eval, innerit_max, Vector{Float64}(), Vector{Float64}(), 
            asyinit, asyinc, asydec, move, move_limit, 
            kkttol, xtol_rel, ftol_rel,
            Vector{AbstractEvaluationFunction}(), Vector{Float64}(), Vector{Function}(), Vector{AbstractDfMethod}(),
            Vector{AbstractEvaluationFunction}(), Vector{Function}(), Vector{AbstractDfMethod}()
        )
    end
end
#--------------------------------------------------------------------------------------------------------
# 目的関数、制約関数の計算+それぞれの感度の計算を行う
#--------------------------------------------------------------------------------------------------------
function compute_inner_loop(
    opt::GCMMA, 
    physics::P, 
    opt_settings::OptimizationSettings,
    opt_filter::OF, 
    opt_step::Int64,
    x::Vector{Float64}
    )::Tuple{Float64, Vector{Float64}} where {P<:AbstractPhysics, OF<:AbstractOptFilter}

    # 目的関数の計算
    f0val = 0.0
    for i in eachindex(opt.objective)
        obj = compute_f(x, opt.objective[i], physics, opt_settings, opt_filter, opt.objective_df_methods[i])
        f0val += opt.objective_weight[i](opt_step) * obj / opt.objective0[i]
    end
    # 制約関数の計算
    fval = Vector{Float64}(undef, length(opt.constraint))
    for i in eachindex(opt.constraint)
        fval[i] = compute_f(x, opt.constraint[i], physics, opt_settings, opt_filter, opt.constraint_df_methods[i]) / opt.constraint_lim[i](opt_step) - 1.0
    end
    return f0val, fval
end
#--------------------------------------------------------------------------------------------------------
# MMAの最適化計算を実行
#--------------------------------------------------------------------------------------------------------
function my_optimize(
    opt::GCMMA, 
    physics::P, 
    opt_settings::OptimizationSettings, 
    opt_filter::OF, 
    xval::Vector{Float64}
    )::Tuple{Vector{MmaHistory}, Vector{Float64}} where {P<:AbstractPhysics, OF<:AbstractOptFilter}

    #
    m = opt.m
    n = opt.n
    epsimin = 1.0e-07
    xold1   = copy(xval)
    xold2   = copy(xval)
    xmin    = opt.lower_bounds
    xmax    = opt.upper_bounds
    low     = Vector{Float64}(undef, n)
    upp     = Vector{Float64}(undef, n)
    c       = fill(1.0e+05, m)
    d       = fill(1.0, m)
    a0      = 1.0
    a       = zeros(Float64, m)
    # for gcmma param--------------
    raa0    = 0.01
    raa     = fill(0.01, m)
    raa0eps = 1.0e-06
    raaeps  = fill(raa0eps, m)
    #------------------------------
    outeriter = 0
    maxoutit  = opt.max_eval
    kkttol  = opt.kkttol
    history = Vector{MmaHistory}(undef, opt.max_eval+1)
    for i in eachindex(history)
        history[i] = MmaHistory(length(opt.objective), length(opt.constraint))
    end

    xmma = Vector{Float64}(undef, n)
    ymma = Vector{Float64}(undef, m)
    zmma = 0.0
    lam = Vector{Float64}(undef, m)
    xsi = Vector{Float64}(undef, n)
    eta = Vector{Float64}(undef, n)
    mu = Vector{Float64}(undef, m)
    zet = 0.0

    # パラメータの設定
    asyinit = opt.asyinit
    asyinc = opt.asyinc
    asydec = opt.asydec
    move = opt.move

    # Outputs status in progress
    println("=====================================================================================")
    println("                                    Basic GCMMA                                      ")
    println("=====================================================================================")

    # 初期値の計算と更新
    f0vals = initial_compute(opt, physics, opt_settings, opt_filter, xval)
    opt.objective0 = abs.(f0vals)

    # 再計算
    f0val, f0vals, df0dx, fval, dfdx = compute(opt, physics, opt_settings, opt_filter, outeriter, xval)

    # Outputs status in progress
    output_progress_info!(opt, outeriter, physics, opt_settings, opt_filter, xval, f0val, f0vals, df0dx, fval, 0, history)

    # The outer iterations start:
    kktnorm = kkttol+10
    while (kktnorm > kkttol) & (outeriter < maxoutit)
        # Update optimization steps
        outeriter = outeriter+1
        opt_settings.opt_steps = outeriter

        # reset output output flag
        opt_settings.output_flag = false
        if outeriter % 10 == 0
            opt_settings.output_flag = true
        end

        # update
        xmin = max.(opt.lower_bounds, xval .- opt.move_limit)
        xmax = min.(opt.upper_bounds, xval .+ opt.move_limit)
        
        # The parameters low, upp, raa0 and raa are calculated:
        low, upp, raa0, raa = asymp(outeriter,n,xval,xold1,xold2,xmin,xmax,low,upp,raa0,raa,raa0eps,raaeps,df0dx,dfdx,asyinit,asyinc,asydec)
    
        # The GCMMA subproblem is solved at the point xval:
        xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp = gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp,raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d,move)

        # The user should now calculate function values (no gradients)
        # of the objective- and constraint functions at the point xmma
        # ( = the optimal solution of the subproblem).
        # The results should be put in f0valnew and fvalnew.
        f0valnew, fvalnew = compute_inner_loop(opt, physics, opt_settings, opt_filter, outeriter, xmma)

        # It is checked if the approximations are conservative:
        conserv = concheck(epsimin,f0app,f0valnew,fapp,fvalnew)
    
        # While the approximations are non-conservative (conserv=0),
        # repeated inner iterations are made:
        innerit=0
        if conserv == 0
            while (conserv == 0) & (innerit <= opt.innerit_max)
                innerit = innerit+1;
                
                # New values on the parameters raa0 and raa are calculated:
                raa0, raa = raaupdate(xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew,f0app,fapp,raa0,raa,raa0eps,raaeps,epsimin)
                
                # The GCMMA subproblem is solved with these new raa0 and raa:
                xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp = gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp,raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d,move)

                # The user should now calculate function values (no gradients)
                # of the objective- and constraint functions at the point xmma
                # ( = the optimal solution of the subproblem).
                # The results should be put in f0valnew and fvalnew:
                f0valnew, fvalnew = compute_inner_loop(opt, physics, opt_settings, opt_filter, outeriter, xmma)
    
                # It is checked if the approximations have become conservative:
                conserv = concheck(epsimin,f0app,f0valnew,fapp,fvalnew)
            end
        end

        # No more inner iterations. Some vectors are updated:
        xold2 .= xold1
        xold1 .= xval
        xval  .= xmma
        
        # The user should now calculate function values and gradients
        # of the objective- and constraint functions at xval.
        # The results should be put in f0val, df0dx, fval and dfdx:
        f0val, f0vals, df0dx, fval, dfdx = compute(opt, physics, opt_settings, opt_filter, outeriter, xval)

        # The residual vector of the KKT conditions is calculated:
        residu, kktnorm, residumax = kktcheck(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,fval,dfdx,a0,a,c,d)

        # Outputs status in progress
        output_progress_info!(opt, outeriter, physics, opt_settings, opt_filter, xval, f0val, f0vals, df0dx, fval, innerit, history)

        # Update filter val
        update_filter_status!(opt, opt_filter, outeriter)
    end

    return history, xval
end
#--------------------------------------------------------------------------------------------------------
# 途中の状態を出力
#--------------------------------------------------------------------------------------------------------
function output_progress_info!(
    ::GCMMA, 
    step::Int64, 
    physics::P, 
    opt_settings::AO, 
    opt_filter::OF, 
    xval::Vector{Float64}, 
    f0val::Float64, 
    f0vals::Vector{Float64}, 
    df0dx::Vector{Float64}, 
    fval::Vector{Float64}, 
    innerit::Int64, 
    history::Vector{MmaHistory}
    ) where {P<:AbstractPhysics, AO<:AbstractOptSettings, OF<:AbstractOptFilter}

    println("--------------------------------------------------------------------------------------")
    @printf("Opt step:    %-3d\n", step)
    @printf("Inner iter:  %-2d\n", innerit)
    @printf("Objective:   %-.6e\n", f0val)
    for iobj in eachindex(f0vals)
        @printf("Objective%1d:  %-.6e\n", iobj, f0vals[iobj])
    end
    for icons in eachindex(fval)
        @printf("Constraint%1d: %-.6e <= 0.0\n", icons, fval[icons])
    end

    # vtkファイルに結果を書き出し
    out_name = physics.output_file_name * "_gcmma_" * string(step)
    output_vtu_for_mma_gcmma(physics, opt_settings, opt_filter, xval, df0dx, out_name, opt_settings.design_space_type)
        
    # 途中の関数値の保存
    add_history_data!(history, step, f0vals, fval)

end