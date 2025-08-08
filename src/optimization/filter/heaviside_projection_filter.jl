

#-------------------------------------------------------------------------------------------------
# トポロジー最適化における密度フィルター
#-------------------------------------------------------------------------------------------------
mutable struct HeavisideProjectionFilter <: AbstractOptFilter
    beta::Float64
    eta::Float64
    x_tilde::Vector{Float64}
    num_svalue_types::Int64            # 設計変数の種類
    num_svalue_domain::Int64           # 設計変数が割り当てられる領域の数   
    dataset::OptFilterDataset          # フィルターに必要なデータのコンテナ
    #----------------------------------------------------------------------------
    #
    #----------------------------------------------------------------------------
    function HeavisideProjectionFilter(beta::Float64=5.0, eta::Float64=0.5)
        return new(beta, eta, Vector{Float64}(), 0, 0)
    end
end
#-------------------------------------------------------------------------------------------------
# データの更新
#-------------------------------------------------------------------------------------------------
function filter_setup!(
    filter::HeavisideProjectionFilter, 
    num_svalue_types::Int64, 
    num_svalue_domain::Int64, 
    physics::P, 
    radius::Float64, 
    design_space_type::D=ElementBase()
    )::Nothing where {P<:AbstractPhysics, D<:AbstractDesignSpace}

    # フィルタリングに必要な情報を作成
    dataset = setup_opt_filter(physics.nodes, physics.elements, radius, design_space_type)
    # パラメータの設定
    filter.num_svalue_types = num_svalue_types
    filter.num_svalue_domain = num_svalue_domain
    filter.dataset = dataset
    filter.x_tilde = Vector{Float64}(undef, num_svalue_domain * num_svalue_types)
    #
    return nothing
end
#-------------------------------------------------------------------------------------------------
# フィルタリングされた設計変数を作成（自動微分用）
#-------------------------------------------------------------------------------------------------
function filtering(filter::HeavisideProjectionFilter, x::Vector{Float64})::Vector{Float64}
    
    # 初期化
    x_tilde2 = zero(x)
    beta = filter.beta
    eta = filter.eta
    
    # domain loop
    Threads.@threads for e in 1 : filter.num_svalue_domain
        # set dataset
        id_list = filter.dataset.id_list[e, :]
        weights = filter.dataset.weights[e, :]
        #
        for itype in 1 : filter.num_svalue_types
            # zero clear
            local A = 0.0
            local B = 0.0
            # フィルター半径に含まれる要素のループ
            for i in 1 : filter.dataset.num_ids[e]
                #
                ID_i = id_list[i]
                # 設計変数の種類itypeと定義の位置id_jの場所を設定
                index_i = make_index(filter.num_svalue_types, itype, ID_i)
                # 重みを取得
                weight_i = weights[i]
                # 和を計算
                A += weight_i * x[index_i]
                B += weight_i
            end
            # フィルタリングされた密度を計算
            dens = A / B
            # 設計変数の種類itypeと定義の位置idの場所を設定
            index_e = make_index(filter.num_svalue_types, itype, e)
            # データの更新
            x_tilde2[index_e] = compute_heviside_function(dens, beta, eta)
        end
    end
    #
    return x_tilde2
end
#-------------------------------------------------------------------------------------------------
# フィルタリングされた設計変数を作成（自動微分用）
#-------------------------------------------------------------------------------------------------
function density_filtering(filter::HeavisideProjectionFilter, x::Vector{Float64})::Vector{Float64}
    # 初期化
    x_tilde2 = zero(x)
    # domain loop
    Threads.@threads for e in 1 : filter.num_svalue_domain
        # set dataset
        id_list = filter.dataset.id_list[e, :]
        weights = filter.dataset.weights[e, :]
        #
        for itype in 1 : filter.num_svalue_types
            # zero clear
            local A = 0.0
            local B = 0.0
            # フィルター半径に含まれる要素のループ
            for i in 1 : filter.dataset.num_ids[e]
                #
                ID_i = id_list[i]
                # 設計変数の種類itypeと定義の位置id_jの場所を設定
                index_i = make_index(filter.num_svalue_types, itype, ID_i)
                # 重みを取得
                weight_i = weights[i]
                # 和を計算
                A += weight_i * x[index_i]
                B += weight_i
            end
            # フィルタリングされた密度を計算
            dens = A / B
            # 設計変数の種類itypeと定義の位置idの場所を設定
            index_e = make_index(filter.num_svalue_types, itype, e)
            # データの更新
            x_tilde2[index_e] = dens
        end
    end
    #
    return x_tilde2
end
#-------------------------------------------------------------------------------------------------
# 感度解析 (随伴変数法用)
#-------------------------------------------------------------------------------------------------
function filtering_df(filter::HeavisideProjectionFilter, df_org::Vector{Float64}, x::Vector{Float64})::Vector{Float64}
    # stock tmp value
    df = zero(df_org)
    beta = filter.beta
    eta = filter.eta
    coeff = tanh(beta * eta) + tanh(beta * (1.0 - eta))
    # 密度フィルター
    x_tilde = density_filtering(filter, x)
    # domain loop
    Threads.@threads for e in 1 : filter.num_svalue_domain
        # set dataset
        id_list_e = filter.dataset.id_list[e, :]
        weights_e = filter.dataset.weights[e, :]
        for itype in 1 : filter.num_svalue_types
            #
            local df_e = 0.0
            # フィルター半径に含まれる要素のループ
            for i in 1 : filter.dataset.num_ids[e]
                # ID eに含まれるIDと重みを取得
                ID_i = id_list_e[i]
                weights_i = filter.dataset.weights[ID_i, :]
                # フィルター半径に含まれる要素のループ
                temp = 0.0
                for j in 1 : filter.dataset.num_ids[ID_i]
                    temp += weights_i[j]
                end
                # 設計変数の種類itypeと定義の位置id_jの場所を設定
                index_i = make_index(filter.num_svalue_types, itype, ID_i)
                # 密度の感度
                ddens_ds = weights_e[i] / temp
                ddens_d_ddens = (beta*(1.0-tanh(beta*(x_tilde[index_i]-eta))^2.0)) / coeff
                # 感度の更新
                df_e += df_org[index_i] * ddens_ds * ddens_d_ddens
            end
            # 設計変数の種類itypeと定義の位置idの場所を設定
            index_e = make_index(filter.num_svalue_types, itype, e)
            df[index_e] = df_e
        end
    end
    return df
end
#-------------------------------------------------------------------------------------------------
# フィルタリングされた感度を作成
#-------------------------------------------------------------------------------------------------
function filtering_grad(::HeavisideProjectionFilter, df_org::Vector{Float64}, ::Vector{Float64})::Vector{Float64}
    return df_org
end
#-------------------------------------------------------------------------------------------------
# update
#-------------------------------------------------------------------------------------------------
function update_status!(filter::HeavisideProjectionFilter, opt_step::Int64)::Nothing
    #beta_new = 20.0 * (Float64(opt_step) / 100.0) 
    #filter.beta = beta_new
    return nothing
end
#-------------------------------------------------------------------------------------------------
# 更新関数
#-------------------------------------------------------------------------------------------------
function compute_heviside_function(dens::Float64, beta::Float64, eta::Float64)::Float64
    return (tanh(beta*eta) + tanh(beta*(dens - eta))) / (tanh(beta * eta) + tanh(beta*(1.0 - eta)))
end
