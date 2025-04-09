function output_est_to_CSV(GP::global_param, dfB::DataFrame, dfC::DataFrame, T::Int, B::Int)
    """
    return a summary table of estimates
    Input:
        - GP: global parameters
        - dfB: a dataframe of estimates for Bertrand
        - dfC: a dataframe of estimates for Cournot
    Output:
        - CSV files of estimates
    """
    dfB0 = copy(dfB)
    dfC0 = copy(dfC)
    σ_true = GP.σ
    C = GP.C
    N = GP.N
    M = GP.M
    dfB0[!, :DGP] .= "Bertrand"
    dfB0[!, :C] .= C
    dfB0[!, :N] .= N
    dfB0[!, :M] .= M
    dfB0[!, :σ_true] .= σ_true
    col_names = ["t", "delta_DSB", "delta_DSC", "var_DSB", "var_DSC", "hausman_stat_DSB", "hausman_stat_DSC", "sigmaD_OLS", "sigmaD_IV", "sigmaS_OLS", "sigmaSB_GMM", "sigmaSC_GMM","sigmaDSB_GMM", "sigmaDSC_GMM", "seSB_GMM", "seSC_GMM","seDSB_GMM", "seDSC_GMM", "HHI","mean_RS0", "DGP", "C", "N", "M", "sigma_true"]
    # rename columns
    rename!(dfB0, Symbol.(col_names))
    # reordering columns
    dfB0 = dfB0[:, [:DGP, :t, :C, :N, :M, :sigma_true, :sigmaD_OLS, :sigmaD_IV, :sigmaS_OLS, :sigmaSB_GMM, :sigmaSC_GMM, :sigmaDSB_GMM, :sigmaDSC_GMM, :seSB_GMM, :seSC_GMM, :seDSB_GMM, :seDSC_GMM, :HHI, :mean_RS0, :hausman_stat_DSB, :hausman_stat_DSC, :delta_DSB, :delta_DSC, :var_DSB, :var_DSC]]
    CSV.write("Data/Out/CES_MC/main_est/main_est_CSV/dfB_C$(C)_N$(N)_M$(M)_MuY$(GP.μ_Y)_VY$(GP.V_Y)_Muxi$(GP.μ_ξ)_Muomega$(GP.μ_ω)_rho$(GP.ρ)_Vomega$(GP.V_ω)_Vxi$(GP.V_ξ)_Mupsi$(GP.μ_psi)_Vpsi$(GP.V_psi)_beta$(GP.β)_sigma$(GP.σ)_T$(T)_B$(B).csv", dfB0)

    dfC0[!, :DGP] .= "Cournot"
    dfC0[!, :C] .= C
    dfC0[!, :N] .= N
    dfC0[!, :M] .= M
    dfC0[!, :σ_true] .= σ_true
    # rename columns
    rename!(dfC0, Symbol.(col_names))
    # reordering columns
    dfC0 = dfC0[:, [:DGP, :t, :C, :N, :M, :sigma_true, :sigmaD_OLS, :sigmaD_IV, :sigmaS_OLS, :sigmaSB_GMM, :sigmaSC_GMM, :sigmaDSB_GMM, :sigmaDSC_GMM, :seSB_GMM, :seSC_GMM, :seDSB_GMM, :seDSC_GMM, :HHI, :mean_RS0, :hausman_stat_DSB, :hausman_stat_DSC, :delta_DSB, :delta_DSC, :var_DSB, :var_DSC]]
    CSV.write("Data/Out/CES_MC/main_est/main_est_CSV/dfC_C$(C)_N$(N)_M$(M)_MuY$(GP.μ_Y)_VY$(GP.V_Y)_Muxi$(GP.μ_ξ)_Muomega$(GP.μ_ω)_rho$(GP.ρ)_Vomega$(GP.V_ω)_Vxi$(GP.V_ξ)_Mupsi$(GP.μ_psi)_Vpsi$(GP.V_psi)_beta$(GP.β)_sigma$(GP.σ)_T$(T)_B$(B).csv", dfC0)
end


function output_est_separate_0_to_CSV(GP::global_param, dfB::DataFrame, dfC::DataFrame, T::Int, B::Int, folder::String)
    """
    return a summary table of estimates
    Input:
        - GP: global parameters
        - dfB: a dataframe of estimates for Bertrand
        - dfC: a dataframe of estimates for Cournot
    Output:
        - CSV files of estimates
    """
    dfB0 = copy(dfB)
    dfC0 = copy(dfC)
    σ_true = GP.σ
    C = GP.C
    N = GP.N
    M = GP.M
    dfB0[!, :DGP] .= "Bertrand"
    dfB0[!, :C] .= C
    dfB0[!, :N] .= N
    dfB0[!, :M] .= M
    dfB0[!, :σ_true] .= σ_true
    col_names = ["t", "delta_DSB", "delta_DSC", "var_DSB", "var_DSC", "hausman_stat_DSB", "hausman_stat_DSC", "sigmaD_OLS", "sigmaD_IV", "sigmaS_OLS", "sigmaSB_GMM", "sigmaSC_GMM","sigmaDSB_GMM", "sigmaDSC_GMM", "seSB_GMM", "seSC_GMM","seDSB_GMM", "seDSC_GMM", "HHI","mean_RS0", "DGP", "C", "N", "M", "sigma_true"]
    # rename columns
    rename!(dfB0, Symbol.(col_names))
    # reordering columns
    dfB0 = dfB0[:, [:DGP, :t, :C, :N, :M, :sigma_true, :sigmaD_OLS, :sigmaD_IV, :sigmaS_OLS, :sigmaSB_GMM, :sigmaSC_GMM, :sigmaDSB_GMM, :sigmaDSC_GMM, :seSB_GMM, :seSC_GMM, :seDSB_GMM, :seDSC_GMM, :HHI, :mean_RS0, :hausman_stat_DSB, :hausman_stat_DSC, :delta_DSB, :delta_DSC, :var_DSB, :var_DSC]]
    filename = folder * "dfB_C$(C)_N$(N)_M$(M)_MuY$(GP.μ_Y)_VY$(GP.V_Y)_Muxi$(GP.μ_ξ)_Muomega$(GP.μ_ω)_rho$(GP.ρ)_Vomega$(GP.V_ω)_Vxi$(GP.V_ξ)_Mupsi$(GP.μ_psi)_Vpsi$(GP.V_psi)_beta$(GP.β)_sigma$(GP.σ)_Muxi0$(GP.μ_ξ0)_Vxi0$(GP.V_ξ0)_T$(T)_B$(B).csv"
    CSV.write(filename, dfB0)

    dfC0[!, :DGP] .= "Cournot"
    dfC0[!, :C] .= C
    dfC0[!, :N] .= N
    dfC0[!, :M] .= M
    dfC0[!, :σ_true] .= σ_true
    # rename columns
    rename!(dfC0, Symbol.(col_names))
    # reordering columns
    dfC0 = dfC0[:, [:DGP, :t, :C, :N, :M, :sigma_true, :sigmaD_OLS, :sigmaD_IV, :sigmaS_OLS, :sigmaSB_GMM, :sigmaSC_GMM, :sigmaDSB_GMM, :sigmaDSC_GMM, :seSB_GMM, :seSC_GMM, :seDSB_GMM, :seDSC_GMM, :HHI, :mean_RS0, :hausman_stat_DSB, :hausman_stat_DSC, :delta_DSB, :delta_DSC, :var_DSB, :var_DSC]]
    filename = folder * "dfC_C$(C)_N$(N)_M$(M)_MuY$(GP.μ_Y)_VY$(GP.V_Y)_Muxi$(GP.μ_ξ)_Muomega$(GP.μ_ω)_rho$(GP.ρ)_Vomega$(GP.V_ω)_Vxi$(GP.V_ξ)_Mupsi$(GP.μ_psi)_Vpsi$(GP.V_psi)_beta$(GP.β)_sigma$(GP.σ)_Muxi0$(GP.μ_ξ0)_Vxi0$(GP.V_ξ0)_T$(T)_B$(B).csv"
    CSV.write(filename, dfC0)
end


function summary_table(GP::global_param, dfB::DataFrame, dfC::DataFrame)
    """
    return a summary table of estimates
    Input:
        - GP: global parameters
        - dfB: a dataframe of estimates for Bertrand
        - dfC: a dataframe of estimates for Cournot
    Output:
        - summary tables of estimates
    """
    σ_true = GP.σ
    C = GP.C
    N = GP.N
    # calculate t statistics (based on estimated se)
    dfB_t_stats_DSB = (dfB[!, :σDSB_GMM].- σ_true)./dfB[!, :seDSB_GMM]
    dfB_t_stats_DSC = (dfB[!, :σDSC_GMM].- σ_true)./dfB[!, :seDSC_GMM]
    dfC_t_stats_DSB = (dfC[!, :σDSB_GMM].- σ_true)./dfC[!, :seDSB_GMM]
    dfC_t_stats_DSC = (dfC[!, :σDSC_GMM].- σ_true)./dfC[!, :seDSC_GMM]
    # construct the full table
    df_H = DataFrame(SB = ["Bertrand"; C; N; mean(dfB[!, :σD_OLS]); mean(dfB[!, :σD_IV]); mean(dfB[!, :σDSB_GMM]); mean(dfB[!, :σDSC_GMM]); mean(dfB_t_stats_DSB); mean(dfB_t_stats_DSC); mean(dfB[!, :hausman_stat_DSB]); mean(dfB[!, :hausman_stat_DSC]); mean(dfB[!, :mean_RS0])],
    SC = ["Cournot"; C; N; mean(dfC[!, :σD_OLS]); mean(dfC[!, :σD_IV]); mean(dfC[!, :σDSB_GMM]); mean(dfC[!, :σDSC_GMM]); mean(dfC_t_stats_DSB); mean(dfC_t_stats_DSC); mean(dfC[!, :hausman_stat_DSB]); mean(dfC[!, :hausman_stat_DSC]); mean(dfC[!, :mean_RS0])])
   
   return df_H 
end            


function t_stats_density_plot(dfB::DataFrame, dfC::DataFrame, label_list::Array{String,1})
    dfB_t_stats_DSB = (dfB[!, :σDSB_GMM].- σ_true)./dfB[!, :seDSB_GMM]
    dfB_t_stats_DSC = (dfB[!, :σDSC_GMM].- σ_true)./dfB[!, :seDSC_GMM]
    dfC_t_stats_DSB = (dfC[!, :σDSB_GMM].- σ_true)./dfC[!, :seDSB_GMM]
    dfC_t_stats_DSC = (dfC[!, :σDSC_GMM].- σ_true)./dfC[!, :seDSC_GMM]
    gr()  
    xmin = minimum([minimum(dfC_t_stats_DSC), minimum(dfB_t_stats_DSB), minimum(dfB_t_stats_DSC), minimum(dfC_t_stats_DSB)]) - 3
    xmax = maximum([maximum(dfC_t_stats_DSC), maximum(dfB_t_stats_DSB), maximum(dfB_t_stats_DSC), maximum(dfC_t_stats_DSB)]) + 4

    # plotly() # use to create interactive graph, not useful here, also latex strings doesn't work
    p = plot(legend = :topleft, xlims = (xmin, xmax))  # create an empty plot
    # add a density plot for each column, using StatsPlots
 
    StatsPlots.density!(dfB_t_stats_DSB,fill=(0, 0.3,:red), label=label_list[1],linecolor=:red)
    StatsPlots.density!(dfB_t_stats_DSC, label=label_list[2], linecolor=:red, linewidth=2)
    StatsPlots.density!(dfC_t_stats_DSB, label=label_list[3], linecolor=:blue,linewidth=2)
    StatsPlots.density!(dfC_t_stats_DSC, fill=(0, 0.3,:blue), label=label_list[4],linecolor=:blue)
    return p
end

function hausman_stats_density_plot(dfB::DataFrame, dfC::DataFrame,label_list::Array{String,1})
    gr()  
    xmin = minimum([minimum(dfC[!, :hausman_stat_DSC]), minimum(dfB[!, :hausman_stat_DSB]), minimum(dfB[!, :hausman_stat_DSC]), minimum(dfC[!, :hausman_stat_DSB])]) - 3
    xmax = maximum([maximum(dfC[!, :hausman_stat_DSC]), maximum(dfB[!, :hausman_stat_DSB]), maximum(dfB[!, :hausman_stat_DSC]), maximum(dfC[!, :hausman_stat_DSB])]) + 4
    # plotly() # use to create interactive graph, not useful here, also latex strings doesn't work
    p = plot(legend = :topright, xlims = (xmin, xmax))  # create an empty plot
    # add a density plot for each column, using StatsPlots
 
    StatsPlots.density!(dfB[!,:hausman_stat_DSB],fill=(0, 0.3,:red), label=label_list[1],linecolor=:red)
    StatsPlots.density!(dfB[!,:hausman_stat_DSC], label=label_list[2], linecolor=:red, linewidth=2)
    StatsPlots.density!(dfC[!,:hausman_stat_DSB], label=label_list[3], linecolor=:blue,linewidth=2)
    StatsPlots.density!(dfC[!,:hausman_stat_DSC], fill=(0, 0.3,:blue), label=label_list[4],linecolor=:blue)
    
    return p
end

