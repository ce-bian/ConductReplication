using PrettyTables, DataFrames, NLsolve, NLSolversBase, RCall, ThreadsX, Revise, Parameters, Random
using CSV, JLD2, LaTeXStrings, Statistics
using Plots, StatsPlots
includet("CES_DGP_separate_draw_structs.jl")
includet("CES_DGP_separate_draw_functions.jl")
includet("CES_output_functions.jl")


# parameters 2
seed = 1234
M = 1000 # number of destination markets. -> total markets = C+M
# Y follows lognormal distribution
μ_Y = 1.0
V_Y = 1.0
μ_ξ = 0.0
μ_ξ0 = -8.0  # further decrease the appeal of the outside good
μ_ω = 0.0
ρ = 0.5 # correlation between ξ and ω
V_ω = 0.5
V_ξ = 30.0
V_ξ0 = 10.0
off1 = ρ*sqrt(V_ξ)*sqrt(V_ω)
off2 = ρ*sqrt(V_ξ)*sqrt(V_ω)
off10 = ρ*sqrt(V_ξ0)*sqrt(V_ω)
off20 = ρ*sqrt(V_ξ0)*sqrt(V_ω)
μ_psi = 1.0
V_psi = 2.0
β = 1.0
σ_true = 5.0 
w = 0.7 # dampening parameter
MAXIT = 2000 # maximum number of iterations
TOL = 1e-8 # tolerance level
MAX_RETRIES = 300 # maximum number of retries for each DGP

T = 50
B = 250


C_N_list = [(2,1)]

df_H_list = []
for combo in C_N_list
    C = combo[1]
    N = combo[2]
    GP = global_param(C=C, N=N, M=M, μ_Y=μ_Y, V_Y=V_Y, μ_ξ=μ_ξ, V_ξ=V_ξ, μ_ξ0=μ_ξ0, V_ξ0=V_ξ0, μ_ω=μ_ω, V_ω=V_ω,ρ=ρ, off1=off1, off2=off2, off10=off10, off20=off20, μ_psi=μ_psi, V_psi=V_psi, β=β, σ=σ_true, w=w, MAXIT=MAXIT, TOL=TOL, MAX_RETRIES=MAX_RETRIES)

    filename = "Data/Out/CES_MC/main_est/dfB_dfC_GP_C$(C)_N$(N)_M$(M)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_T$(T)_B$(B)_supplement.jld2"
    dfB, dfC, GP = load(filename, "dfB", "dfC", "GP")
    # # output CSV
    # output_est_separate_0_to_CSV(GP, dfB, dfC, T, B)
    # calculate hausman stat
    df_H  = summary_table(GP, dfB, dfC)
    R"""
    library(HeadR)
    df = $df_H
    df <- t(df)
    df <- as.data.frame(df)
    setDT(df)
    df[, (4:ncol(df)) := lapply(.SD, as.numeric), .SDcols = 4:ncol(df)]
    # round columns 4 to 11 with 3 decimal places
    df[, (4:ncol(df)) := lapply(.SD, function(x) round(x, 3)), .SDcols = 4:ncol(df)]
    # set the first column to be character
    setnames(df, 1, "DGP")
    df[, DGP := as.character(DGP)]
    """
    df_H = @rget df
    push!(df_H_list, df_H)

    # # plot density
    # label_list = ["Betrand DGP, Bertrand est", "Bertrand DGP, Cournot est", "Cournot DGP, Bertrand est", "Cournot DGP, Cournot est"]
    # p = hausman_stats_density_plot(dfB, dfC, label_list)
    # filename = "Figures/CES_MC/main_est_varyCN/hausman_density_plot_supply_demand_C$(C)_N$(N)_M$(M)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_T$(T)_B$(B).pdf"
    # StatsPlots.savefig(p, filename)

    # p = t_stats_density_plot(dfB, dfC, label_list)
    # filename = "Figures/CES_MC/main_est_varyCN/t_density_plot_supply_demand_C$(C)_N$(N)_M$(M)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_T$(T)_B$(B).pdf"
    # StatsPlots.savefig(p, filename)
end
df_H = vcat(df_H_list...)

# problem: μ_ξ0 = -10.0 leads to insignificant hausman stat for Betrand (3.978 only)


