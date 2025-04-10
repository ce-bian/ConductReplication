# same distributional parameters as in parameter 3 for C2N1 case
using Parameters, Optim, ForwardDiff, LinearAlgebra, Distributions, Random
using PrettyTables, DataFrames, NLsolve, NLSolversBase, RCall, ThreadsX, Revise
using CSV, JLD2, LaTeXStrings, Statistics
using Plots, StatsPlots

includet("../CES_simDT_structs.jl")
# includet("../../../CES_MC/CES_est_diagnostic_functions.jl")
includet("../../../CES_MC/CES_est_functions.jl")
includet("../CES_simDT_estimation_diagnostic_functions.jl")
includet("../CES_simDT_varyM_functions.jl")


############################## parameter 3 ################################
seed = 1224
C = 3
N = 1  # change here
M = 1000
# Y follows lognormal distribution
μ_Y = 1.0
V_Y = 1.0
μ_ξ = 0.0
μ_ξ0 = -3.0  # make the mean of ξ0 small
μ_ω = 0.0
μ_ω0 = 1.5   # make the mean of ω0 large  
ρ = 0.5 # correlation between ξ and ω
V_ω = 0.5
V_ω0 = 0.1  # add
V_ξ = 30.0
V_ξ0 = 1.0 # also decrease the variance of ξ0
off1 = ρ*sqrt(V_ξ)*sqrt(V_ω)
off2 = ρ*sqrt(V_ξ)*sqrt(V_ω)
off10 = ρ*sqrt(V_ξ0)*sqrt(V_ω0)
off20 = ρ*sqrt(V_ξ0)*sqrt(V_ω0)
μ_psi = 1.0
V_psi = 2.0
β = 1.0
σ_true = 4.0  # modified 
w = 0.7 # dampening parameter
MAXIT = 2000 # maximum number of iterations
TOL = 1e-8 # tolerance level
MAX_RETRIES = 500 # maximum number of retries for each DGP

T = 50
B = 250



M_list = [50, 100, 200, 500, 800, 1000]
# M_list = [1000, 800, 500]

for selected_M in M_list
    InputDTfile = "Data/Out/CES_MC/main_est/simDT_output/eqbaB_eqbaC_GP_seed$(seed)_C$(C)_N$(N)_M$(M)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_filter_extremes.jld2"    
    # eqbm_output objects 
    eqbm_output_B, eqbm_output_C, GP, seed = load(InputDTfile, "eqbm_output_B", "eqbm_output_C", "GP", "seed");
    eqbm_subset_B = sample_mkt_f(eqbm_output_B, selected_M);
    eqbm_subset_C = sample_mkt_f(eqbm_output_C, selected_M);
    
    # length(eqbm_subset_B.eqbm_t_list)
    # length(eqbm_subset_B.eqbm_t_list[1].eqbm_m_list)
    # update M in GP
    GP_new = global_param(C=C, N=N, M=selected_M, μ_Y=μ_Y, V_Y=V_Y, μ_ξ=μ_ξ, V_ξ=V_ξ, μ_ξ0=μ_ξ0, V_ξ0=V_ξ0, μ_ω=μ_ω, V_ω=V_ω,ρ=ρ, off1=off1, off2=off2, off10=off10, off20=off20, μ_psi=μ_psi, V_psi=V_psi, β=β, σ=σ_true, w=w, MAXIT=MAXIT, TOL=TOL, MAX_RETRIES=MAX_RETRIES)
    
    dfB = simDT_bootstrap_MC_safe("B", GP_new, eqbm_subset_B, B, full = true)
    dfC = simDT_bootstrap_MC_safe("C", GP_new, eqbm_subset_C, B, full = true)
    # save estimation results
    OutputEstFile = "Data/Out/CES_MC/main_est/generalized_varyM/generalized_varyM_est/new/dfB_dfC_GP_C$(C)_N$(N)_M$(selected_M)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_case1_retry.jld2"

    save(OutputEstFile, "dfB", dfB, "dfC", dfC, "GP", GP_new)
end