# seems no additional functions need to be changed to suit the new data structure for DGP. Already run (2024.7.26)
using Parameters, Optim, ForwardDiff, LinearAlgebra, Distributions, Random
using PrettyTables, DataFrames, NLsolve, NLSolversBase, RCall, ThreadsX, Revise, DataFramesMeta
using CSV, JLD2, Test, Plots, StatsPlots, LaTeXStrings

includet("../CES_simDT_structs.jl")
includet("../../CES_solve_eqba_m_functions.jl")
includet("../../CES_eqba_properties_m_functions.jl")
includet("../../CES_welfare_m_functions.jl")
includet("../../CES_optimal_s_det_m_functions.jl")
includet("../CES_policy_table_functions.jl")
includet("../CES_government_estimation_functions.jl")
includet("../CES_simDT_varyM_functions.jl")


############################## parameter 3 ################################
seed = 1224
C = 2
N = 1
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


# w = 0.7 # dampening parameter
w = 0.95
# MAXIT = 2000 # maximum number of iterations
MAXIT = 5000
TOL = 1e-8 # tolerance level
MAX_RETRIES = 500 # maximum number of retries for each DGP

T = 50
B = 250

# M_list = [10, 20, 50, 100, 200, 500, 800, 1000]
M_list = [50, 100, 200]
# M_list = [50, 100, 200]

c = 1 # first market


## based on Hausman test
# critical_value_list = [2.706, 3.841, 6.635] # 90%, 95%, 99% confidence intervals
critical_value_list = [3.841, 6.635] 
# critical_value_list = [6.635] 
for critical_value in critical_value_list
    for M_select in M_list
        # GP = global_param(C=C, N=N, M=M, μ_Y=μ_Y, V_Y=V_Y, μ_ξ=μ_ξ, V_ξ=V_ξ, μ_ξ0=μ_ξ0, V_ξ0=V_ξ0, μ_ω=μ_ω, V_ω=V_ω, ρ=ρ, off1=off1, off2=off2, off10=off10, off20=off20, μ_psi=μ_psi, V_psi=V_psi, β=β, σ=σ_true, w=w, MAXIT=MAXIT, TOL=TOL, MAX_RETRIES=MAX_RETRIES)

        InputDTFile = "Data/Out/CES_MC/main_est/simDT_output/eqbaB_eqbaC_GP_seed$(seed)_C$(C)_N$(N)_M$(M)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_filter_extremes.jld2"

        InputEstFile = "Data/Out/CES_MC/main_est/C2N1_varyM/C2N1_varyM_est/new/dfB_dfC_GP_C$(C)_N$(N)_M$(M_select)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_parameter3_retry.jld2"

        W_dfB, W_dfC = W_hausman_inferred_true_varyM_full_f(InputDTFile, InputEstFile, critical_value, c, M)
        OutputFile = "Data/Out/CES_MC/main_est/policy_table_ad_valorem/C2N1_varyM/WdfB_WdfC_C$(C)_N$(N)_M$(M_select)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_hausman_critical_value$(critical_value)_parameter3_retry_new.jld2"
        save(OutputFile, "W_dfB", W_dfB, "W_dfC", W_dfC)

    end
end



## based on closeness of the estimated parameters
for M_select in M_list 

    InputDTFile = "Data/Out/CES_MC/main_est/simDT_output/eqbaB_eqbaC_GP_seed$(seed)_C$(C)_N$(N)_M$(M)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_filter_extremes.jld2"

    InputEstFile = "Data/Out/CES_MC/main_est/C2N1_varyM/C2N1_varyM_est/new/dfB_dfC_GP_C$(C)_N$(N)_M$(M_select)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_parameter3_retry.jld2"

    W_dfB, W_dfC = W_closest_inferred_true_varyM_full_f(InputDTFile, InputEstFile, c, M)
    OutputFile = "Data/Out/CES_MC/main_est/policy_table_ad_valorem/C2N1_varyM/WdfB_WdfC_C$(C)_N$(N)_M$(M_select)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_closest_parameter3_retry_new.jld2"

    save(OutputFile, "W_dfB", W_dfB, "W_dfC", W_dfC)
end

