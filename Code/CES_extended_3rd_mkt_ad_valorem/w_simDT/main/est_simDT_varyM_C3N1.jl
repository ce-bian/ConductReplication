# output estimation results (estimates + Hausman stats etc) for different M, given simulated data. Already run (2024.7.26)
# seems no additional functions need to be changed to suit the new data structure for DGP. 
using Parameters, Optim, ForwardDiff, LinearAlgebra, Distributions, Random
using PrettyTables, DataFrames, NLsolve, NLSolversBase,ThreadsX, Revise
using CSV, JLD2, LaTeXStrings, Statistics
using Plots, StatsPlots

# mkpath("log")
# log_file = open("log/simDT_varyM_parameters3_modified.txt", "w")

# # Redirect stdout and stderr to the log file
# redirect_stdout(log_file)
# redirect_stderr(log_file)


# includet("../CES_simDT_structs.jl")
includet("../../../CES_structs.jl")
# includet("../../../CES_MC/CES_est_diagnostic_functions.jl")
includet("../../../CES_MC/CES_est_functions.jl")
includet("../CES_simDT_estimation_functions.jl")
includet("../CES_simDT_varyM_functions.jl")


# input main parameters
includet("../../../main_parameters.jl")

C = 3
N = 1
M = 1000

w = 0.7 # dampening parameter
MAXIT = 2000 # maximum number of iterations
TOL = 1e-8 # tolerance level
MAX_RETRIES = 500 # maximum number of retries for each DGP


# M_list = [10, 20, 50, 100, 200, 500, 800, 1000]
M_list = [50, 100, 200, 500, 800, 1000]

for selected_M in M_list
    InputDTfile = "Data/eqbaB_eqbaC_GP_C$(C)_N$(N)_M$(M).jld2"  
    # eqbm_output objects 
    eqbm_output_B, eqbm_output_C, GP, seed = load(InputDTfile, "eqbm_output_B", "eqbm_output_C", "GP", "seed");
    eqbm_subset_B = sample_mkt_f(eqbm_output_B, selected_M);
    eqbm_subset_C = sample_mkt_f(eqbm_output_C, selected_M);
    
    GP_new = global_param(C=C, N=N, M=selected_M, μ_Y=μ_Y, V_Y=V_Y, μ_ξ=μ_ξ, V_ξ=V_ξ, μ_ξ0=μ_ξ0, V_ξ0=V_ξ0, μ_ω=μ_ω, V_ω=V_ω, μ_ω0=μ_ω0, V_ω0=V_ω0, ρ=ρ, off1=off1, off2=off2, off10=off10, off20=off20, μ_psi=μ_psi, V_psi=V_psi, β=β, σ=σ_true, w=w, MAXIT=MAXIT, TOL=TOL, MAX_RETRIES=MAX_RETRIES)
    
    dfB = simDT_bootstrap_MC_safe("B", GP_new, eqbm_subset_B, B, full = true)
    dfC = simDT_bootstrap_MC_safe("C", GP_new, eqbm_subset_C, B, full = true)
    # save estimation results
    OutputEstFile = "Data/dfB_dfC_GP_C$(C)_N$(N)_subset_M$(selected_M).jld2"

    save(OutputEstFile, "dfB", dfB, "dfC", dfC, "GP", GP_new)
end

# close(log_file)


# InputDTfile = "Data/eqbaB_eqbaC_GP_C$(C)_N$(N)_M$(M).jld2"  
# # eqbm_output objects 
# eqbm_output_B, eqbm_output_C, GP, seed = load(InputDTfile, "eqbm_output_B", "eqbm_output_C", "GP", "seed");
# eqbm_output_B.eqbm_t_list[2].eqbm_m_list[3]

# C = 3
# InputDTfile = "Data/eqbaB_eqbaC_GP_C$(C)_N$(N)_M$(M).jld2"  
# # eqbm_output objects 
# eqbm_output_B, eqbm_output_C, GP, seed = load(InputDTfile, "eqbm_output_B", "eqbm_output_C", "GP", "seed");
# eqbm_output_B.eqbm_t_list[2].eqbm_m_list[3]