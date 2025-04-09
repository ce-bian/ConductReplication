# both DGP and simulation, based on modified structs and DGP functions
using Parameters, Optim, ForwardDiff, LinearAlgebra, Distributions, Random
using PrettyTables, DataFrames, NLsolve, NLSolversBase, RCall, ThreadsX, Revise, DataFramesMeta
using CSV, JLD2,  LaTeXStrings, Glob
includet("../CES_DGP_separate_draw_structs_modified.jl")
includet("../CES_DGP_separate_draw_modified_functions_2.jl")
includet("../CES_est_functions.jl")
includet("../CES_MC_functions.jl")
includet("../CES_MC_DGP_simDT_output_separate_draw_functions.jl")



# parameters 3 as current C2N1 case
seed = 1224
M = 1000 # number of destination markets. -> total markets = C+M
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

# C_N_list = [(4,1), (2,2)]
# C_N_list = [(3,3)]
# C_N_list = [(4,2)]
# C_N_list = [(4,4)]
C_N_list = [(3,1)]


for combo in C_N_list
    C = combo[1]
    N = combo[2]
    GP = global_param(C=C, N=N, M=M, μ_Y=μ_Y, V_Y=V_Y, μ_ξ=μ_ξ, V_ξ=V_ξ, μ_ξ0=μ_ξ0, V_ξ0=V_ξ0, μ_ω=μ_ω, V_ω=V_ω,μ_ω0=μ_ω0, V_ω0=V_ω0, ρ=ρ, off1=off1, off2=off2, off10=off10, off20=off20, μ_psi=μ_psi, V_psi=V_psi, β=β, σ=σ_true, w=w, MAXIT=MAXIT, TOL=TOL, MAX_RETRIES=MAX_RETRIES)

    dfB, dfC = bootstrap_MC_safe(GP, B, T, seed) 
    filename = "Data/Out/CES_MC/main_est/dfB_dfC_GP_C$(C)_N$(N)_M$(M)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_filter_extremes.jld2"
    save(filename, "dfB",dfB, "dfC",dfC, "GP", GP)
end    

for combo in C_N_list
    C = combo[1]
    N = combo[2]
    GP = global_param(C=C, N=N, M=M, μ_Y=μ_Y, V_Y=V_Y, μ_ξ=μ_ξ, V_ξ=V_ξ, μ_ξ0=μ_ξ0, V_ξ0=V_ξ0, μ_ω=μ_ω, V_ω=V_ω,μ_ω0=μ_ω0, V_ω0=V_ω0, ρ=ρ, off1=off1, off2=off2, off10=off10, off20=off20, μ_psi=μ_psi, V_psi=V_psi, β=β, σ=σ_true, w=w, MAXIT=MAXIT, TOL=TOL, MAX_RETRIES=MAX_RETRIES)

    inputFile = "Data/Out/CES_MC/main_est/dfB_dfC_GP_C$(C)_N$(N)_M$(M)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_filter_extremes.jld2"

    outputFile = "Data/Out/CES_MC/main_est/simDT_output/eqbaB_eqbaC_GP_seed$(seed)_C$(C)_N$(N)_M$(M)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_filter_extremes.jld2"

    output_JLD2_file(seed, inputFile, outputFile)
end  