using Parameters, Optim, ForwardDiff, LinearAlgebra, Distributions, Random
using PrettyTables, DataFrames, NLsolve, NLSolversBase, ThreadsX, Revise, DataFramesMeta
using CSV, JLD2,  LaTeXStrings, Glob
includet("../../CES_structs.jl")
includet("../CES_DGP_functions.jl")
includet("../CES_est_functions.jl")
includet("../CES_MC_functions.jl")
includet("../CES_MC_DGP_simDT_output_functions.jl")

# input main parameters
includet("../../main_parameters.jl")
M = 1000 # number of destination markets. -> total markets = C+M

w = 0.7 # dampening parameter
MAXIT = 2000 # maximum number of iterations
TOL = 1e-8 # tolerance level
MAX_RETRIES = 500 # maximum number of retries for each DGP

T = 50
B = 250

C_N_list = [(2,1), (3,1)]

# estimate the parameters for each combination of C and N with M = 1000
for combo in C_N_list
    local C = combo[1]
    local N = combo[2]
    GP = global_param(C=C, N=N, M=M, μ_Y=μ_Y, V_Y=V_Y, μ_ξ=μ_ξ, V_ξ=V_ξ, μ_ξ0=μ_ξ0, V_ξ0=V_ξ0, μ_ω=μ_ω, V_ω=V_ω,μ_ω0=μ_ω0, V_ω0=V_ω0, ρ=ρ, off1=off1, off2=off2, off10=off10, off20=off20, μ_psi=μ_psi, V_psi=V_psi, β=β, σ=σ_true, w=w, MAXIT=MAXIT, TOL=TOL, MAX_RETRIES=MAX_RETRIES)

    dfB, dfC = bootstrap_MC_safe(GP, B, T, seed) 
    filename = "Data/dfB_dfC_GP_C$(C)_N$(N)_M$(M).jld2"
    save(filename, "dfB",dfB, "dfC",dfC, "GP", GP)
end    

# output simulated data
for combo in C_N_list
    local C = combo[1]
    local N = combo[2]
    GP = global_param(C=C, N=N, M=M, μ_Y=μ_Y, V_Y=V_Y, μ_ξ=μ_ξ, V_ξ=V_ξ, μ_ξ0=μ_ξ0, V_ξ0=V_ξ0, μ_ω=μ_ω, V_ω=V_ω,μ_ω0=μ_ω0, V_ω0=V_ω0, ρ=ρ, off1=off1, off2=off2, off10=off10, off20=off20, μ_psi=μ_psi, V_psi=V_psi, β=β, σ=σ_true, w=w, MAXIT=MAXIT, TOL=TOL, MAX_RETRIES=MAX_RETRIES)

    inputFile = "Data/dfB_dfC_GP_C$(C)_N$(N)_M$(M).jld2"

    outputFile = "Data/eqbaB_eqbaC_GP_C$(C)_N$(N)_M$(M).jld2"

    output_JLD2_file(seed, inputFile, outputFile)
end  



