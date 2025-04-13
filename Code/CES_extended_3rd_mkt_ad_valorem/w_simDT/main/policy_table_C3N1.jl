using Parameters, Optim, ForwardDiff, LinearAlgebra, Distributions, Random
using PrettyTables, DataFrames, NLsolve, NLSolversBase, ThreadsX, Revise, DataFramesMeta
using CSV, JLD2, Test, Plots, StatsPlots, LaTeXStrings

# includet("../CES_simDT_structs.jl")
includet("../../../CES_structs.jl")
includet("../../CES_solve_eqba_m_functions.jl")
includet("../../CES_eqba_properties_m_functions.jl")
includet("../../CES_welfare_m_functions.jl")
includet("../../CES_optimal_s_det_m_functions.jl")
includet("../CES_policy_table_functions.jl")
includet("../CES_government_estimation_functions.jl")
includet("../CES_simDT_varyM_functions.jl")




# input main parameters
includet("../../../main_parameters.jl")

C = 3
N = 1
M = 1000

# w = 0.7 # dampening parameter
w = 0.95
# MAXIT = 2000 # maximum number of iterations
MAXIT = 5000
TOL = 1e-8 # tolerance level
MAX_RETRIES = 500 # maximum number of retries for each DGP


M_list = [50, 100, 200, 500, 800, 1000]
# M_list = [50, 100, 200]
# M_list = [50, 100, 200]

c = 1 # first market


## based on Hausman test
critical_value_list = [6.635, 3.841, 2.706] # 90%, 95%, 99% confidence intervals
# critical_value_list = [3.841, 6.635] 
# critical_value_list = [6.635] 
for critical_value in critical_value_list
    for M_select in M_list
        # GP = global_param(C=C, N=N, M=M, μ_Y=μ_Y, V_Y=V_Y, μ_ξ=μ_ξ, V_ξ=V_ξ, μ_ξ0=μ_ξ0, V_ξ0=V_ξ0, μ_ω=μ_ω, V_ω=V_ω, ρ=ρ, off1=off1, off2=off2, off10=off10, off20=off20, μ_psi=μ_psi, V_psi=V_psi, β=β, σ=σ_true, w=w, MAXIT=MAXIT, TOL=TOL, MAX_RETRIES=MAX_RETRIES)

        InputDTFile = "Data/eqbaB_eqbaC_GP_C$(C)_N$(N)_M$(M).jld2"

        InputEstFile = "Data/dfB_dfC_GP_C$(C)_N$(N)_subset_M$(M_select).jld2"

        W_dfB, W_dfC = W_hausman_inferred_true_varyM_full_f(InputDTFile, InputEstFile, critical_value, c, M_select)
        OutputFile = "Data/WdfB_WdfC_C$(C)_N$(N)_subset_M$(M_select)_hausman_critical_value$(critical_value).jld2"
        save(OutputFile, "W_dfB", W_dfB, "W_dfC", W_dfC)

    end
end



## based on closeness of the estimated parameters
for M_select in M_list 

    InputDTFile = "Data/eqbaB_eqbaC_GP_C$(C)_N$(N)_M$(M).jld2"

    InputEstFile = "Data/dfB_dfC_GP_C$(C)_N$(N)_subset_M$(M_select).jld2"

    W_dfB, W_dfC = W_closest_inferred_true_varyM_full_f(InputDTFile, InputEstFile, c, M_select)
    OutputFile = "Data/WdfB_WdfC_C$(C)_N$(N)_subset_M$(M_select)_closest.jld2"

    save(OutputFile, "W_dfB", W_dfB, "W_dfC", W_dfC)
end

