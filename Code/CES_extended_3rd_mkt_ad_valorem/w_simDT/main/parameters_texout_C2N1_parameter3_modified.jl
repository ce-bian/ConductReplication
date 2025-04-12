using Parameters, DataFrames, RCall, Revise, DataFramesMeta, Statistics
using JLD2, LaTeXStrings, PrettyTables

includet("../texout_f.jl")
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
w = 0.7 # dampening parameter
MAXIT = 2000 # maximum number of iterations
TOL = 1e-8 # tolerance level
MAX_RETRIES = 500 # maximum number of retries for each DGP

T = 50
B = 250


# output a prettytable tex file
# first column: variable name, second column: mu, third column: sigma

df = DataFrame(variable = [L"Y_m", L"\xi_{icm}", L"\xi_{0mm}", L"\omega_{icm}", L"\omega_{0mm}", L"\psi_{cm}", L"\sigma", L"\rho(\xi_{icm}, \omega_{icm})", L"\rho(\xi_{0mm}, \omega_{0mm})"],
    value = ["", "", "", "", "", "", σ_true,  ρ, ρ],
    mu = [μ_Y, μ_ξ, μ_ξ0, μ_ω, μ_ω, μ_psi, "", "", ""],
    sigma = [V_Y, V_ξ, V_ξ0, V_ω, V_ω, V_psi, "", "", ""],
    note = ["Market size", "Demand shifter of importing goods", "Demand shifter of local goods", "Cost shifter of exporters", "Cost shifter of local firms", "Trade cost shifter", "Elasticity", "Demand-cost correlation of importing goods", "Demand-cost correlation of local goods"])
# transform the dataframe to a matrix
# df = Matrix(df)
open("Tables/CES_MC/main_est_varyCN/policy_table_ad_valorem/C2N1_parameters_table_parameter3_modified.tex", "w") do io
    pretty_table(io, df, header=["variable", "value", "mean", "variance", "definition"],  backend = Val(:latex))
end