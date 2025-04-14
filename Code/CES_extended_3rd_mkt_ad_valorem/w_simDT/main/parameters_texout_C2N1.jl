using Parameters, DataFrames, Revise, DataFramesMeta, Statistics
using JLD2, LaTeXStrings, PrettyTables

# input main parameters
includet("../../../main_parameters.jl")

C = 3
N = 1
M = 1000

w = 0.7 # dampening parameter
MAXIT = 2000 # maximum number of iterations
TOL = 1e-8 # tolerance level
MAX_RETRIES = 500 # maximum number of retries for each DGP



# output a prettytable tex file
# first column: variable name, second column: mu, third column: sigma

df = DataFrame(variable = [L"Y_m", L"\xi_{icm}", L"\xi_{0mm}", L"\omega_{icm}", L"\omega_{0mm}", L"\psi_{cm}", L"\sigma", L"\rho(\xi_{icm}, \omega_{icm})", L"\rho(\xi_{0mm}, \omega_{0mm})"],
    value = ["", "", "", "", "", "", σ_true,  ρ, ρ],
    mu = [μ_Y, μ_ξ, μ_ξ0, μ_ω, μ_ω0, μ_psi, "", "", ""],
    sigma = [V_Y, V_ξ, V_ξ0, V_ω, V_ω0, V_psi, "", "", ""],
    note = ["Market size", "Demand shifter of importing goods", "Demand shifter of local goods", "Cost shifter of exporters", "Cost shifter of local firms", "Trade cost shifter", "Elasticity", "Demand-cost correlation of importing goods", "Demand-cost correlation of local goods"])
# transform the dataframe to a matrix
# df = Matrix(df)
open("Tables/Table1.tex", "w") do io
    pretty_table(io, df, header=["variable", "value", "mean", "variance", "definition"],  backend = Val(:latex))
end