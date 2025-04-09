# incomplete
# output simDT for policy analysis
# correspond to table 1 

# baseline estimates 
seed = 1234
M = 1000 # number of destination markets. -> total markets = C+M
# Y follows lognormal distribution
μ_Y = 1.0
V_Y = 1.0
μ_ξ = 0.0
μ_ξ0 = -2.0
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

C_N_list = [(2,1), (2,2), (2,3), (3,2), (3,3)]