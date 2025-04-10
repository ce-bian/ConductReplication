seed = 1224
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
σ_true = 4.0  

T = 50
B = 250
