# No S_m involved, with outside good

Σ(x)= sum(x)

"""
Input:
    - m: market index
    - c: origin index
    - i: firm index
    - P_m: C*N array of prices conditional on market m
    - exp_ξ_m: C*N array of demand shifters conditional on market m
    - τ: C*(C+M) matrix of trade costs
    - exp_ω: C*N*(C+M) array of orgin-firm-market specific cost shifters
    - P0_m: price of outside good in market m
    - exp_ξ0_m: demand shifter of outside good of market m
    - exp_ω0_m: cost shifter of outside good of market m
    - σ: elasticity of substitution
"""

# market shares
RS_m(P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, σ) = Σ((P_m.^(1-σ)).*exp_ξ_m)
RS_icm(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, σ) = P_m[c,i]^(1-σ)*exp_ξ_m[c,i]/(P0_m^(1-σ)*exp_ξ0_m+RS_m(P_m,exp_ξ_m,σ))
RS_0mm(P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, σ) = P0_m^(1-σ)*exp_ξ0_m/(P0_m^(1-σ)*exp_ξ0_m+RS_m(P_m,exp_ξ_m,σ))
# markups
μB_icm(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, σ) = 1/((σ-1)*(1-RS_icm(c,i,P_m,exp_ξ_m, P0_m,exp_ξ0_m, σ))) + 1
μC_icm(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, σ) = σ/((σ-1)*(1-RS_icm(c,i,P_m,exp_ξ_m,P0_m,exp_ξ0_m,σ)))
μB_0mm(P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, σ) = 1/((σ-1)*(1-RS_0mm(P_m,exp_ξ_m,P0_m,exp_ξ0_m,σ))) + 1
μC_0mm(P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, σ) = σ/((σ-1)*(1-RS_0mm(P_m,exp_ξ_m,P0_m,exp_ξ0_m,σ)))
# marginal costs (modified)
MC_icm(c::Int,i::Int,m::Int,τ::Array{Float64,2},exp_ω_m::Array{Float64,2}) = τ[c,m]*exp_ω_m[c,i]
MC_0mm(τ::Array{Float64,2},exp_ω0_m::Float64) = 1*exp_ω0_m # τ[m,m] = 1 in current setting
# markups given market shares
μB_icm_R(c::Int,i::Int,RS::Array{Float64,2},σ) = 1/((σ-1)*(1-RS[c,i])) + 1
μC_icm_R(c::Int,i::Int,RS::Array{Float64,2},σ) = σ/((σ-1)*(1-RS[c,i]))
μB_0mm_R(RS0::Float64,σ) = 1/((σ-1)*(1-RS0)) + 1
μC_0mm_R(RS0::Float64,σ) = σ/((σ-1)*(1-RS0))
