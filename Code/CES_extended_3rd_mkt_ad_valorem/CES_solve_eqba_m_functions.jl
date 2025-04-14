# no need to simulate any data; input the data directly
# Given a market and parameters, solve for equilibrium prices and quantities in the CES model, with outside good 
# each country (C)'s subsidies/taxes are market specific
# solve for individual market equilibrium

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
"""
# market shares
RS_m(P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, σ) = Σ((P_m.^(1-σ)).*exp_ξ_m)
RS_icm(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, σ) = P_m[c,i]^(1-σ)*exp_ξ_m[c,i]/(P0_m^(1-σ)*exp_ξ0_m+RS_m(P_m,exp_ξ_m,σ))
RS_0mm(P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, σ) = P0_m^(1-σ)*exp_ξ0_m/(P0_m^(1-σ)*exp_ξ0_m+RS_m(P_m,exp_ξ_m,σ))
q_icm(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, Y_m::Float64, σ) = Y_m* P_m[c,i]^(-σ)*exp_ξ_m[c,i]/(P0_m^(1-σ)*exp_ξ0_m+RS_m(P_m,exp_ξ_m,σ))
# markups, updated
μB_icm(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, S_m::Vector{Float64}, σ) = (σ + (1-σ)*RS_icm(c,i,P_m,exp_ξ_m, P0_m,exp_ξ0_m, σ))/((σ-1)*(1-RS_icm(c,i,P_m,exp_ξ_m, P0_m,exp_ξ0_m, σ))*(1+S_m[c]))
μC_icm(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, S_m::Vector{Float64}, σ) = σ/((σ-1)*(1-RS_icm(c,i,P_m,exp_ξ_m,P0_m,exp_ξ0_m,σ))*(1+S_m[c]))
μB_0mm(P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, σ) = 1/((σ-1)*(1-RS_0mm(P_m,exp_ξ_m,P0_m,exp_ξ0_m,σ))) + 1
μC_0mm(P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, σ) = σ/((σ-1)*(1-RS_0mm(P_m,exp_ξ_m,P0_m,exp_ξ0_m,σ)))
# marginal costs
MC_icm(c::Int,i::Int,m::Int,τ::Array{Float64,2},exp_ω_m::Array{Float64,2}) = τ[c,m]*exp_ω_m[c,i]
MC_0mm(τ::Array{Float64,2},exp_ω0_m::Float64) = 1*exp_ω0_m # τ[m,m] = 1 in current setting
# prices, take into account ad valorem subsidy 
PB_icm(c::Int, i::Int,m::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, exp_ω_m::Array{Float64,2},τ::Array{Float64,2},S_m::Vector{Float64}, σ) = μB_icm(c,i,P_m,exp_ξ_m,P0_m,exp_ξ0_m,S_m,σ)*MC_icm(c,i,m,τ,exp_ω_m)
PC_icm(c::Int, i::Int, m::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, exp_ω_m::Array{Float64,2}, τ::Array{Float64,2},S_m::Vector{Float64}, σ) = μC_icm(c,i,P_m,exp_ξ_m,P0_m,exp_ξ0_m,S_m,σ)*MC_icm(c,i,m,τ,exp_ω_m)


QS_m(Q_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, σ) = Σ((Q_m.^((σ-1)/σ)).*(exp_ξ_m).^(1/σ))
P_icm(c::Int, i::Int, Q_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, Q0_m::Float64, exp_ξ0_m::Float64, Y_m::Float64, σ) = Y_m* Q_m[c,i]^(-1/σ)*(exp_ξ_m[c,i])^(1/σ)/(Q0_m^((σ-1)/σ)*exp_ξ0_m^(1/σ) + QS_m(Q_m,exp_ξ_m,σ))
P0_icm(Q_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, Q0_m::Float64, exp_ξ0_m::Float64, σ) = Q0_m^(-1/σ)*exp_ξ0_m^(1/σ)/(Q0_m^((σ-1)/σ)*exp_ξ0_m^(1/σ) + QS_m(Q_m,exp_ξ_m,σ))



function solve_P_m_R_nolog(R::String, m::Int, C::Int, N::Int, exp_ξ_m::Array{Float64,2}, exp_ω_m::Array{Float64,2}, exp_ξ0_m::Float64, exp_ω0_m::Float64, S_m::Vector{Float64}, τ::Array{Float64,2}, σ::Float64, w::Float64, MAXIT::Int64, TOL::Float64)
    """
    Fixed point iteration to solve for P_m
        Input:
        - R: "B" (Bertrand) or "C" (Cournot)
        - m: market index
        - C: number of origins
        - N: number of firms
        - exp_ξ_m: C*N array of demand shifters of market m
        - exp_ξ0_m: C array of demand shifters of outside good
        - exp_ω_m: C*N array of cost shifters of market m
        - exp_ω0_m: C array of cost shifters of outside good
        - τ: C*(C+M) matrix of trade costs
        - σ: elasticity of substitution
        - w: dampening parameter
        - MAXIT: maximum number of iterations
        - TOL: tolerance level
    Output:
        - P_m_g: C*N matrix of prices conditional on market m
        - P0_m_g: float of prices of outside good
    """

    # # define μR_icm based on input R
    μR_icm = R == "B" ? μB_icm : μC_icm   
    μR_0mm = R == "B" ? μB_0mm : μC_0mm   
    #pR_icm = R == "B" ? PB_icm_safe : PC_icm_safe

    # initial guess: P_m_g: C*N matrix
    P_m_g = [σ/(σ-1)*MC_icm(c,i,m,τ,exp_ω_m)/(1+S_m[c]) for c in 1:C, i in 1:N]
    P0_m_g = σ/(σ-1)*MC_0mm(τ,exp_ω0_m)
    k = 1
    diff = 1.0
        while diff > TOL && k<MAXIT
            P_m_N = [μR_icm(c,i,P_m_g,exp_ξ_m,P0_m_g,exp_ξ0_m,S_m,σ)*MC_icm(c,i,m,τ,exp_ω_m) for c in 1:C, i in 1:N] 
            P0_m_N = μR_0mm(P_m_g,exp_ξ_m,P0_m_g,exp_ξ0_m,σ)*MC_0mm(τ,exp_ω0_m)
            diff = maximum(max.(abs.(P_m_N-P_m_g),abs(P0_m_N-P0_m_g)))
            P_m_g = w.*P_m_g .+ (1-w).*P_m_N
            P0_m_g = w*P0_m_g + (1-w)*P0_m_N
            k += 1
        end # end while loop  
   
    @assert k < MAXIT "No convergence after $(MAXIT) iterations"
    @assert all(P_m_g .> 0) "Negative prices found"
    @assert all(P0_m_g .> 0) "Negative prices found"

    return P_m_g, P0_m_g
end
 

function DGP_m(R::String, IDT_m::input_dt_m)
    """
    Output the simulated data
    Input:
        - R: "B" (Bertrand) or "C" (Cournot)
        - IDT_m: input data of market m 
    Output:
        - eqbm_m: equilibrium object for market m
    """
    C = IDT_m.C
    N = IDT_m.N
    m = IDT_m.m   # market index
    σ = IDT_m.σ
    # w = IDT_m.w
    w = 0.95
    TOL = IDT_m.TOL
    # MAXIT = IDT_m.MAXIT
    MAXIT = 5000
    
    exp_ξ_m = IDT_m.exp_ξ_m
    exp_ω_m = IDT_m.exp_ω_m
    exp_ξ0_m = IDT_m.exp_ξ0_m
    exp_ω0_m = IDT_m.exp_ω0_m
    S_m = IDT_m.S_m
    τ = IDT_m.τ
    Y_m = IDT_m.Y_m
    
    MC_m = [MC_icm(c,i,m,τ,exp_ω_m) for c in 1:C, i in 1:N]
    MC0_m = MC_0mm(τ,exp_ω0_m)

    P_m, P0_m = solve_P_m_R_nolog(R, m, C, N, exp_ξ_m, exp_ω_m, exp_ξ0_m, exp_ω0_m, S_m, τ, σ, w, MAXIT, TOL)
    RS_mat = [RS_icm(c,i,P_m,exp_ξ_m,P0_m, exp_ξ0_m, σ) for c in 1:C, i in 1:N]  
    RS0_m = RS_0mm(P_m,exp_ξ_m,P0_m,exp_ξ0_m,σ)
    herfindahl = sum(RS_mat.^2) + RS0_m^2
    
    RV_m = [RS_mat[c,i]*Y_m for c in 1:C, i in 1:N]
    RV0_m = RS0_m*Y_m
    Q_m = [RS_mat[c,i]*Y_m/P_m[c,i] for c in 1:C, i in 1:N]
    Q0_m = RS0_m*Y_m/P0_m

    return eqbm_m(IDT_m, R, P_m, P0_m, Q_m, Q0_m, RS_mat, RS0_m, MC_m, MC0_m, herfindahl, RV_m, RV0_m)
end

