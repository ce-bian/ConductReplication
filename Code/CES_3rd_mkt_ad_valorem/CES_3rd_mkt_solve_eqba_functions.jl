# no need to simulate any data, but input the data directly
# Given parameters, solve for equilibrium prices and quantities in the CES model

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
RS_icm(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, σ) = P_m[c,i]^(1-σ)*exp_ξ_m[c,i]/RS_m(P_m,exp_ξ_m,σ)
q_icm(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, Y_m::Float64, σ) = Y_m* P_m[c,i]^(-σ)*exp_ξ_m[c,i]/RS_m(P_m,exp_ξ_m,σ)
# markups, modified with subsidy/tax
μB_icm(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, S_m::Array{Float64,2},σ) = (σ + (1-σ)*RS_icm(c,i,P_m,exp_ξ_m, σ))/((σ-1)*(1-RS_icm(c,i,P_m,exp_ξ_m, σ))*(1+S_m[c,i]))
μC_icm(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, S_m::Array{Float64,2}, σ) = σ/((σ-1)*(1-RS_icm(c,i,P_m,exp_ξ_m,σ))*(1+S_m[c,i]))
# marginal costs 
MC_icm(c::Int,i::Int,m::Int,τ::Array{Float64,2},exp_ω_m::Array{Float64,2}) = τ[c,m]*exp_ω_m[c,i]


QS_m(Q_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, σ) = Σ((Q_m.^((σ-1)/σ)).*(exp_ξ_m).^(1/σ))
P_icm(c::Int, i::Int, Q_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, σ) = Q_m[c,i]^(-1/σ)*(exp_ξ_m[c,i])^(1/σ)/QS_m(Q_m,exp_ξ_m,σ)



function solve_P_m_R_nolog(R::String, m::Int, C::Int, N::Int, exp_ξ_m::Array{Float64,2}, exp_ω_m::Array{Float64,2}, S_m::Array{Float64,2}, τ::Array{Float64,2}, σ::Float64, w::Float64, MAXIT::Int64, TOL::Float64)
    """
    Fixed point iteration to solve for P_m
        Input:
        - R: "B" (Bertrand) or "C" (Cournot)
        - m: market index
        - C: number of origins
        - N: number of firms
        - exp_ξ_m: C*N array of demand shifters of market m
        - exp_ω_m: C*N array of cost shifters of market m
        - τ: C*(C+M) matrix of trade costs
        - σ: elasticity of substitution
        - w: dampening parameter
        - MAXIT: maximum number of iterations
        - TOL: tolerance level
    Output:
        - P_m_g: C*N matrix of prices conditional on market m
    """

    # # define μR_icm based on input R
    μR_icm = R == "B" ? μB_icm : μC_icm   
    #pR_icm = R == "B" ? PB_icm_safe : PC_icm_safe

    # initial guess: P_m_g: C*N matrix
    P_m_g = [σ/(σ-1)*MC_icm(c,i,m,τ,exp_ω_m)/(1+S_m[c,i]) for c in 1:C, i in 1:N]
    k = 1
    diff = 1.0
        while diff > TOL && k<MAXIT
            P_m_N = [μR_icm(c,i,P_m_g,exp_ξ_m,S_m,σ)*MC_icm(c,i,m,τ,exp_ω_m) for c in 1:C, i in 1:N] 
            # normalize prices, with first firm's price as 1
            # P_m_N = [pR_icm(c,i,m,P_m_g,exp_ξ_m,exp_ω_m,τ,S_m,σ) for c in 1:C, i in 1:N]
            # P_m_N = P_m_N./P_m_N[1,1]
            diff = maximum(abs.(P_m_N-P_m_g))
            P_m_g = w.*P_m_g .+ (1-w).*P_m_N

            k += 1
        end # end while loop  
   
    @assert k < MAXIT "No convergence after $(MAXIT) iterations"
   #  @assert all(P_m_g .> 0) "Negative prices found"

    return P_m_g
end




function DGP_m(R::String, IDT_m::input_dt_m_no_OG)
    """
    Output the simulated data in array format
    Input:
        - R: "B" (Bertrand) or "C" (Cournot)
        - IDT_m: input data of market m 
    Output:
        - equilibrium object for market m
    """
    C = IDT_m.C
    N = IDT_m.N
    m = IDT_m.m   # market index
    σ = IDT_m.σ
    w = IDT_m.w
    TOL = IDT_m.TOL
    MAXIT = IDT_m.MAXIT
    
    exp_ξ_m = IDT_m.exp_ξ_m
    exp_ω_m = IDT_m.exp_ω_m
    S_m = IDT_m.S_m
    τ = IDT_m.τ
    Y_m = IDT_m.Y_m
    
    MC_m = [MC_icm(c,i,m,τ,exp_ω_m) for c in 1:C, i in 1:N]

    P_m = solve_P_m_R_nolog(R, m, C, N, exp_ξ_m, exp_ω_m, S_m, τ, σ, w, MAXIT, TOL)
    RS_m = [RS_icm(c,i,P_m,exp_ξ_m, σ) for c in 1:C, i in 1:N]  

    herfindahl = sum(RS_m.^2)
    
    RV_m = [RS_m[c,i]*Y_m for c in 1:C, i in 1:N]
    Q_m = [RS_m[c,i]*Y_m/P_m[c,i] for c in 1:C, i in 1:N]
    return eqbm_m_no_OG(IDT_m, R, P_m, Q_m, RS_m, MC_m, Y_m, τ, herfindahl, exp_ξ_m, exp_ω_m, S_m, RV_m)
end
