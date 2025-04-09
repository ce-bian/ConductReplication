# no need to simulate any data, but input the data directly
# Given parameters, solve for equilibrium prices and quantities in the CES model

@with_kw mutable struct input_dt_no_OG #data
    C::Int64
    N::Int64
    M::Int64
    σ::Float64
    w::Float64 = 0.7
    MAXIT::Int64 = 2000
    TOL::Float64 = 1e-8
    MAX_RETRIES::Int64 = 300
    Y::Array{Float64,1}
    exp_ξ::Array{Float64,3}
    exp_ω::Array{Float64,3}
    τ::Array{Float64,2}
    S::Array{Float64,3} # subsidy/tax
end

@with_kw mutable struct input_dt_m_no_OG #data
    C::Int64
    N::Int64
    m::Int64
    σ::Float64
    w::Float64 = 0.7
    MAXIT::Int64 = 2000
    TOL::Float64 = 1e-8
    MAX_RETRIES::Int64 = 300
    Y_m::Float64
    exp_ξ_m::Array{Float64,2}
    exp_ω_m::Array{Float64,2}
    τ::Array{Float64,2}
    S_m::Array{Float64,2} # subsidy/tax
end


struct eqbm_no_OG #equilibrium objects + input data
    IDT::input_dt_no_OG
    R::String
    P :: Array{Float64,3}
    Q :: Array{Float64,3}
    RS :: Array{Float64,3}
    S :: Array{Float64,3}
    MC:: Array{Float64,3}
    Y :: Array{Float64,1}
    τ :: Array{Float64,2}
    herfindahl :: Float64
end


# memory heavy
struct eqbm_full_no_OG #equilibrium objects + input data
    IDT::input_dt_no_OG
    R::String
    P :: Array{Float64,3}
    Q :: Array{Float64,3}
    RS :: Array{Float64,3}
    MC:: Array{Float64,3}
    Y :: Array{Float64,1}
    τ :: Array{Float64,2}
    herfindahl :: Float64
    exp_ξ ::Array{Float64,3}
    exp_ω :: Array{Float64,3}
    S :: Array{Float64,3}
    RV :: Array{Float64,3}
end

# equilibrium object for a single market
struct eqbm_m_no_OG #equilibrium objects + input data
    IDT::input_dt_m_no_OG
    R::String
    P_m :: Array{Float64,2}
    Q_m :: Array{Float64,2}
    RS_m :: Array{Float64,2}
    MC_m:: Array{Float64,2}
    Y_m :: Float64
    τ :: Array{Float64,2}
    herfindahl :: Float64
    exp_ξ_m ::Array{Float64,2}
    exp_ω_m :: Array{Float64,2}
    S_m :: Array{Float64,2}
    RV_m :: Array{Float64,2}
end




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
# μB_icm_safe(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, σ) = (1/((σ-1)*(1-RS_icm(c,i,P_m,exp_ξ_m, σ))) + 1)/1e16
# μC_icm_safe(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, σ) = σ/(1e16*(σ-1)*(1-RS_icm(c,i,P_m,exp_ξ_m,σ)))
# marginal costs 
MC_icm(c::Int,i::Int,m::Int,τ::Array{Float64,2},exp_ω_m::Array{Float64,2}) = τ[c,m]*exp_ω_m[c,i]
# prices (not used at this point)
PB_icm(c::Int, i::Int,m::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, exp_ω_m::Array{Float64,2},τ::Array{Float64,2},S_m::Array{Float64,2}, σ) = μB_icm(c,i,P_m,exp_ξ_m,S_m,σ)*MC_icm(c,i,m,τ,exp_ω_m)
PC_icm(c::Int, i::Int, m::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, exp_ω_m::Array{Float64,2}, τ::Array{Float64,2},S_m::Array{Float64,2}, σ) = μC_icm(c,i,P_m,exp_ξ_m,S_m,σ)*MC_icm(c,i,m,τ,exp_ω_m)
# PB_icm_safe(c::Int, i::Int,m::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, exp_ω_m::Array{Float64,2},τ::Array{Float64,2}, S_m::Array{Float64,2},σ) = 1e16*μB_icm_safe(c,i,P_m,exp_ξ_m,σ)*MC_icm(c,i,m,τ,exp_ω_m,S_m)
# PC_icm_safe(c::Int, i::Int, m::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, exp_ω_m::Array{Float64,2}, τ::Array{Float64,2}, S_m::Array{Float64,2}, σ) = 1e16*μC_icm_safe(c,i,P_m,exp_ξ_m,σ)*MC_icm(c,i,m,τ,exp_ω_m, S_m)

QS_m(Q_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, σ) = Σ((Q_m.^((σ-1)/σ)).*(exp_ξ_m).^(1/σ))
P_icm(c::Int, i::Int, Q_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, σ) = Q_m[c,i]^(-1/σ)*(exp_ξ_m[c,i])^(1/σ)/QS_m(Q_m,exp_ξ_m,σ)


function solve_P_m_R(R::String, m::Int, C::Int, N::Int, exp_ξ_m::Array{Float64,2}, exp_ω_m::Array{Float64,2}, S_m::Array{Float64,2}, τ::Array{Float64,2}, σ::Float64, w::Float64, MAXIT::Int64, TOL::Float64, log_file::String)
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
    open(log_file, "a") do io
        while diff > TOL && k<MAXIT
            P_m_N = [μR_icm(c,i,P_m_g,exp_ξ_m,σ)*MC_icm(c,i,m,τ,exp_ω_m) for c in 1:C, i in 1:N] 
            # normalize prices, with first firm's price as 1
            # P_m_N = [pR_icm(c,i,m,P_m_g,exp_ξ_m,exp_ω_m,τ,S_m,σ) for c in 1:C, i in 1:N]
            # P_m_N = P_m_N./P_m_N[1,1]
            diff = maximum(abs.(P_m_N-P_m_g))
            for c in 1:C
                for i in 1:N
                    # println(io, "Iteration:", k, "; Price [$(c), $(i), $(m)]: ", P_m_N[c, i], 
                    # "; mu: ", μR_icm(c,i,P_m_g,exp_ξ_m,σ), "; MC: ", MC_icm(c,i,m,τ,exp_ω_m,S_m),
                    # "; RS: ", RS_icm(c,i,P_m_g,exp_ξ_m,σ), "; q: ", q_icm(c,i,P_m_g,exp_ξ_m,σ),
                    # "; Norm: ",diff)
                    println(io, "Iteration:", k, "; Price [$(c), $(i), $(m)]: ", P_m_N[c, i], 
                    "; mu: ", μR_icm(c,i,P_m_g,exp_ξ_m,S_m,σ), "; MC: ", MC_icm(c,i,m,τ,exp_ω_m),
                    "; RS: ", RS_icm(c,i,P_m_g,exp_ξ_m,σ), 
                    "; Norm: ",diff)
                    #println("Iteration: ",k,"; Norm: ",diff)
                end
            end
            P_m_g = w.*P_m_g .+ (1-w).*P_m_N

            k += 1
        end   
    end
   
    @assert k < MAXIT "No convergence after $(MAXIT) iterations"
   #  @assert all(P_m_g .> 0) "Negative prices found"

    return P_m_g
end




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



function DGP(R::String, IDT::input_dt_no_OG, log_file::String, full::Bool)
    """
    Output the simulated data in array format - modified assertion error
    Input:
        - R: "B" (Bertrand) or "C" (Cournot)
        - t: task index
        - IDT: input data
        - full: whether to output full equilibrium objects
    Output:
        - P: C*N*(C+M) array of prices
        - RS: C*N*(C+M) array of revenue shares
        - RV: N*C*M array of revenues (optional)
        - Y: C+M vector of market sizes
        - exp_ξ: C*N*(C+M) array of orgin-firm-market specific demand shifters (optional)
        - exp_ω: C*N*(C+M) array of orgin-firm-market specific cost shifters (optional)
        - τ: C*(C+M) matrix of trade costs
        - herfindahl: Herfindahl index
    """
    C = IDT.C
    N = IDT.N
    M = IDT.M
    σ = IDT.σ
    w = IDT.w
    TOL = IDT.TOL
    MAXIT = IDT.MAXIT
    
    exp_ξ = IDT.exp_ξ
    exp_ω = IDT.exp_ω
    τ = IDT.τ
    Y = IDT.Y
    
    MC = [MC_icm(c,i,m,τ,exp_ω[:,:,m]) for c in 1:C, i in 1:N, m in 1:C+M] # already take into account the subsidy/tax
    # initialize P, RS, RV
    P = zeros(C, N, C+M)
    RS = zeros(C, N, C+M)

    m = 1

    # clear log file
    open(log_file, "w") do io
        println(io, "Start solving for equilibrium prices and quantities")
    end

    while m <= C+M
       P[:,:,m]= solve_P_m_R(R, m, C, N, exp_ξ[:,:,m], exp_ω[:,:,m], S[:,:,m], τ, σ, w, MAXIT, TOL, log_file)
       RS[:,:,m] = [RS_icm(c,i,P[:,:,m],exp_ξ[:,:,m], σ) for c in 1:C, i in 1:N]  
       m += 1
    end
    # add
    herfindahl = sum(RS.^2)
    
    if full
        RV = [RS[c,i,m]*Y[m] for c in 1:C, i in 1:N, m in 1:C+M]
        Q = [RS[c,i,m]*Y[m]/P[c,i,m] for c in 1:C, i in 1:N, m in 1:C+M]
        return eqbm_full_no_OG(IDT, R, P, Q, RS,MC, Y, τ, herfindahl, exp_ξ, exp_ω, S, RV)
    else
        return eqbm(IDT, R, P, Q, RS, S, MC, Y, τ, herfindahl)
    end
end



function DGP_nolog(R::String, IDT::input_dt_no_OG, full::Bool)
    """
    Output the simulated data in array format - modified assertion error
    Input:
        - R: "B" (Bertrand) or "C" (Cournot)
        - t: task index
        - IDT: input data
        - full: whether to output full equilibrium objects
    Output:
        - P: C*N*(C+M) array of prices
        - RS: C*N*(C+M) array of revenue shares
        - RV: N*C*M array of revenues (optional)
        - Y: C+M vector of market sizes
        - exp_ξ: C*N*(C+M) array of orgin-firm-market specific demand shifters (optional)
        - exp_ω: C*N*(C+M) array of orgin-firm-market specific cost shifters (optional)
        - τ: C*(C+M) matrix of trade costs
        - herfindahl: Herfindahl index
    """
    C = IDT.C
    N = IDT.N
    M = IDT.M
    σ = IDT.σ
    w = IDT.w
    TOL = IDT.TOL
    MAXIT = IDT.MAXIT
    
    exp_ξ = IDT.exp_ξ
    exp_ω = IDT.exp_ω
    τ = IDT.τ
    Y = IDT.Y
    
    MC = [MC_icm(c,i,m,τ,exp_ω[:,:,m]) for c in 1:C, i in 1:N, m in 1:C+M]
    # initialize P, RS, RV
    P = zeros(C, N, C+M)
    RS = zeros(C, N, C+M)

    m = 1

    while m <= C+M
       P[:,:,m]= solve_P_m_R_nolog(R, m, C, N, exp_ξ[:,:,m], exp_ω[:,:,m], S[:,:,m], τ, σ, w, MAXIT, TOL)
       RS[:,:,m] = [RS_icm(c,i,P[:,:,m],exp_ξ[:,:,m], σ) for c in 1:C, i in 1:N]  
       m += 1
    end
    # add
    herfindahl = sum(RS.^2)
    
    if full
        RV = [RS[c,i,m]*Y[m] for c in 1:C, i in 1:N, m in 1:C+M]
        Q = [RS[c,i,m]*Y[m]/P[c,i,m] for c in 1:C, i in 1:N, m in 1:C+M]
        return eqbm_full_no_OG(IDT, R, P, Q, RS,MC, Y, τ, herfindahl, exp_ξ, exp_ω, S, RV)
    else
        return eqbm_no_OG(IDT, R, P, Q, RS, S, MC, Y, τ, herfindahl)
    end
end



function DGP_m(R::String, IDT_m::input_dt_m_no_OG)
    """
    Output the simulated data in array format - modified assertion error
    Input:
        - R: "B" (Bertrand) or "C" (Cournot)
        - t: task index
        - IDT_m: input data of market m 
    Output:
        - P_m: C*N mat of prices
        - RS_m: C*N mat of revenue shares
        - RV_m: C*N mat of revenues
        - Y_m: market size
        - exp_ξ: C*N mat of orgin-firm specific demand shifters
        - exp_ω: C*N array of orgin-firm specific cost shifters
        - τ: C*(C+M) matrix of trade costs
        - herfindahl: Herfindahl index
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



# transform equilibrium objects to dataframe
function DT_to_df(eq::Union{eqbm_no_OG, eqbm_full_no_OG}, full::Bool)
    """
    Generate a dataframe for regression
    Input:
       - eq: equilibrium object, can be either eqbm or eqbm_full
       - full: whether to input & output full equilibrium objects
    Output:
        - DT: a dataframe with C*N*(C+M) observations. 
              Each observation is a row with the following columns:
            - c: origin index
            - i: firm index
            - m: market index
            - R_normalized: normalized revenue, LHS variable
            - RS: revenue share
            - RV: revenue (optional)
            - P: relative price
            - τ: relative trade cost 
            - μ_R_true: markup (optional)
            - ξ: demand shifter (optional)
            - ω: cost shifter (optional)
            - v_true: residual (optional)
    """
    C = eq.IDT.C
    N = eq.IDT.N
    M = eq.IDT.M
    R = eq.R
    P = eq.P
    Q = eq.Q
    RS = eq.RS
    MC = eq.MC
    τ = eq.τ
    σ = eq.IDT.σ
    R = eq.R
    Y = eq.Y

    # reshape RS to a C*N*(C+M) vector
    RS_vec = reduce(vcat,reduce(vcat,permutedims(RS, [2,1,3]))) # absolute revenue shares!
    # reshape MC to a C*N*(C+M) vector
    MC_vec = reduce(vcat,reduce(vcat,permutedims(MC, [2,1,3]))) # absolute marginal costs!
    # reshape P to a C*N*(C+M) vector
    P_vec = reduce(vcat,reduce(vcat,permutedims(P, [2,1,3])))
    # reshape Q to a C*N*(C+M) vector
    Q_vec = reduce(vcat,reduce(vcat,permutedims(Q, [2,1,3])))
    # reshape τ to a C*(C+M) vector
    τ_vec = reduce(vcat,τ)
    τ_long = repeat(τ_vec, inner = N) # C*N*(C+M) vector
    c = repeat(repeat(1:C, inner=N), outer=C + M)
    i = repeat(1:N, outer=C * (C + M))
    m = repeat(1:C+M, inner=N * C)
    Y = repeat(Y, inner = N*C)
    S = eq.S
    S_true = reduce(vcat, reduce(vcat, permutedims(S, [2, 1, 3])))
    if full
        RV = eq.RV
        exp_ξ = eq.exp_ξ
        exp_ω = eq.exp_ω
        exp_ξ_true = reduce(vcat, reduce(vcat, permutedims(exp_ξ, [2, 1, 3])))
        ξ = log.(exp_ξ_true)
        exp_ω_true = reduce(vcat, reduce(vcat, permutedims(exp_ω, [2, 1, 3])))
        ω = log.(exp_ω_true)
        RV = reduce(vcat, reduce(vcat, permutedims(RV, [2, 1, 3])))
        μR_icm = R == "B" ? μB_icm : μC_icm     
        μ_array = [μR_icm(c,i,P[:,:,m],exp_ξ[:,:,m],σ) for c in 1:C, i in 1:N, m in 1:C+M]    
        μ_R_true = reduce(vcat, reduce(vcat, permutedims(μ_array, [2, 1, 3])))
        # v_true = ξ .- ξ0_long .+ (1 - σ) * (ω .- ω0_long)
        DT = DataFrame(c=c, i=i, Y_m=Y, m=m, RS=RS_vec,P=P_vec, Q=Q_vec, τ=τ_long, MC=MC_vec, μ_R_true=μ_R_true, ξ=ξ, ω=ω, S=S_true) 
    else
        DT = DataFrame(c=c, i=i, m=m, Y_m=Y, RS=RS_vec, P=P_vec, Q=Q_vec, τ=τ_long, MC= MC_vec, S=S_true)
    end
    return DT
end


# not used?
# function sim_dt_f2(R::String, t::Int64, IDT::input_dt_no_OG, log_file::String, full::Bool)
#     """
#     Return a dataframe of simulated data
#     """
#     eq = DGP(R, IDT, log_file, full)
#     DT = DT_to_df(eq, full)
#     C = IDT.C
#     # add column t
#     DT[!, :t] .= t
#     # only return DT with m > C
#     return DT[DT.m .> C, :]
# end