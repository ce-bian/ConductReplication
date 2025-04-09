# this version still assumes N is common across all origins.
# using Parameters, LinearAlgebra, Distributions, Random, DataFrames, ThreadsX

Σ(x)= sum(x)

function set_θ(GP::global_param, task_rng::AbstractRNG)
    """
    Get the rest of parameters in "θ"
    Input: global parameters (GP) and task-specific random number generator (task_rng)
    Output:
        - Y: (C+M) vector of incomes
        - exp_ξ: C*N*(C+M) array of orgin-firm-market specific demand shifters
        - exp_ω: C*N*(C+M) array of orgin-firm-market specific cost shifters
        - τ: C*(C+M) matrix of orgin-market trade costs
    """
    C = GP.C
    N = GP.N
    M = GP.M
    μ_Y = GP.μ_Y
    V_Y = GP.V_Y
    μ_ξ = GP.μ_ξ
    V_ξ = GP.V_ξ
    μ_ξ0 = GP.μ_ξ0
    V_ξ0 = GP.V_ξ0
    μ_ω = GP.μ_ω
    V_ω = GP.V_ω
    off1 = GP.off1
    off2 = GP.off2
    off10 = GP.off10
    off20 = GP.off20
    μ_psi = GP.μ_psi
    V_psi = GP.V_psi
    β = GP.β
    Y = rand(task_rng,LogNormal(μ_Y, V_Y),C+M)
    # V_mvnormal = MvNormal([1.0;1.0],[V_ξ_ω ρ*V_ξ_ω; ρ*V_ξ_ω V_ξ_ω])
    V_mvnormal = MvNormal([μ_ξ;μ_ω],[V_ξ off1; off2 V_ω])
    temp = rand(task_rng, V_mvnormal,N*C*(C+M))'
    ξ_long = temp[:,1]
    ω_long = temp[:,2]
    ξ = zeros(C,N,(C+M))
    ω = zeros(C,N,(C+M))
    # reshape ξ and ω
    for m in 1:C+M
        for c in 1:C
            for i in 1:N
                ξ[c,i,m] = ξ_long[C*N*(m-1) + N*(c-1) + i]
                ω[c,i,m] = ω_long[C*N*(m-1) + N*(c-1) + i]
            end
        end
    end 
    exp_ξ = exp.(ξ)
    exp_ω = exp.(ω)
    ψ = rand(task_rng, Normal(μ_psi, V_psi),C,C+M)
    τ = ones(C,C+M)
    for i in 1:C
        for j in 1:C+M
            if i != j
                τ[i,j] = (1+exp(ψ[i,j]))^β
            end
        end
    end
    # outside good
    V0_mvnormal = MvNormal([μ_ξ0;μ_ω],[V_ξ0 off10; off20 V_ω])
    temp = rand(task_rng, V0_mvnormal,C+M)'
    ξ_long = temp[:,1]
    ω_long = temp[:,2]
    exp_ξ0 = exp.(ξ_long)
    exp_ω0 = exp.(ω_long)
   return sim_dt(GP,Y,exp_ξ,exp_ω,exp_ξ0,exp_ω0,τ)
end



function set_θ_update_m(m::Int, SDT_old::sim_dt, GP::global_param, task_rng::AbstractRNG)
    """
    Get the rest of parameters in "θ", update m-th market
    Input: 
    Output:
        - Y: (C+M) vector of incomes
        - exp_ξ: C*N*(C+M) array of orgin-firm-market specific demand shifters
        - exp_ω: C*N*(C+M) array of orgin-firm-market specific cost shifters
        - τ: C*(C+M) matrix of orgin-market trade costs
    """
    C = GP.C
    N = GP.N
    # M = GP.M
    μ_Y = GP.μ_Y
    V_Y = GP.V_Y
    μ_ξ = GP.μ_ξ
    V_ξ = GP.V_ξ
    μ_ξ0 = GP.μ_ξ0
    V_ξ0 = GP.V_ξ0
    μ_ω = GP.μ_ω
    V_ω = GP.V_ω
    off1 = GP.off1
    off2 = GP.off2
    off10 = GP.off10
    off20 = GP.off20
   
    Y = SDT_old.Y
    ξ = log.(SDT_old.exp_ξ)
    ω = log.(SDT_old.exp_ω)
    ξ0 = log.(SDT_old.exp_ξ0) # add
    ω0 = log.(SDT_old.exp_ω0)
    τ = SDT_old.τ
    Y_m = rand(task_rng,LogNormal(μ_Y, V_Y))
    # update Y
    Y[m] = Y_m

    V_mvnormal = MvNormal([μ_ξ;μ_ω],[V_ξ off1; off2 V_ω])
    temp_m = rand(task_rng, V_mvnormal, N*C)' 
    ξ_long = temp_m[:,1]
    ω_long = temp_m[:,2]

    # update ξ and ω
    for c in 1:C
        for i in 1:N
            ξ[c,i,m] = ξ_long[N*(c-1) + i]
            ω[c,i,m] = ω_long[N*(c-1) + i]
        end
    end
    exp_ξ = exp.(ξ)
    exp_ω = exp.(ω)

    # update outside good for m-th market
    V0_mvnormal = MvNormal([μ_ξ0;μ_ω],[V_ξ0 off10; off20 V_ω])
    temp0_m = rand(task_rng, V0_mvnormal)
    ξ0[m] = temp0_m[1]
    ω0[m] = temp0_m[2]
    exp_ξ0 = exp.(ξ0)
    exp_ω0 = exp.(ω0)

    return sim_dt(GP,Y,exp_ξ,exp_ω,exp_ξ0,exp_ω0,τ)
end


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
# markups
μB_icm(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, σ) = 1/((σ-1)*(1-RS_icm(c,i,P_m,exp_ξ_m, P0_m,exp_ξ0_m, σ))) + 1
μC_icm(c::Int, i::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, σ) = σ/((σ-1)*(1-RS_icm(c,i,P_m,exp_ξ_m,P0_m,exp_ξ0_m,σ)))
μB_0mm(P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, σ) = 1/((σ-1)*(1-RS_0mm(P_m,exp_ξ_m,P0_m,exp_ξ0_m,σ))) + 1
μC_0mm(P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, σ) = σ/((σ-1)*(1-RS_0mm(P_m,exp_ξ_m,P0_m,exp_ξ0_m,σ)))
# marginal costs (modified)
MC_icm(c::Int,i::Int,m::Int,τ::Array{Float64,2},exp_ω_m::Array{Float64,2}) = τ[c,m]*exp_ω_m[c,i]
MC_0mm(τ::Array{Float64,2},exp_ω0_m::Float64) = 1*exp_ω0_m # τ[m,m] = 1 in current setting
# prices (not used at this point)
PB_icm(c::Int, i::Int,m::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, exp_ω_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, τ::Array{Float64,2},σ) = μB_icm(c,i,P_m,exp_ξ_m,P0_m,exp_ξ0_m,σ)*MC_icm(c,i,m,τ,exp_ω_m)
PC_icm(c::Int, i::Int, m::Int, P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, exp_ω_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, τ::Array{Float64,2},σ) = μC_icm(c,i,P_m,exp_ξ_m,P0_m,exp_ξ0_m,σ)*MC_icm(c,i,m,τ,exp_ω_m)
PB_0mm(P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, exp_ω0_m::Float64, τ::Array{Float64,2},σ) = μB_0mm(P_m,exp_ξ_m,P0_m,exp_ξ0_m,σ)*MC_0mm(τ,exp_ω0_m)
PC_0mm(P_m::Array{Float64,2}, exp_ξ_m::Array{Float64,2}, P0_m::Float64, exp_ξ0_m::Float64, exp_ω0_m::Float64, τ::Array{Float64,2},σ) = μC_0mm(P_m,exp_ξ_m,P0_m,exp_ξ0_m,σ)*MC_0mm(τ,exp_ω0_m)
# markups given market shares
μB_icm_R(c::Int,i::Int,RS::Array{Float64,2},σ) = 1/((σ-1)*(1-RS[c,i])) + 1
μC_icm_R(c::Int,i::Int,RS::Array{Float64,2},σ) = σ/((σ-1)*(1-RS[c,i]))
μB_0mm_R(RS0::Float64,σ) = 1/((σ-1)*(1-RS0)) + 1
μC_0mm_R(RS0::Float64,σ) = σ/((σ-1)*(1-RS0))



function solve_P_m_R(R::String, m::Int, C::Int, N::Int, exp_ξ_m::Array{Float64,2}, exp_ω_m::Array{Float64,2}, exp_ξ0_m::Float64, exp_ω0_m::Float64, τ::Array{Float64,2}, σ::Float64, w::Float64, MAXIT::Int64, TOL::Float64)
    """
    Fixed point iteration to solve for P_m
        Input:
        - R: "B" (Bertrand) or "C" (Cournot)
        - m: market index
        - C: number of origins
        - N: number of firms
        - exp_ξ_m: C*N array of demand shifters of market m
        - exp_ω_m: C*N array of cost shifters of market m
        - exp_ξ0_m: demand shifter of outside good of market m
        - exp_ω0_m: cost shifter of outside good of market m
        - τ: C*(C+M) matrix of trade costs
        - σ: elasticity of substitution
        - w: dampening parameter
        - MAXIT: maximum number of iterations
        - TOL: tolerance level
    Output:
        - P_m_g: C*N matrix of prices conditional on market m
        - P0_m: price of outside good in market m
    """

    # # define μR_icm based on input R
    μR_icm = R == "B" ? μB_icm : μC_icm   
    μR_0mm = R == "B" ? μB_0mm : μC_0mm 
    #PR_icm = R == "B" ? PB_icm : PC_icm
    #PR_0mm = R == "B" ? PB_0mm : PC_0mm

    # initial guess: P_m_g: C*N matrix
    P_m_g = [σ/(σ-1)*MC_icm(c,i,m,τ,exp_ω_m) for c in 1:C, i in 1:N]
    P0_m_g = σ/(σ-1)*MC_0mm(τ,exp_ω0_m)
    k = 1
    diff = 1.0
    while diff > TOL && k<MAXIT
        P_m_N = [μR_icm(c,i,P_m_g,exp_ξ_m,P0_m_g,exp_ξ0_m,σ)*MC_icm(c,i,m,τ,exp_ω_m) for c in 1:C, i in 1:N] 
        # P_m_N = [PR_icm(c,i,m, P_m_g, exp_ξ[:,:,m],exp_ω,P0_m_g,exp_ξ[m], τ,σ) for c in 1:C, i in 1:N]
        P0_m_N = μR_0mm(P_m_g,exp_ξ_m,P0_m_g,exp_ξ0_m,σ)*MC_0mm(τ,exp_ω0_m)
        # P0_m_N = PR_0mm(m, P_m_g,exp_ξ[:,:,m],P0_m_g,exp_ξ0[m],exp_ω0,τ,σ)
        diff = maximum(max.(abs.(P_m_N-P_m_g),abs(P0_m_N-P0_m_g)))
        P_m_g = w.*P_m_g .+ (1-w).*P_m_N
        P0_m_g = w*P0_m_g + (1-w)*P0_m_N
        #println("Iteration: ",k,"; Norm: ",diff)
        k += 1
    end   
    @assert k < MAXIT "No convergence after $(MAXIT) iterations"
   #  @assert all(P_m_g .> 0) "Negative prices found"

    return P_m_g, P0_m_g
end


function DGP_t(t::Int64, R::String, GP::global_param, global_seed::Int64, full::Bool)
    """
    Output the simulated data in array format - modified assertion error
    Input:
        - R: "B" (Bertrand) or "C" (Cournot)
        - t: task index
        - GP: global parameters
        - global_seed: global seed for random number generator
        - full: whether to output full equilibrium objects
    Output:
        - P: C*N*(C+M) array of prices
        - P0: C+M vector of prices of outside goods
        - RS: C*N*(C+M) array of revenue shares
        - RS0: C+M vector of revenue shares of outside goods
        - RV: N*C*M array of revenues (optional)
        - RV0: C+M vector of revenues of outside goods (optional)
        - Y: C+M vector of market sizes
        - exp_ξ: C*N*(C+M) array of orgin-firm-market specific demand shifters (optional)
        - exp_ω: C*N*(C+M) array of orgin-firm-market specific cost shifters (optional)
        - exp_ξ0: C+M vector of demand shifters for outside goods (optional)
        - exp_ω0: C+M vector of cost shifters for outside goods (optional)
        - τ: C*(C+M) matrix of trade costs
        - herfindahl: Herfindahl index
    """
    C = GP.C
    N = GP.N
    M = GP.M
    σ = GP.σ
    MAX_RETRIES = GP.MAX_RETRIES
    w = GP.w
    TOL = GP.TOL
    MAXIT = GP.MAXIT


    unique_seed = global_seed + t # Create a unique seed for this task
    task_rng = MersenneTwister(unique_seed)  # Initialize RNG with the task-specific seed
    # Y, exp_ξ, exp_ω, τ = set_θ(GP, task_rng)
    SDT = set_θ(GP, task_rng)
    exp_ξ = SDT.exp_ξ
    exp_ω = SDT.exp_ω
    exp_ξ0 = SDT.exp_ξ0
    exp_ω0 = SDT.exp_ω0
    τ = SDT.τ
    Y = SDT.Y
    # initialize P, RS, RV
    P = zeros(C, N, C+M)
    P0 = zeros(C+M)
    RS = zeros(C, N, C+M)
    RS0 = zeros(C+M)

    m = 1
    num_retry = 0

    while m <= C+M
        retries = 0
    while true
        try
            P[:,:,m], P0[m] = solve_P_m_R(R, m, C, N, exp_ξ[:,:,m], exp_ω[:,:,m], exp_ξ0[m], exp_ω0[m], τ, σ, w, MAXIT, TOL)
            break # Break the loop if the above line executes successfully
        catch e
            num_retry += 1
            if isa(e, AssertionError)
                retries += 1
                # println("AssertionError caught. Retry $retries...")
                if retries > MAX_RETRIES
                    error("Maximum retries reached. Unable to find a valid solution.")
                end
                unique_seed_new = global_seed + t + 1000*num_retry # Create a unique seed for this task
                task_rng_new = MersenneTwister(unique_seed_new)  # Initialize RNG with the task-specific seed
                # Y, exp_ξ, exp_ω, τ = set_θ(GP, task_rng)
                SDT = set_θ_update_m(m, SDT, GP, task_rng_new) # # update original SDT
                exp_ξ = SDT.exp_ξ
                exp_ω = SDT.exp_ω
                exp_ξ0 = SDT.exp_ξ0
                exp_ω0 = SDT.exp_ω0
                τ = SDT.τ
                Y = SDT.Y
            else
                rethrow()
            end
        end # end for try-catch
    end # end for while true
        RS[:,:,m] = [RS_icm(c,i,P[:,:,m],exp_ξ[:,:,m],P0[m],exp_ξ0[m], σ) for c in 1:C, i in 1:N]  
        RS0[m] = RS_0mm(P[:,:,m],exp_ξ[:,:,m],P0[m],exp_ξ0[m],σ)
        m += 1
    end
    # add
    herfindahl = sum(RS.^2) + sum(RS0.^2)
    
    if full
        RV = [RS[c,i,m]*Y[m] for c in 1:C, i in 1:N, m in 1:C+M]
        RV0 = RS0.*Y
        return eqbm_full(GP, R, P, P0, RS, RS0, Y, τ, herfindahl, num_retry, exp_ξ, exp_ω, exp_ξ0, exp_ω0, RV, RV0)
    else
        return eqbm(GP, R, P, P0, RS, RS0, Y, τ, herfindahl, num_retry)
    end
end



function DT_for_reg(eq::Union{eqbm, eqbm_full}, full::Bool)
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
    C = eq.GP.C
    N = eq.GP.N
    M = eq.GP.M
    R = eq.R
    P = eq.P
    P0 = eq.P0
    RS = eq.RS
    RS0 = eq.RS0
    τ = eq.τ
    σ = eq.GP.σ
    R = eq.R
    
    # inside good relative revenues (identical to relative shares)
    RV_normalized = zeros(C,N,C+M)
    P_relative = zeros(C,N,C+M)
    τ_relative = zeros(C,C+M)
    for m in 1:C+M
        RV_normalized[:,:,m] = RS[:,:,m]./RS0[m]   # change!
        P_relative[:,:,m] = P[:,:,m]./P0[m]
        τ_relative[:,m] = τ[:,m] # τ[m,m] = 1 in current setting
    end
    # reshape RV_normalized to a C*N*(C+M) vector
    RV_normalized_vec = reduce(vcat,reduce(vcat,permutedims(RV_normalized, [2,1,3])))
    RS_vec = reduce(vcat,reduce(vcat,permutedims(RS, [2,1,3]))) # absolute revenue shares!
    # reshape P_relative to a C*N*(C+M) vector
    P_vec = reduce(vcat,reduce(vcat,permutedims(P_relative, [2,1,3])))
    # reshape τ_relative to a C*(C+M) vector
    τ_vec = reduce(vcat,τ_relative)
    τ_long = repeat(τ_vec, inner = N) # C*N*(C+M) vector
    c = repeat(repeat(1:C, inner=N), outer=C + M)
    i = repeat(1:N, outer=C * (C + M))
    m = repeat(1:C+M, inner=N * C)
    RS0_long = repeat(RS0, inner = N*C)
    if full
        RV = eq.RV
        exp_ξ = eq.exp_ξ
        exp_ω = eq.exp_ω
        exp_ξ_true = reduce(vcat, reduce(vcat, permutedims(exp_ξ, [2, 1, 3])))
        ξ = log.(exp_ξ_true)
        exp_ω_true = reduce(vcat, reduce(vcat, permutedims(exp_ω, [2, 1, 3])))
        ω = log.(exp_ω_true)
        RV = reduce(vcat, reduce(vcat, permutedims(RV, [2, 1, 3])))
        exp_ξ0 = eq.exp_ξ0
        ξ0 = log.(exp_ξ0)
        ξ0_long = repeat(ξ0, inner = N*C)
        exp_ω0 = eq.exp_ω0
        ω0 = log.(exp_ω0)
        ω0_long = repeat(ω0, inner = N*C)
        μR_icm = R == "B" ? μB_icm : μC_icm     
        μ_array = [μR_icm(c,i,P[:,:,m],exp_ξ[:,:,m],P0[m],exp_ξ0[m],σ) for c in 1:C, i in 1:N, m in 1:C+M]    
        μR_0mm = R == "B" ? μB_0mm : μC_0mm
        μ0_array = [μR_0mm(P[:,:,m],exp_ξ[:,:,m],P0[m],exp_ξ0[m],σ) for m in 1:C+M]
        μ0_long = repeat(μ0_array, inner = N*C)
        μ_R_true = reduce(vcat, reduce(vcat, permutedims(μ_array, [2, 1, 3])))
        v_true = ξ .- ξ0_long .+ (1 - σ) * (ω .- ω0_long)
        P_true = reduce(vcat, reduce(vcat, permutedims(P, [2, 1, 3])))
        P0_vec = repeat(P0, inner = N*C)
        DT = DataFrame(c=c, i=i, m=m, R_normalized=RV_normalized_vec,RS=RS_vec,RS0 = RS0_long, P=P_vec, τ=τ_long, μ_R_true=μ_R_true, μ0 = μ0_long, ξ=ξ, ξ0 = ξ0_long, ω=ω,  ω0 = ω0_long, v_true=v_true, P_true = P_true, P0 = P0_vec) # P are relative prices; τ are relative τ
    else
        DT = DataFrame(c=c, i=i, m=m, R_normalized=RV_normalized_vec,RS=RS_vec, RS0 = RS0_long, P=P_vec, τ=τ_long)
    end
    return DT
end
