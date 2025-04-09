# given data, estimate the parameters. Different from MC_functions where the first step is to simulate data

Σ(x)= sum(x)
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

function simDT_for_reg(eq_input::eqbm_t, full::Bool)
    """
    Create a DataFrame for regression
    """
    C = eq_input.eqbm_m_list[1].IDT.C
    N = eq_input.eqbm_m_list[1].IDT.N
    M = length(eq_input.eqbm_m_list) - C
    R = eq_input.eqbm_m_list[1].R

    P = zeros(C,N,C+M)
    for m in 1:C+M
        P[:,:,m] = eq_input.eqbm_m_list[m].P_m
    end
    P0 = [eq_input.eqbm_m_list[m].P0_m for m in 1:C+M]
    RS = zeros(C,N,C+M)
    for m in 1:C+M
        RS[:,:,m] = eq_input.eqbm_m_list[m].RS_m
    end
    RS0 = [eq_input.eqbm_m_list[m].RS0_m for m in 1:C+M]
    τ = eq_input.eqbm_m_list[1].IDT.τ
    σ = eq_input.eqbm_m_list[1].IDT.σ

    t = eq_input.t

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
        RV = zeros(C,N,C+M)
        for m in 1:C+M
            RV[:,:,m] = eq_input.eqbm_m_list[m].RV_m
        end
        exp_ξ = zeros(C,N,C+M)
        for m in 1:C+M
            exp_ξ[:,:,m] = eq_input.eqbm_m_list[m].IDT.exp_ξ_m
        end
        exp_ω = zeros(C,N,C+M)
        for m in 1:C+M
            exp_ω[:,:,m] = eq_input.eqbm_m_list[m].IDT.exp_ω_m
        end
        exp_ξ_true = reduce(vcat, reduce(vcat, permutedims(exp_ξ, [2, 1, 3])))
        ξ = log.(exp_ξ_true)
        exp_ω_true = reduce(vcat, reduce(vcat, permutedims(exp_ω, [2, 1, 3])))
        ω = log.(exp_ω_true)
        RV = reduce(vcat, reduce(vcat, permutedims(RV, [2, 1, 3])))
        exp_ξ0 = [eq_input.eqbm_m_list[m].IDT.exp_ξ0_m for m in 1:C+M]
        ξ0 = log.(exp_ξ0)
        ξ0_long = repeat(ξ0, inner = N*C)
        exp_ω0 = [eq_input.eqbm_m_list[m].IDT.exp_ω0_m for m in 1:C+M]
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
        DT = DataFrame(t=t,c=c, i=i, m=m, R_normalized=RV_normalized_vec,RS=RS_vec,RS0 = RS0_long, P=P_vec, τ=τ_long, μ_R_true=μ_R_true, μ0 = μ0_long, ξ=ξ, ξ0 = ξ0_long, ω=ω,  ω0 = ω0_long, v_true=v_true, P_true = P_true, P0 = P0_vec) # P are relative prices; τ are relative τ
    else
        DT = DataFrame(t=t,c=c, i=i, m=m, R_normalized=RV_normalized_vec,RS=RS_vec, RS0 = RS0_long, P=P_vec, τ=τ_long)
    end
    return DT
end


function store_est(R::String, eq_input::eqbm_t, GP::global_param, full::Bool)
    """
    return a vector of estimates using different methods
    Input:
        - R: "B" (Bertrand) or "C" (Cournot)
        - eq_input: eqbm_t object, simulated data for a given t/task
        - global_seed: random seed
        - GP: global parameters
        - full: whether full versions of the estimates are needed
    Output:
        - dataframe of σD_OLS, σD_IV, σS_OLS, σSB_GMM, σSC_GMM, σDSB_GMM, σDSC_GMM, seSB_GMM, seSC_GMM, seDSB_GMM, seDSC_GMM, herfindahl
    """
    
    RS0 = [eq_input.eqbm_m_list[m].RS0_m for m in 1:GP.C+GP.M]
    mean_RS0 = mean(RS0)
    DT = simDT_for_reg(eq_input, full) # one period data
    σD_OLS = est_demand_OLS(DT)
    σD_IV = est_demand_IV(DT)
    σS_OLS = est_supply_OLS(DT)
    σSB_GMM, seSB_GMM = est_supply_GMM_twostep(GP,μB_icm_R,μB_0mm_R,DT)
    σSC_GMM, seSC_GMM = est_supply_GMM_twostep(GP,μC_icm_R,μC_0mm_R,DT)
    σDSB_GMM, seDSB_GMM = est_supply_demand_GMM_twostep(GP,μB_icm_R,μB_0mm_R,DT)
    σDSC_GMM, seDSC_GMM = est_supply_demand_GMM_twostep(GP,μC_icm_R,μC_0mm_R,DT)
    df = DataFrame(t=t, R=R, σD_OLS = σD_OLS, σD_IV = σD_IV, σS_OLS = σS_OLS, 
    σSB_GMM = σSB_GMM[1], σSC_GMM = σSC_GMM[1], σDSB_GMM = σDSB_GMM[1], σDSC_GMM = σDSC_GMM[1],
    seSB_GMM = seSB_GMM[1],seSC_GMM = seSC_GMM[1], seDSB_GMM=seDSB_GMM[1], seDSC_GMM=seDSC_GMM[1],
    herfindahl = eq.herfindahl, mean_RS0 = mean_RS0)
    return df
end


function simDT_bootstrap_single_iteration(t::Int64, iteration::Int64, GP::global_param, DT::DataFrame, size::Int64, global_seed::Int64, iteration_seed::Int64)
    """
    Return estimates for a single iteration in the bootstrap, for demand and supply GMM only
    """
    unique_seed = global_seed + 1000*t + iteration + iteration_seed # Create a unique seed for this task
    task_rng = MersenneTwister(unique_seed)  

    selected_market = sample(task_rng, 1:GP.M+GP.C, size, replace=true) # with replacement
    DT_temp = []
    for s in selected_market
        temp = filter(row -> row.m == s, DT)
        push!(DT_temp, temp)
    end
    DT_temp = vcat(DT_temp...)

    σD_IV = est_demand_IV(DT_temp)
    σDSB_GMM, seDSB_GMM = est_supply_demand_GMM_twostep(GP,μB_icm_R,μB_0mm_R,DT_temp,market_list_input = selected_market, seTRUE = false) # no need to calculate se
    σDSC_GMM, seDSC_GMM = est_supply_demand_GMM_twostep(GP,μC_icm_R,μC_0mm_R,DT_temp,market_list_input = selected_market, seTRUE = false) # no need to calculate se
    return  DataFrame(ΔDSB = σDSB_GMM[1] - σD_IV ,ΔDSC = σDSC_GMM[1] - σD_IV)
end



# var(itr; corrected::Bool=true, mean=nothing[, dims]) If corrected is true, then the sum is scaled with n-1, whereas the sum is scaled with n if corrected is false where n is the number of elements in itr.
function simDT_bootstrap_t(R::String, global_seed::Int64, GP::global_param, eq_input::eqbm_t, B::Int64, t::Int64, full::Bool, B_seed ::Int64; outputDT::Bool = false, folder_path::String = "")
    """
    Return hausman test for a single t, including DGP
    """

    RS0 = [eq_input.eqbm_m_list[m].RS0_m for m in 1:GP.C+GP.M]
    mean_RS0 = mean(RS0)
    herfindahl_list = [eq_input.eqbm_m_list[m].herfindahl for m in 1:GP.C+GP.M]
    herfindahl = sum(herfindahl_list)

    DT = simDT_for_reg(eq_input, full) # one period data
    if outputDT
        filename = folder_path*"temp_Conduct$(R)_C$(GP.C)_N$(GP.N)_M$(GP.M)_t$(t).csv"
        CSV.write(filename, DT)
    end

    σD_OLS = est_demand_OLS(DT)
    σD_IV = est_demand_IV(DT)
    σS_OLS = est_supply_OLS(DT)
    σSB_GMM, seSB_GMM = est_supply_GMM_twostep(GP,μB_icm_R,μB_0mm_R,DT)
    σSC_GMM, seSC_GMM = est_supply_GMM_twostep(GP,μC_icm_R,μC_0mm_R,DT)
    σDSB_GMM, seDSB_GMM = est_supply_demand_GMM_twostep(GP,μB_icm_R,μB_0mm_R,DT)
    σDSC_GMM, seDSC_GMM = est_supply_demand_GMM_twostep(GP,μC_icm_R,μC_0mm_R,DT)

    market_size = GP.C+GP.M # current setup
    # bootstrap with replacement
    # without parallelization
    result = []
    for iteration in (1+B_seed):(B+B_seed)
        retry = 0 
        while true
           try 
               push!(result, simDT_bootstrap_single_iteration(t, iteration, GP, DT, market_size, global_seed, retry*1000))
               break
           catch e
                retry += 1
                if retry > 100
                    println("Too many retries for a single iteration")
                    error("Too many retries for a single iteration")
                end
            end
        end
    #    push!(result, simDT_bootstrap_single_iteration(t, iteration, GP, DT, market_size, global_seed))
    end
    df_boot = vcat(result...)
    var_DSB = var(df_boot[!, :ΔDSB])
    var_DSC = var(df_boot[!, :ΔDSC])
    
    df = DataFrame(R=R, t=t, ΔDSB = σDSB_GMM[1] - σD_IV, ΔDSC = σDSC_GMM[1] - σD_IV, var_DSB = var_DSB, var_DSC = var_DSC, hausman_stat_DSB = (σDSB_GMM[1] - σD_IV)*inv(var_DSB)*(σDSB_GMM[1] - σD_IV), hausman_stat_DSC = (σDSC_GMM[1] - σD_IV)*inv(var_DSC)*(σDSC_GMM[1] - σD_IV), 
    σD_OLS = σD_OLS, σD_IV = σD_IV, σS_OLS = σS_OLS, 
    σSB_GMM = σSB_GMM[1], σSC_GMM = σSC_GMM[1], σDSB_GMM = σDSB_GMM[1], σDSC_GMM = σDSC_GMM[1],
    seSB_GMM = seSB_GMM[1],seSC_GMM = seSC_GMM[1], seDSB_GMM=seDSB_GMM[1], seDSC_GMM=seDSC_GMM[1],
    herfindahl = herfindahl, mean_RS0 = mean_RS0)
    
    # add: throw an error if var_DSB or var_DSC is too large
    # @assert var_DSB < 1e3 && var_DSC < 1e3
    println("period $t done.")
    return df
end



# var(itr; corrected::Bool=true, mean=nothing[, dims]) If corrected is true, then the sum is scaled with n-1, whereas the sum is scaled with n if corrected is false where n is the number of elements in itr.
function simDT_t(R::String, global_seed::Int64, GP::global_param, eq_input::eqbm_t, B::Int64, t::Int64, full::Bool; outputDT::Bool = false, folder_path::String = "")
    """
    estimation for a single t, including DGP
    """

    RS0 = [eq_input.eqbm_m_list[m].RS0_m for m in 1:GP.C+GP.M]
    mean_RS0 = mean(RS0)
    herfindahl_list = [eq_input.eqbm_m_list[m].herfindahl for m in 1:GP.C+GP.M]
    herfindahl = sum(herfindahl_list)

    DT = simDT_for_reg(eq_input, full) # one period data
    if outputDT
        filename = folder_path*"temp_Conduct$(R)_C$(GP.C)_N$(GP.N)_M$(GP.M)_t$(t).csv"
        CSV.write(filename, DT)
    end

    σD_OLS = est_demand_OLS(DT)
    σD_IV = est_demand_IV(DT)
    σS_OLS = est_supply_OLS(DT)
    σSB_GMM, seSB_GMM = est_supply_GMM_twostep(GP,μB_icm_R,μB_0mm_R,DT)
    σSC_GMM, seSC_GMM = est_supply_GMM_twostep(GP,μC_icm_R,μC_0mm_R,DT)
    σDSB_GMM, seDSB_GMM = est_supply_demand_GMM_twostep(GP,μB_icm_R,μB_0mm_R,DT)
    σDSC_GMM, seDSC_GMM = est_supply_demand_GMM_twostep(GP,μC_icm_R,μC_0mm_R,DT)

    
    df = DataFrame(R=R, t=t, ΔDSB = σDSB_GMM[1] - σD_IV, ΔDSC = σDSC_GMM[1] - σD_IV,
    σD_OLS = σD_OLS, σD_IV = σD_IV, σS_OLS = σS_OLS, 
    σSB_GMM = σSB_GMM[1], σSC_GMM = σSC_GMM[1], σDSB_GMM = σDSB_GMM[1], σDSC_GMM = σDSC_GMM[1],
    seSB_GMM = seSB_GMM[1],seSC_GMM = seSC_GMM[1], seDSB_GMM=seDSB_GMM[1], seDSC_GMM=seDSC_GMM[1],
    herfindahl = herfindahl, mean_RS0 = mean_RS0)

    return df
end


function simDT_t_modified(R::String, global_seed::Int64, GP::global_param, eq_input::eqbm_t, B::Int64, t::Int64, full::Bool; outputDT::Bool = false, folder_path::String = "")
    """
    estimation for a single t, including DGP
    """

    RS0 = [eq_input.eqbm_m_list[m].RS0_m for m in 1:GP.C+GP.M]
    mean_RS0 = mean(RS0)
    herfindahl_list = [eq_input.eqbm_m_list[m].herfindahl for m in 1:GP.C+GP.M]
    herfindahl = sum(herfindahl_list)

    DT = simDT_for_reg(eq_input, full) # one period data
    if outputDT
        filename = folder_path*"temp_Conduct$(R)_C$(GP.C)_N$(GP.N)_M$(GP.M)_t$(t).csv"
        CSV.write(filename, DT)
    end

    σD_OLS = est_demand_OLS(DT)
    σD_IV = est_demand_IV(DT)
    σS_OLS = est_supply_OLS(DT)
    σSB_GMM, seSB_GMM = est_supply_GMM_twostep(GP,μB_icm_R,μB_0mm_R,DT)
    σSC_GMM, seSC_GMM = est_supply_GMM_twostep(GP,μC_icm_R,μC_0mm_R,DT)
    σDSB_GMM, seDSB_GMM = est_supply_demand_GMM_twostep_modified(GP,μB_icm_R,μB_0mm_R,DT)
    σDSC_GMM, seDSC_GMM = est_supply_demand_GMM_twostep_modified(GP,μC_icm_R,μC_0mm_R,DT)

    
    df = DataFrame(R=R, t=t, ΔDSB = σDSB_GMM[1] - σD_IV, ΔDSC = σDSC_GMM[1] - σD_IV,
    σD_OLS = σD_OLS, σD_IV = σD_IV, σS_OLS = σS_OLS, 
    σSB_GMM = σSB_GMM[1], σSC_GMM = σSC_GMM[1], σDSB_GMM = σDSB_GMM[1], σDSC_GMM = σDSC_GMM[1],
    seSB_GMM = seSB_GMM[1],seSC_GMM = seSC_GMM[1], seDSB_GMM=seDSB_GMM[1], seDSC_GMM=seDSC_GMM[1],
    herfindahl = herfindahl, mean_RS0 = mean_RS0)

    return df
end



function simDT_bootstrap_MC_safe(R::String, GP::global_param, eq_input::eqbm_output, B::Int64; full = false)
    """
    Monte-Carlo simulation, returns Hausman stats etc. Safe version. Deal with unexpected errors like "ArgumentError: matrix contains Infs or NaNs"
        - GP: global parameters
        - global_seed: two random seeds
        - T: number of simulations
        - B: number of bootstrap iterations
        - full: whether full versions of the estimates are needed
    Output:
        - dfB: a dataframe of statistics for DGP based on Bertrand
        - dfC: a dataframe of statistics for DGP based on Cournot
    """ 

    function try_bootstrap(t,R,B_seed)
        try
            return (idx=t, success=true, result = simDT_bootstrap_t(R, seed, GP, eqbm_t_list[t], B, t, full, B_seed))
        catch e
            println("Error in task $t: $e")
            return (idx=t, success=false, result=nothing)
        end
    end  
    eqbm_t_list = eq_input.eqbm_t_list
    T = length(eqbm_t_list)

    results = ThreadsX.map(t -> try_bootstrap(t, R, 0), 1:T)
    # Filter successful and unsuccessful results
    successes = filter(r -> r.success, results)
    df = reduce(vcat, [r.result for r in successes], cols=:union)

    # # Add: remove rows with large variances
    # df = filter(row -> row.var_DSB < 1e3 && row.var_DSC < 1e3, df)

    # count number of successful tasks
    n_success = size(df, 1)
    times = 1
    while(n_success < T && times <=50)
        println("Retry for $T - $n_success unsuccessful tasks")
        # retry for unsuccessful tasks
        unsuccessful = filter(r -> !r.success, results)
        # update 
        results = ThreadsX.map(r -> try_bootstrap(r.idx, R, times*1000*B), unsuccessful)
        successes = filter(r -> r.success, results)
        if size(successes, 1) > 0
            df_temp = reduce(vcat, [r.result for r in successes], cols=:union)
            df = vcat(df, df_temp)
            n_success = size(df, 1)    
        end
        times += 1
    end

    return df # only return successful results
end



function simDT_MC_safe(R::String, GP::global_param, eq_input::eqbm_output, B::Int64; full = false)
    """
    Monte-Carlo simulation, safe version. Deal with unexpected errors like "ArgumentError: matrix contains Infs or NaNs"
        - GP: global parameters
        - global_seed: two random seeds
        - T: number of simulations
        - B: number of bootstrap iterations
        - full: whether full versions of the estimates are needed
    Output:
        - dfB: a dataframe of statistics for DGP based on Bertrand
        - dfC: a dataframe of statistics for DGP based on Cournot
    """ 

    function try_est(t, R)
        try
            return (idx=t, success=true, result = simDT_t(R, seed, GP, eqbm_t_list[t], B, t, full))
        catch e
            println("Error in estimation task $t: $e")
            return (idx=t, success=false, result=nothing)
        end
    end  
    eqbm_t_list = eq_input.eqbm_t_list
    T = length(eqbm_t_list)

    results = ThreadsX.map(t -> try_est(t, R), 1:T)
    # Filter successful and unsuccessful results
    successes = filter(r -> r.success, results)
    df = reduce(vcat, [r.result for r in successes], cols=:union)

    # # Add: remove rows with large variances
    # df = filter(row -> row.var_DSB < 1e3 && row.var_DSC < 1e3, df)

    return df # only return successful results
end




# function simDT_MC_safe_modified(R::String, GP::global_param, eq_input::eqbm_output, B::Int64; full = false)
#     """
#     Monte-Carlo simulation, safe version. Deal with unexpected errors like "ArgumentError: matrix contains Infs or NaNs"
#         - GP: global parameters
#         - global_seed: two random seeds
#         - T: number of simulations
#         - B: number of bootstrap iterations
#         - full: whether full versions of the estimates are needed
#     Output:
#         - dfB: a dataframe of statistics for DGP based on Bertrand
#         - dfC: a dataframe of statistics for DGP based on Cournot
#     """ 

#     function try_est(t, R)
#         try
#             return (idx=t, success=true, result = simDT_t_modified(R, seed, GP, eqbm_t_list[t], B, t, full))
#         catch e
#             println("Error in estimation task $t: $e")
#             return (idx=t, success=false, result=nothing)
#         end
#     end  
#     eqbm_t_list = eq_input.eqbm_t_list
#     T = length(eqbm_t_list)

#     results = ThreadsX.map(t -> try_est(t, R), 1:T)
#     # Filter successful and unsuccessful results
#     successes = filter(r -> r.success, results)
#     df = reduce(vcat, [r.result for r in successes], cols=:union)

#     # # Add: remove rows with large variances
#     # df = filter(row -> row.var_DSB < 1e3 && row.var_DSC < 1e3, df)

#     return df # only return successful results
# end