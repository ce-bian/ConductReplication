# 2024.6.23: need to update the code to keep track of simulated data and output it

# function store_est(R::String, t::Int64, GP::global_param, global_seed::Int64, full::Bool)
#     """
#     return a vector of estimates using different methods
#     Input:
#         - R: "B" (Bertrand) or "C" (Cournot)
#         - t: simulation index
#         - global_seed: random seed
#         - GP: global parameters
#         - full: whether full versions of the estimates are needed
#     Output:
#         - dataframe of σD_OLS, σD_IV, σS_OLS, σSB_GMM, σSC_GMM, σDSB_GMM, σDSC_GMM, seSB_GMM, seSC_GMM, seDSB_GMM, seDSC_GMM, herfindahl
#     """
#     # P, RS,RV, Y, exp_ξ, exp_ω, τ = DGP(R,σ_true)
#     eq = DGP_t(t, R,GP, global_seed, full) # modified here
#     mean_RS0 = mean(eq.RS0)
#     DT = DT_for_reg(eq, full) # one period data
#     σD_OLS = est_demand_OLS(DT)
#     σD_IV = est_demand_IV(DT)
#     σS_OLS = est_supply_OLS(DT)
#     σSB_GMM, seSB_GMM = est_supply_GMM_twostep(GP,μB_icm_R,μB_0mm_R,DT)
#     σSC_GMM, seSC_GMM = est_supply_GMM_twostep(GP,μC_icm_R,μC_0mm_R,DT)
#     σDSB_GMM, seDSB_GMM = est_supply_demand_GMM_twostep(GP,μB_icm_R,μB_0mm_R,DT)
#     σDSC_GMM, seDSC_GMM = est_supply_demand_GMM_twostep(GP,μC_icm_R,μC_0mm_R,DT)
#     df = DataFrame(t=t,σD_OLS = σD_OLS, σD_IV = σD_IV, σS_OLS = σS_OLS, 
#     σSB_GMM = σSB_GMM[1], σSC_GMM = σSC_GMM[1], σDSB_GMM = σDSB_GMM[1], σDSC_GMM = σDSC_GMM[1],
#     seSB_GMM = seSB_GMM[1],seSC_GMM = seSC_GMM[1], seDSB_GMM=seDSB_GMM[1], seDSC_GMM=seDSC_GMM[1],
#     herfindahl = eq.herfindahl, mean_RS0 = mean_RS0)
#     return df
# end


# function MC_simulation(GP::global_param, global_seed::Int64, T::Int64, full::Bool)
#     """
#     Monte-carlo simulation, returns the estimated sigma and se
#     Input:
#         - GP: global parameters
#         - global_seed: random seed
#         - T: number of simulations
#         - full: whether full versions of the estimates are needed
#     Output:
#         - dfB: a dataframe of estimates for DGP based on Bertrand
#         - dfC: a dataframe of estimates for DGP based on Cournot
#     """
#      result = ThreadsX.map(t->store_est("B",t, GP, global_seed, full),1:T)
#      dfB = reduce(vcat,result,cols=:union)
#      result = ThreadsX.map(t->store_est("C",t, GP, global_seed, full),1:T)
#      dfC = reduce(vcat,result,cols=:union)
#     return dfB, dfC
# end


# function MC_simulation_twoseeds(GP::global_param, global_seed::Vector{Int64}, T::Int64, full::Bool)
#     """
#     Monte-carlo simulation, returns the estimated sigma and se
#     Input:
#         - GP: global parameters
#         - global_seed: two random seeds
#         - T: number of simulations
#         - full: whether full versions of the estimates are needed
#     Output:
#         - dfB: a dataframe of estimates for DGP based on Bertrand
#         - dfC: a dataframe of estimates for DGP based on Cournot
#     """
#      result = ThreadsX.map(t->store_est("B",t, GP, global_seed[1], full),1:T)
#      dfB = reduce(vcat,result,cols=:union)
#      result = ThreadsX.map(t->store_est("C",t, GP, global_seed[2], full),1:T)
#      dfC = reduce(vcat,result,cols=:union)
#     return dfB, dfC
# end


# function sim_dt_f(R::String, t::Int64, GP::global_param, global_seed::Int64, full::Bool)
#     """
#     Return a dataframe of simulated data
#     """
#     eq = DGP_t(t, R, GP, global_seed, full)
#     DT = DT_for_reg(eq, full)
#     # add column t
#     DT[!, :t] .= t
#     return DT
# end



function bootstrap_single_iteration(t::Int64, iteration::Int64, GP::global_param, DT::DataFrame, size::Int64, global_seed::Int64)
    """
    Return estimates for a single iteration in the bootstrap, for demand and supply GMM only
    """
    unique_seed = global_seed + 1000*t + iteration # Create a unique seed for this task
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
function bootstrap_t(R::String, t::Int64, GP::global_param, B::Int64, global_seed::Int64, full::Bool; folder_path::String = "")
    """
    Return hausman test for a single t, including DGP
    """
    eq = DGP_t(t, R, GP, global_seed, full) 
    mean_RS0 = mean(eq.RS0)
    DT = DT_for_reg(eq, full) # one period data
    if full
        # add t to DT
        DT[!, :t] .= t
        # modify later
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
    for iteration in 1:B
       push!(result, bootstrap_single_iteration(t, iteration, GP, DT, market_size, global_seed))
    end
    df_boot = vcat(result...)
    var_DSB = var(df_boot[!, :ΔDSB])
    var_DSC = var(df_boot[!, :ΔDSC])
    
    df = DataFrame(t=t, ΔDSB = σDSB_GMM[1] - σD_IV, ΔDSC = σDSC_GMM[1] - σD_IV, var_DSB = var_DSB, var_DSC = var_DSC, hausman_stat_DSB = (σDSB_GMM[1] - σD_IV)*inv(var_DSB)*(σDSB_GMM[1] - σD_IV), hausman_stat_DSC = (σDSC_GMM[1] - σD_IV)*inv(var_DSC)*(σDSC_GMM[1] - σD_IV), 
    σD_OLS = σD_OLS, σD_IV = σD_IV, σS_OLS = σS_OLS, 
    σSB_GMM = σSB_GMM[1], σSC_GMM = σSC_GMM[1], σDSB_GMM = σDSB_GMM[1], σDSC_GMM = σDSC_GMM[1],
    seSB_GMM = seSB_GMM[1],seSC_GMM = seSC_GMM[1], seDSB_GMM=seDSB_GMM[1], seDSC_GMM=seDSC_GMM[1],
    herfindahl = eq.herfindahl, mean_RS0 = mean_RS0)

    return df
end


function bootstrap_MC_safe(GP::global_param, B::Int64, T::Int64, global_seed::Int64; full = false)
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
    function try_bootstrap(t, type)
        try
            return (idx=t, success=true, result=bootstrap_t(type, t, GP, B, global_seed, full))
        catch e
            return (idx=t, success=false, result=nothing)
        end
    end  
    successful_results = []
    remaining_indices = 1:T  # Start with all indices
    while length(successful_results) < T
        results = ThreadsX.map(t -> try_bootstrap(t, "B"), remaining_indices)
        # Filter successful and unsuccessful results
        successes = filter(r -> r.success, results)
        unsuccessful = filter(r -> !r.success, results)
        append!(successful_results, successes)
        # If there are more results needed, prepare for another round
        if length(successful_results) < T
            # Prepare new indices for the next batch
            remaining_indices = [r.idx + T for r in unsuccessful]  
        end
    end
    dfB = reduce(vcat, [r.result for r in successful_results], cols=:union)

    # Repeat for type "C"
    successful_results = []
    remaining_indices = 1:T  # Reset for "C"
    while length(successful_results) < T
        results = ThreadsX.map(t -> try_bootstrap(t, "C"), remaining_indices)
        # Filter successful and unsuccessful results
        successes = filter(r -> r.success, results)
        unsuccessful = filter(r -> !r.success, results)
        append!(successful_results, successes)
        if length(successful_results) < T
            remaining_indices = [r.idx + T for r in unsuccessful]  
        end
    end
    dfC = reduce(vcat, [r.result for r in successful_results], cols=:union)
    return dfB, dfC
end
