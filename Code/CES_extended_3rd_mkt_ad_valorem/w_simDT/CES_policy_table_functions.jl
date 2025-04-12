# version with full grid dataframe, fill in NaN  
# old functions are moved to Old/CES_policy_table_functions_20240527.jl

# need modification -> should use demand IV estimate as the benchmark??? UPDATED (2024.05.28)
function inferred_conduct_t_hausman_f(critical_value:: Float64, σD_IV, σDSB_GMM::Float64, σDSC_GMM::Float64, hausman_stat_DSB::Float64, hausman_stat_DSC::Float64; close = false)
    """
    Determine the optimal conduct for a given t based on the Hausman test statistics
    Input:
        - σD_IV: demand IV estimate of σ
        - σDSB_GMM: GMM estimate of σ for Bertrand
        - σDSC_GMM: GMM estimate of σ for Cournot
        - hausman_stat_DSB: Hausman statistic for Bertrand
        - hausman_stat_DSC: Hausman statistic for Cournot
        - close: if true, return the closest estimate to the true σ if both hausman tests are insignificant
    Output:
        - optimal_conduct: "B" or "C"
    """
    if hausman_stat_DSB > critical_value && hausman_stat_DSC < critical_value # reject H0 for Bertrand, fail to reject H0 for Cournot
        return "C"
    elseif hausman_stat_DSB < critical_value && hausman_stat_DSC > critical_value # reject H0 for Cournot, fail to reject H0 for Bertrand
        return "B"
    elseif hausman_stat_DSB < critical_value && hausman_stat_DSC < critical_value # fail to reject H0 for both
        if close
            if abs(σDSB_GMM - σD_IV) < abs(σDSC_GMM - σD_IV)
                return "B"
            else
                return "C"
            end
        else 
            return "N"
        end
    else # reject both, so can't determine
        return "N"
    end
end


function inferred_conduct_t_closest_f(σD_IV::Float64, σDSB_GMM::Float64, σDSC_GMM::Float64)
    """
    Determine the optimal conduct for a given t based on the estimated GMM σ, see which one is closer to the demand IV σ
    Input:
        - σD_IV: demand IV estimate of σ
        - σDSB_GMM: GMM estimate of σ for Bertrand
        - σDSC_GMM: GMM estimate of σ for Cournot
    Output:
        - optimal_conduct: "B" or "C"
    """
    DSB_diff = abs(σDSB_GMM - σD_IV)
    DSC_diff = abs(σDSC_GMM - σD_IV)
    if DSB_diff <= DSC_diff
        return "B"
    else
        return "C"
    end
end


function optimal_s_inferred_true_full_t(R_true::String, R_inferred::String, eqbm_input::eqbm_t, eqbm_true::eqbm_t, c::Int64, C::Int64, M::Int64, t::Int64)
    """
    Determine the optimal s and calculate welfare for a given t, both for R_true and R_inferred
    Input:
        - R_true: true conduct: "B" (Bertrand) or "C" (Cournot)
        - R: inferred conduct: "B" (Bertrand) or "C" (Cournot) or "N" (None)
        -  
    Output: 
        - W_df_m: DataFrame containing the optimal s for each third market
    """
    eqbm_m_list = eqbm_input.eqbm_m_list
    S_m_list = Array{Float64,1}(undef, M)  # M third markets 
    W_0_list = Array{Float64,1}(undef, M)  # M third markets
    W_optimal_m_list = Array{Float64,1}(undef, M)  # M third markets
    delta_m_list = Array{Float64,1}(undef, M)  # M third markets
    pct_m_list = Array{Float64,1}(undef, M)  # M third markets
    
    # eqba_m_true_R_list = eqbm_true_R_input.eqbm_m_list
    S_m_true_list = Array{Float64,1}(undef, M)  # M third markets 
    W_0_true_list = Array{Float64,1}(undef, M)  # M third markets, should be the same as W_0_list
    W_optimal_m_true_list = Array{Float64,1}(undef, M)  # M third markets
    delta_m_true_list = Array{Float64,1}(undef, M)  # M third markets
    pct_m_true_list = Array{Float64,1}(undef, M)  # M third markets

    i = 1
    for m in C+1:size(eqbm_m_list, 1)  # C+M
        eq_m = eqbm_m_list[m]
        IDT_m = eq_m.IDT
        # eq_true_R_m = eqba_m_true_R_list[m]
        # IDT_true_R_m = eq_true_R_m.IDT
        eq_true_m = eqbm_true.eqbm_m_list[m]
        IDT_true_m = eq_true_m.IDT

        W_0_list[i] = W_s_c_given_eq(R_true, eq_true_m)

        if R_inferred == "B"
            S_m_list[i] = Bertrand_optimal_s_m_c_f(IDT_m, c)
            # check if the Cournot optimal s is NaN, i.e. no convergence
            if isnan(S_m_list[i])
                W_optimal_m_list[i] = NaN
                delta_m_list[i] = NaN
                pct_m_list[i] = NaN
            else
                W_optimal_m_list[i] = W_s_c(R_true, IDT_true_m, S_m_list[i]) # use IDT_true_m (underlying information)
                delta_m_list[i] = W_optimal_m_list[i] - W_0_list[i]
                pct_m_list[i] = (W_optimal_m_list[i] - W_0_list[i]) / W_0_list[i]
            end
        elseif R_inferred == "C"
            S_m_list[i] = Cournot_optimal_s_m_c_f(IDT_m, c)
            # check if the Cournot optimal s is NaN, i.e. no convergence
            if isnan(S_m_list[i])
                W_optimal_m_list[i] = NaN
                delta_m_list[i] = NaN
                pct_m_list[i] = NaN
            else
                W_optimal_m_list[i] = W_s_c(R_true, IDT_true_m, S_m_list[i]) # use IDT_true_m (underlying information)
                delta_m_list[i] = W_optimal_m_list[i] - W_0_list[i]
                pct_m_list[i] = (W_optimal_m_list[i] - W_0_list[i]) / W_0_list[i]
            end
        elseif R_inferred == "N"
            # change here
            S_m_list[i] = 0.0  # s = 0 since can't determine the conduct
            W_optimal_m_list[i] = W_0_list[i] 
            delta_m_list[i] = 0.0
            pct_m_list[i] = 0.0 # no change
        end
        

        # this part should be independent of the inferred conduct
        eq_true_m = eqbm_true.eqbm_m_list[m]
        IDT_true_m = eq_true_m.IDT
        W_0_true_list[i] = W_s_c_given_eq(R_true, eq_true_m) # same as W_0_list[i]
        if R_true == "B" 
            # if the true conduct is Bertrand
            S_m_true_list[i] = Bertrand_optimal_s_m_c_f(IDT_true_m, c)
            # check if the Bertrand optimal s is NaN, i.e. no convergence
            if isnan(S_m_true_list[i])
                W_optimal_m_true_list[i] = NaN
                delta_m_true_list[i] = NaN
                pct_m_true_list[i] = NaN
            else
                W_optimal_m_true_list[i] = W_s_c(R_true, IDT_true_m, S_m_true_list[i])
                delta_m_true_list[i] = W_optimal_m_true_list[i] - W_0_true_list[i]
                pct_m_true_list[i] = (W_optimal_m_true_list[i] - W_0_true_list[i]) / W_0_true_list[i]
            end
        elseif R_true == "C" 
            # if the true conduct is Cournot
            S_m_true_list[i] = Cournot_optimal_s_m_c_f(IDT_true_m, c)
            # check if the Cournot optimal s is NaN
            if isnan(S_m_true_list[i])
                W_optimal_m_true_list[i] = NaN
                delta_m_true_list[i] = NaN
                pct_m_true_list[i] = NaN
            else
                W_optimal_m_true_list[i] = W_s_c(R_true, IDT_true_m, S_m_true_list[i])
                delta_m_true_list[i] = W_optimal_m_true_list[i] - W_0_true_list[i]
                pct_m_true_list[i] = (W_optimal_m_true_list[i] - W_0_true_list[i]) / W_0_true_list[i]
            end
        end
        i += 1
    end
    W_df_m = DataFrame(true_conduct = R_true, inferred_conduct = R_inferred, t = t, c = c, m = C+1:C+M, 
    S_optimal_m = S_m_list, W0_m = W_0_list, Woptimal_m = W_optimal_m_list, deltaW_m = delta_m_list, pctW_m = pct_m_list, 
    S_optimal_true_m = S_m_true_list, W0_true_m = W_0_true_list, Woptimal_true_m = W_optimal_m_true_list, deltaW_true_m = delta_m_true_list, pctW_true_m = pct_m_true_list)

    # check the size of the dataframe
    if size(W_df_m, 1) != M
        error("Error: the size of the dataframe is not M")
    end

    return W_df_m    
end





# key, need to return a full grid dataframe, fill in NaN, modified to parallel version
function W_hausman_inferred_true_varyM_full_f(InputDTFile::String, InputEstFile::String, critical_value::Float64, c::Int64, selected_M::Int64; Output = false, OutputFile::String = "")
    dfB, dfC, GP = load(InputEstFile, "dfB", "dfC", "GP")
    eqbm_output_B, eqbm_output_C, _ = load(InputDTFile, "eqbm_output_B", "eqbm_output_C", "GP")
    eqbm_subset_B = sample_mkt_f(eqbm_output_B, selected_M)
    eqbm_subset_C = sample_mkt_f(eqbm_output_C, selected_M)
    eqbm_t_list_B = eqbm_subset_B.eqbm_t_list
    eqbm_t_list_C = eqbm_subset_C.eqbm_t_list

    t_list_B = unique(dfB.t)
    t_list_C = unique(dfC.t)
    T_length_B = length(t_list_B)
    T_length_C = length(t_list_C)
    T_dict_B = Dict(i => t_list_B[i] for i in 1:T_length_B)
    T_dict_C = Dict(i => t_list_C[i] for i in 1:T_length_C)

    W_dfB_chunks = Vector{Any}(undef, T_length_B)
    dfB_inferred_conduct_data = Vector{Tuple{String, Int64, String}}(undef, T_length_B)
    
    pb_B = Progress(T_length_B, desc="Processing B", dt=0.5)

    Threads.@threads for t in 1:T_length_B
        dfB_t = @subset(dfB, :t .== T_dict_B[t])
        inferred = inferred_conduct_t_hausman_f(critical_value, dfB_t.σD_IV[1], dfB_t.σDSB_GMM[1], dfB_t.σDSC_GMM[1], dfB_t.hausman_stat_DSB[1], dfB_t.hausman_stat_DSC[1])
        dfB_inferred_conduct_data[t] = ("B", T_dict_B[t], inferred)
        est_eqbm_t = if inferred == "B"
            estimated_eqbm_t_f("B", eqbm_t_list_B[t], dfB_t.σDSB_GMM[1])
        elseif inferred == "C"
            estimated_eqbm_t_f("C", eqbm_t_list_B[t], dfB_t.σDSC_GMM[1])
        else
            eqbm_t_list_B[t]
        end
        W_df_m = optimal_s_inferred_true_full_t("B", inferred, est_eqbm_t, eqbm_t_list_B[t], c, GP.C, selected_M, T_dict_B[t])
        W_dfB_chunks[t] = W_df_m
        next!(pb_B)
    end

    W_dfC_chunks = Vector{Any}(undef, T_length_C)
    dfC_inferred_conduct_data = Vector{Tuple{String, Int64, String}}(undef, T_length_C)
    
    pb_C = Progress(T_length_C, desc="Processing C", dt=0.5)

    Threads.@threads for t in 1:T_length_C
        dfC_t = @subset(dfC, :t .== T_dict_C[t])
        inferred = inferred_conduct_t_hausman_f(critical_value, dfC_t.σD_IV[1], dfC_t.σDSB_GMM[1], dfC_t.σDSC_GMM[1], dfC_t.hausman_stat_DSB[1], dfC_t.hausman_stat_DSC[1])
        dfC_inferred_conduct_data[t] = ("C", T_dict_C[t], inferred)
        est_eqbm_t = if inferred == "B"
            estimated_eqbm_t_f("B", eqbm_t_list_C[t], dfC_t.σDSB_GMM[1])
        elseif inferred == "C"
            estimated_eqbm_t_f("C", eqbm_t_list_C[t], dfC_t.σDSC_GMM[1])
        else
            eqbm_t_list_C[t]
        end
        W_df_m = optimal_s_inferred_true_full_t("C", inferred, est_eqbm_t, eqbm_t_list_C[t], c, GP.C, selected_M, T_dict_C[t])
        W_dfC_chunks[t] = W_df_m
        next!(pb_C)
    end

    dfB_inferred_conduct = DataFrame(DGP = [x[1] for x in dfB_inferred_conduct_data],
                                     t = [x[2] for x in dfB_inferred_conduct_data],
                                     optimal_conduct = [x[3] for x in dfB_inferred_conduct_data])
    dfC_inferred_conduct = DataFrame(DGP = [x[1] for x in dfC_inferred_conduct_data],
                                     t = [x[2] for x in dfC_inferred_conduct_data],
                                     optimal_conduct = [x[3] for x in dfC_inferred_conduct_data])

    W_dfB = vcat(W_dfB_chunks...)
    W_dfC = vcat(W_dfC_chunks...)

    if Output
        save(OutputFile, "dfB_inferred_conduct", dfB_inferred_conduct, "dfC_inferred_conduct", dfC_inferred_conduct, "W_dfB", W_dfB, "W_dfC", W_dfC)
    end

    if size(W_dfB, 1) != selected_M * T_length_B
        error("Error: the size of the dataframe W_dfB is not selected_M * T_length_B")
    end
    if size(W_dfC, 1) != selected_M * T_length_C
        error("Error: the size of the dataframe W_dfC is not selected_M * T_length_C")
    end

    return W_dfB, W_dfC
end





# key, need to return a full grid dataframe, fill in NaN
function W_hausman_inferred_R_only_f( InputEstFile::String, critical_value::Float64, selected_M::Int64)
    """
    Determine the optimal conduct for each t based on the Hausman test statistics
    Input:
        - InputDTFile: the file containing the simulated data
        - InputEstFile: the file containing the estimates
        - critical_value: critical value for the Hausman test
    Output:
        - dfB_optimal_conduct: DataFrame containing the optimal conduct for Bertrand
        - dfC_optimal_conduct: DataFrame containing the optimal conduct for Cournot 
    """
    dfB, dfC, GP = load(InputEstFile, "dfB", "dfC", "GP")
    # get the unique t_list 
    t_list_B = unique(dfB.t)
    t_list_C = unique(dfC.t)
    T_length_B = length(t_list_B)
    T_length_C = length(t_list_C)
    T_dict_B = Dict(i => t_list_B[i] for i in 1:T_length_B)
    T_dict_C = Dict(i => t_list_C[i] for i in 1:T_length_C)
    # # sort T_dict_C
    # T_dict_C = sort(T_dict_C)
    # T_dict_B = sort(T_dict_B)
    dfB_inferred_conduct = DataFrame(DGP = String[], t = Int64[], M = Int64[], optimal_conduct = String[], hausman_stat_DSB = Float64[], hausman_stat_DSC = Float64[], critical_value = Float64[])
    dfC_inferred_conduct = DataFrame(DGP = String[], t = Int64[], 
    M= Int64[], optimal_conduct = String[], hausman_stat_DSB = Float64[], hausman_stat_DSC = Float64[], critical_value = Float64[])

    for t in 1:T_length_B
        dfB_t = @subset(dfB, :t .== T_dict_B[t])
        inferred_conduct_t = inferred_conduct_t_hausman_f(critical_value, dfB_t.σD_IV[1], dfB_t.σDSB_GMM[1], dfB_t.σDSC_GMM[1], dfB_t.hausman_stat_DSB[1], dfB_t.hausman_stat_DSC[1])
        push!(dfB_inferred_conduct, ("B", T_dict_B[t], selected_M, inferred_conduct_t, dfB_t.hausman_stat_DSB[1], dfB_t.hausman_stat_DSC[1], critical_value))
    end
    for t in 1:T_length_C
        dfC_t = @subset(dfC, :t .== T_dict_C[t])
        inferred_conduct_t = inferred_conduct_t_hausman_f(critical_value, dfC_t.σD_IV[1], dfC_t.σDSB_GMM[1], dfC_t.σDSC_GMM[1], dfC_t.hausman_stat_DSB[1], dfC_t.hausman_stat_DSC[1])
        push!(dfC_inferred_conduct, ("C", T_dict_C[t], selected_M, inferred_conduct_t, dfC_t.hausman_stat_DSB[1], dfC_t.hausman_stat_DSC[1], critical_value))
    end


    # Check the size of the dataframes
    # if size(dfB_inferred_conduct, 1) != selected_M * T_length_B
    #     error("Error: the size of the dataframe W_dfB is not selected_M * T_length_B")
    # end
    # if size(dfC_inferred_conduct, 1) != selected_M * T_length_C
    #     error("Error: the size of the dataframe W_dfC is not selected_M * T_length_C")
    # end

    return dfB_inferred_conduct, dfC_inferred_conduct
end



function W_closest_inferred_true_varyM_full_f(InputDTFile::String, InputEstFile::String, c::Int64, selected_M::Int64; Output = false, OutputFile::String = "")
    dfB, dfC, GP = load(InputEstFile, "dfB", "dfC", "GP")
    eqbm_output_B, eqbm_output_C, _ = load(InputDTFile, "eqbm_output_B", "eqbm_output_C", "GP")
    eqbm_subset_B = sample_mkt_f(eqbm_output_B, selected_M)
    eqbm_subset_C = sample_mkt_f(eqbm_output_C, selected_M)
    eqbm_t_list_B = eqbm_subset_B.eqbm_t_list
    eqbm_t_list_C = eqbm_subset_C.eqbm_t_list

    t_list_B = unique(dfB.t)
    t_list_C = unique(dfC.t)
    T_length_B = length(t_list_B)
    T_length_C = length(t_list_C)
    T_dict_B = Dict(i => t_list_B[i] for i in 1:T_length_B)
    T_dict_C = Dict(i => t_list_C[i] for i in 1:T_length_C)

    W_dfB_chunks = Vector{Any}(undef, T_length_B)
    dfB_inferred_conduct_data = Vector{Tuple{String, Int64, String}}(undef, T_length_B)

    pb_B = Progress(T_length_B, desc = "Processing B", dt = 0.5)

    Threads.@threads for t in 1:T_length_B
        dfB_t = @subset(dfB, :t .== T_dict_B[t])
        inferred = inferred_conduct_t_closest_f(dfB_t.σD_IV[1], dfB_t.σDSB_GMM[1], dfB_t.σDSC_GMM[1])
        dfB_inferred_conduct_data[t] = ("B", T_dict_B[t], inferred)
        est_eqbm_t = if inferred == "B"
            estimated_eqbm_t_f("B", eqbm_t_list_B[t], dfB_t.σDSB_GMM[1])
        elseif inferred == "C"
            estimated_eqbm_t_f("C", eqbm_t_list_B[t], dfB_t.σDSC_GMM[1])
        else
            eqbm_t_list_B[t]  # use true equilibrium, doesn't matter
        end
        W_dfB_chunks[t] = optimal_s_inferred_true_full_t("B", inferred, est_eqbm_t, eqbm_t_list_B[t], c, GP.C, selected_M, T_dict_B[t])
        next!(pb_B)
    end

    W_dfC_chunks = Vector{Any}(undef, T_length_C)
    dfC_inferred_conduct_data = Vector{Tuple{String, Int64, String}}(undef, T_length_C)

    pb_C = Progress(T_length_C, desc = "Processing C", dt = 0.5)

    Threads.@threads for t in 1:T_length_C
        dfC_t = @subset(dfC, :t .== T_dict_C[t])
        inferred = inferred_conduct_t_closest_f(dfC_t.σD_IV[1], dfC_t.σDSB_GMM[1], dfC_t.σDSC_GMM[1])
        dfC_inferred_conduct_data[t] = ("C", T_dict_C[t], inferred)
        est_eqbm_t = if inferred == "B"
            estimated_eqbm_t_f("B", eqbm_t_list_C[t], dfC_t.σDSB_GMM[1])
        elseif inferred == "C"
            estimated_eqbm_t_f("C", eqbm_t_list_C[t], dfC_t.σDSC_GMM[1])
        else
            eqbm_t_list_C[t]  # use true equilibrium, doesn't matter
        end
        W_dfC_chunks[t] = optimal_s_inferred_true_full_t("C", inferred, est_eqbm_t, eqbm_t_list_C[t], c, GP.C, selected_M, T_dict_C[t])
        next!(pb_C)
    end

    dfB_inferred_conduct = DataFrame(DGP = [x[1] for x in dfB_inferred_conduct_data],
                                     t = [x[2] for x in dfB_inferred_conduct_data],
                                     optimal_conduct = [x[3] for x in dfB_inferred_conduct_data])
    dfC_inferred_conduct = DataFrame(DGP = [x[1] for x in dfC_inferred_conduct_data],
                                     t = [x[2] for x in dfC_inferred_conduct_data],
                                     optimal_conduct = [x[3] for x in dfC_inferred_conduct_data])

    W_dfB = vcat(W_dfB_chunks...)
    W_dfC = vcat(W_dfC_chunks...)

    if Output
        save(OutputFile, "dfB_inferred_conduct", dfB_inferred_conduct, "dfC_inferred_conduct", dfC_inferred_conduct, "W_dfB", W_dfB, "W_dfC", W_dfC)
    end

    return W_dfB, W_dfC
end




function policy_table_inferred_true_w_numT_diff_samples_f(dfB, dfC, total_TB::Int64, total_TC::Int64, M::Int64)
    """
    Create the policy table, check the number of T
    Input:
        - W_dfB: DataFrame containing the optimal conduct for Bertrand
        - W_dfC: DataFrame containing the optimal conduct for Cournot
    Output:
        - B_policy_table: DataFrame containing the policy table for Bertrand
        - C_policy_table: DataFrame containing the policy table for Cournot
   """
    W_dfB0 = deepcopy(dfB)
    W_dfC0 = deepcopy(dfC)
    total_mktB = total_TB * M
    total_mktC = total_TC * M
    B_policy_table = DataFrame()
    C_policy_table = DataFrame()
    if isempty(W_dfB0)
        B_policy_table = DataFrame(
            T=total_TB, inferred_B_pct=0.0,
            inferred_C_pct=0.0,
            subsidy_pct=NaN,
            tax_pct=NaN,
            none_pct=1.0,
            NaN_pct=0.0,
            avg_subsidy=NaN,
            avg_tax=NaN,
            avg_W_pct=NaN, 
            true_inferrred_B_pct=missing,
            true_inferrred_C_pct=missing,
            true_subsidy_pct=missing,
            true_tax_pct=missing,
            true_none_pct=missing,
            true_NaN_pct=missing,
            avg_true_subsidy=missing,
            avg_true_tax=missing,
            avg_true_W_pct=missing
        )
    end
    if isempty(W_dfC0)
        C_policy_table = DataFrame(
            T=total_TC, inferred_B_pct=0.0,
            inferred_C_pct=0.0,
            subsidy_pct=NaN,
            tax_pct=NaN,
            none_pct=1.0,
            NaN_pct=0.0,
            avg_subsidy=NaN,
            avg_tax=NaN,
            avg_W_pct=NaN, 
            true_inferrred_B_pct=missing,
            true_inferred_C_pct=missing,
            true_subsidy_pct=missing,
            true_tax_pct=missing,
            true_none_pct=missing,
            true_NaN_pct=missing,
            avg_true_subsidy=missing,
            avg_true_tax=missing,
            avg_true_W_pct=missing
        )
    end
    if !isempty(W_dfB0)
        # check NaN, no convergence
        B_S_NaN = sum(isnan.(W_dfB0.S_optimal_m))
        B_S_true_NaN = sum(isnan.(W_dfB0.S_optimal_true_m))
        # drop NaN
        W_dfB = @subset(W_dfB0, .!isnan.(:S_optimal_m))
        W_true_dfB = @subset(W_dfB0, .!isnan.(:S_optimal_true_m))
        
        B_S_optimal_subsidy = sum(W_dfB.S_optimal_m .> 0)
        B_S_optimal_tax = sum(W_dfB.S_optimal_m .< 0)
        # B_S_optimal_none = total_mktB - B_S_optimal_subsidy - B_S_optimal_tax - B_S_NaN
        B_S_optimal_none = sum(W_dfB.inferred_conduct .== "N")

        B_inferred_conduct_B = sum(W_dfB.inferred_conduct .== "B")
        B_inferred_conduct_C = sum(W_dfB.inferred_conduct .== "C")

        B_S_optimal_true_subsidy = sum(W_true_dfB.S_optimal_true_m .> 0)
        B_S_optimal_true_tax = sum(W_true_dfB.S_optimal_true_m .< 0)
        # B_S_optimal_true_none = total_mktB - B_S_optimal_true_subsidy - B_S_optimal_true_tax - B_S_true_NaN
        B_S_optimal_true_none = sum(W_true_dfB.S_optimal_true_m .== 0) # should be 0

        B_inferred_conduct_true_B = total_mktB - B_S_optimal_true_none - B_S_true_NaN  
        B_inferred_conduct_true_C = 0.0

        B_policy_table = DataFrame(
            T=total_TB, inferred_B_pct=round(B_inferred_conduct_B / total_mktB, digits=3),
            inferred_C_pct=round(B_inferred_conduct_C / total_mktB, digits=3),
            subsidy_pct=round(B_S_optimal_subsidy / total_mktB, digits=3),
            tax_pct=round(B_S_optimal_tax / total_mktB, digits=3),
            none_pct=round(B_S_optimal_none / total_mktB, digits=3),
            NaN_pct=round(B_S_NaN / total_mktB, digits=3),
            avg_subsidy=round(mean(W_dfB.S_optimal_m[W_dfB.S_optimal_m.>0]), digits=3),
            avg_tax=round(mean(W_dfB.S_optimal_m[W_dfB.S_optimal_m.<0]), digits=3),
            avg_W_pct=round(mean(W_dfB.pctW_m[.!isnan.(W_dfB.S_optimal_m)]), digits=3), true_inferrred_B_pct=round(B_inferred_conduct_true_B / total_mktB, digits=3),
            true_inferrred_C_pct=round(B_inferred_conduct_true_C / total_mktB, digits=3),
            true_subsidy_pct=round(B_S_optimal_true_subsidy / total_mktB, digits=3),
            true_tax_pct=round(B_S_optimal_true_tax / total_mktB, digits=3),
            true_none_pct=round(B_S_optimal_true_none / total_mktB, digits=3),
            true_NaN_pct=round(B_S_true_NaN / total_mktB, digits=3),
            avg_true_subsidy=round(mean(W_true_dfB.S_optimal_true_m[W_true_dfB.S_optimal_true_m.>0]), digits=3),
            avg_true_tax=round(mean(W_true_dfB.S_optimal_true_m[W_true_dfB.S_optimal_true_m.<0]), digits=3),
            avg_true_W_pct=round(mean(W_true_dfB.pctW_true_m[.!isnan.(W_true_dfB.S_optimal_true_m)]), digits=3)
        )
    end
    if !isempty(W_dfC0)
        # check NaN
        C_S_NaN = sum(isnan.(W_dfC0.S_optimal_m))
        C_S_true_NaN = sum(isnan.(W_dfC0.S_optimal_true_m))
        # drop NaN
        W_dfC = @subset(W_dfC0, .!isnan.(:S_optimal_m))
        W_true_dfC = @subset(W_dfC0, .!isnan.(:S_optimal_true_m))
        C_S_optimal_subsidy = sum(W_dfC.S_optimal_m .> 0)
        C_S_optimal_tax = sum(W_dfC.S_optimal_m .< 0)
        # C_S_optimal_none = total_mktC - C_S_optimal_subsidy - C_S_optimal_tax - C_S_NaN
        C_S_optimal_none = sum(W_dfC.inferred_conduct .== "N") 
        C_inferred_conduct_B = sum(W_dfC.inferred_conduct .== "B")
        C_inferred_conduct_C = sum(W_dfC.inferred_conduct .== "C")

        C_S_optimal_true_subsidy = sum(W_true_dfC.S_optimal_true_m .> 0)
        C_S_optimal_true_tax = sum(W_true_dfC.S_optimal_true_m .< 0)
        # C_S_optimal_true_none = total_mktC - C_S_optimal_true_subsidy - C_S_optimal_true_tax - C_S_true_NaN
        C_S_optimal_true_none = sum(W_true_dfC.S_optimal_true_m .== 0) # should be 0
        C_inferred_conduct_true_B = 0.0 # should be 0
        C_inferred_conduct_true_C = total_mktC - C_S_optimal_true_none - C_S_true_NaN

        C_policy_table = DataFrame(
            T=total_TC, inferred_B_pct=round(C_inferred_conduct_B / total_mktC, digits=3),
            inferred_C_pct=round(C_inferred_conduct_C / total_mktC, digits=3),
            subsidy_pct=round(C_S_optimal_subsidy / total_mktC, digits=3),
            tax_pct=round(C_S_optimal_tax / total_mktC, digits=3),
            none_pct=round(C_S_optimal_none / total_mktC, digits=3),
            NaN_pct=round(C_S_NaN / total_mktC, digits=3),
            avg_subsidy=round(mean(W_dfC.S_optimal_m[W_dfC.S_optimal_m.>0]), digits=3),
            avg_tax=round(mean(W_dfC.S_optimal_m[W_dfC.S_optimal_m.<0]), digits=3),
            avg_W_pct=round(mean(W_dfC.pctW_m[.!isnan.(W_dfC.S_optimal_m)]), digits=3), true_inferrred_B_pct=round(C_inferred_conduct_true_B / total_mktC, digits=3),
            true_inferred_C_pct=round(C_inferred_conduct_true_C / total_mktC, digits=3),
            true_subsidy_pct=round(C_S_optimal_true_subsidy / total_mktC, digits=3),
            true_tax_pct=round(C_S_optimal_true_tax / total_mktC, digits=3),
            true_none_pct=round(C_S_optimal_true_none / total_mktC, digits=3),
            true_NaN_pct=round(C_S_true_NaN / total_mktC, digits=3),
            avg_true_subsidy=round(mean(W_true_dfC.S_optimal_true_m[W_true_dfC.S_optimal_true_m.>0]), digits=3),
            avg_true_tax=round(mean(W_true_dfC.S_optimal_true_m[W_true_dfC.S_optimal_true_m.<0]), digits=3),
            avg_true_W_pct=round(mean(W_true_dfC.pctW_true_m[.!isnan.(W_true_dfC.S_optimal_true_m)]), digits=3)
        )
    end
    return B_policy_table, C_policy_table

end



# at this point simply keep a same sample that have non-NaN for both inferred and true conduct
function policy_table_inferred_true_w_numT_f(dfB, dfC, total_TB::Int64, total_TC::Int64, M::Int64)
    """
    Create the policy table, check the number of T
    Input:
        - W_dfB: DataFrame containing the optimal conduct for Bertrand
        - W_dfC: DataFrame containing the optimal conduct for Cournot
    Output:
        - B_policy_table: DataFrame containing the policy table for Bertrand
        - C_policy_table: DataFrame containing the policy table for Cournot
   """
    W_dfB0 = deepcopy(dfB)
    W_dfC0 = deepcopy(dfC)
    total_mktB = total_TB * M
    total_mktC = total_TC * M
    B_policy_table = DataFrame()
    C_policy_table = DataFrame()
    if isempty(W_dfB0)
        B_policy_table = DataFrame(
            T=total_TB, 
            total_nonNaN_obs = missing,
            NaN_pct = missing,
            inferred_B_pct= missing,
            inferred_C_pct= missing,
            subsidy_pct= missing,
            tax_pct= missing,
            none_pct= missing,
            avg_subsidy=missing,
            avg_tax=missing,
            avg_W_pct=missing,
            min_W_pct = missing,
            pct_10_W_pct = missing,
            pct_25_W_pct = missing,
            pct_75_W_pct = missing,
            pct_90_W_pct = missing,
            max_W_pct = missing,

            true_subsidy_pct=missing,
            true_tax_pct=missing,
            avg_true_W_pct=missing
        )
    end
    if isempty(W_dfC0)
        C_policy_table = DataFrame(
            T=total_TC, 
            total_nonNaN_obs = missing,
            NaN_pct = missing,
            inferred_B_pct= missing,
            inferred_C_pct= missing,
            subsidy_pct= missing,
            tax_pct= missing,
            none_pct= missing,
            avg_subsidy=missing,
            avg_tax=missing,
            avg_W_pct=missing,
            min_W_pct = missing,
            pct_10_W_pct = missing,
            pct_25_W_pct = missing,
            pct_75_W_pct = missing,
            pct_90_W_pct = missing,
            max_W_pct = missing,

            true_subsidy_pct=missing,
            true_tax_pct=missing,
            avg_true_W_pct=missing
        )
    end
    if !isempty(W_dfB0)

        # drop NaN based on both inferred and true
        W_dfB = @subset(W_dfB0, .!isnan.(:S_optimal_m))
        W_dfB = @subset(W_dfB, .!isnan.(:S_optimal_true_m))
        W_dfB = @subset(W_dfB, .!isnan.(:Woptimal_m)) 
        W_dfB = @subset(W_dfB, .!isnan.(:Woptimal_true_m)) 

        # W_dfB_subset = @subset(W_dfB, .!isnan.(:S_optimal_true_m))
        
        # check the number of dropped observations
        B_S_union_NaN = total_mktB - size(W_dfB, 1)

        B_S_optimal_subsidy = sum(W_dfB.S_optimal_m .> 0)
        B_S_optimal_tax = sum(W_dfB.S_optimal_m .< 0)
        # B_S_optimal_none = total_mktB - B_S_optimal_subsidy - B_S_optimal_tax - B_S_NaN
        B_S_optimal_none = sum(W_dfB.inferred_conduct .== "N")

        B_inferred_conduct_B = sum(W_dfB.inferred_conduct .== "B")
        B_inferred_conduct_C = sum(W_dfB.inferred_conduct .== "C")

        B_S_optimal_true_subsidy = sum(W_dfB.S_optimal_true_m .> 0)
        B_S_optimal_true_tax = sum(W_dfB.S_optimal_true_m .< 0)
        # # B_S_optimal_true_none = total_mktB - B_S_optimal_true_subsidy - B_S_optimal_true_tax - B_S_true_NaN
        # B_S_optimal_true_none = sum(W_true_dfB.S_optimal_true_m .== 0) # should be 0， in other words, the sum of subsidy and tax percentage should be 1

        # B_inferred_conduct_true_B = total_mktB - B_S_optimal_true_none - B_S_true_NaN  
        # B_inferred_conduct_true_C = 0.0

        B_policy_table = DataFrame(
            T=total_TB, 
            total_nonNaN_obs = size(W_dfB, 1),
            NaN_pct=round(B_S_union_NaN / total_mktB, digits=3),
            inferred_B_pct=round(B_inferred_conduct_B / total_mktB, digits=3),
            inferred_C_pct=round(B_inferred_conduct_C / total_mktB, digits=3),
            subsidy_pct=round(B_S_optimal_subsidy / total_mktB, digits=3),
            tax_pct=round(B_S_optimal_tax / total_mktB, digits=3),
            none_pct=round(B_S_optimal_none / total_mktB, digits=3),
            avg_subsidy=round(mean(W_dfB.S_optimal_m[W_dfB.S_optimal_m.>0]), digits=3),
            avg_tax=round(mean(W_dfB.S_optimal_m[W_dfB.S_optimal_m.<0]), digits=3),
            avg_W_pct=round(mean(W_dfB.pctW_m), digits=3),
            min_W_pct=round(minimum(W_dfB.pctW_m), digits=3),
            #pct_5_W_pct=round(quantile(W_dfB.pctW_m, 0.05), digits=3),
            pct_10_W_pct=round(quantile(W_dfB.pctW_m, 0.10), digits=3),
            pct_25_W_pct=round(quantile(W_dfB.pctW_m, 0.25), digits=3),
            pct_75_W_pct=round(quantile(W_dfB.pctW_m, 0.75), digits=3), 
            pct_90_W_pct=round(quantile(W_dfB.pctW_m, 0.90), digits=3),
            # pct_95_W_pct=round(quantile(W_dfB.pctW_m, 0.95), digits=3),
            max_W_pct=round(maximum(W_dfB.pctW_m), digits=3),   
            #true_inferrred_B_pct=round(B_inferred_conduct_true_B / total_mktB, digits=3),
            #true_inferrred_C_pct=round(B_inferred_conduct_true_C / total_mktB, digits=3),
            true_subsidy_pct=round(B_S_optimal_true_subsidy / total_mktB, digits=3),
            true_tax_pct=round(B_S_optimal_true_tax / total_mktB, digits=3),
            # true_none_pct=round(B_S_optimal_true_none / total_mktB, digits=3),
            # true_NaN_pct=round(B_S_true_NaN / total_mktB, digits=3),
            # avg_true_subsidy=round(mean(W_true_dfB.S_optimal_true_m), digits=3),
            # avg_true_tax=round(mean(W_true_dfB.S_optimal_true_m), digits=3),
            avg_true_W_pct=round(mean(W_dfB.pctW_true_m), digits=3)
        )
    end
    if !isempty(W_dfC0)
        
        #  # drop NaN based on both inferred and true
        W_dfC = @subset(W_dfC0, .!isnan.(:S_optimal_m))
        W_dfC = @subset(W_dfC, .!isnan.(:S_optimal_true_m)) 
        W_dfC = @subset(W_dfC, .!isnan.(:Woptimal_m)) 
        W_dfC = @subset(W_dfC, .!isnan.(:Woptimal_true_m)) 
         # check the number of dropped observations
        C_S_union_NaN = total_mktC - size(W_dfC, 1)
        
        C_S_optimal_subsidy = sum(W_dfC.S_optimal_m .> 0)
        C_S_optimal_tax = sum(W_dfC.S_optimal_m .< 0)
        # C_S_optimal_none = total_mktC - C_S_optimal_subsidy - C_S_optimal_tax - C_S_NaN
        C_S_optimal_none = sum(W_dfC.inferred_conduct .== "N") 
        C_inferred_conduct_B = sum(W_dfC.inferred_conduct .== "B")
        C_inferred_conduct_C = sum(W_dfC.inferred_conduct .== "C")

        C_S_optimal_true_subsidy = sum(W_dfC.S_optimal_true_m .> 0)
        C_S_optimal_true_tax = sum(W_dfC.S_optimal_true_m .< 0)
        # C_S_optimal_true_none = total_mktC - C_S_optimal_true_subsidy - C_S_optimal_true_tax - C_S_true_NaN
        # C_S_optimal_true_none = sum(W_true_dfC.S_optimal_true_m .== 0) # should be 0
        # C_inferred_conduct_true_B = 0.0 # should be 0, in other words, the sum of subsidy and tax percentage should be 1
        # C_inferred_conduct_true_C = total_mktC - C_S_optimal_true_none - C_S_true_NaN
        
        # # calculate the distribution of W_dfC.pctW_m
        # sum(isnan.(W_dfC.pctW_m))
        # W_dfC_nan = @subset(W_dfC, isnan.(:pctW_m))

        C_policy_table = DataFrame(
            T=total_TC, 
            total_nonNaN_obs = size(W_dfC, 1),
            NaN_pct=round(C_S_union_NaN / total_mktC, digits=3),
            inferred_B_pct=round(C_inferred_conduct_B / total_mktC, digits=3),
            inferred_C_pct=round(C_inferred_conduct_C / total_mktC, digits=3),
            subsidy_pct=round(C_S_optimal_subsidy / total_mktC, digits=3),
            tax_pct=round(C_S_optimal_tax / total_mktC, digits=3),
            none_pct=round(C_S_optimal_none / total_mktC, digits=3),
            avg_subsidy=round(mean(W_dfC.S_optimal_m[W_dfC.S_optimal_m.>0]), digits=3),
            avg_tax=round(mean(W_dfC.S_optimal_m[W_dfC.S_optimal_m.<0]), digits=3),
            avg_W_pct=round(mean(W_dfC.pctW_m), digits=3), 
            min_W_pct=round(minimum(W_dfC.pctW_m), digits=3),
            #pct_5_W_pct=round(quantile(W_dfC.pctW_m, 0.05), digits=3),
            pct_10_W_pct=round(quantile(W_dfC.pctW_m, 0.10), digits=3),
            pct_25_W_pct=round(quantile(W_dfC.pctW_m, 0.25), digits=3),
            pct_75_W_pct=round(quantile(W_dfC.pctW_m, 0.75), digits=3),
            pct_90_W_pct=round(quantile(W_dfC.pctW_m, 0.90), digits=3),
            # pct_95_W_pct=round(quantile(W_dfC.pctW_m, 0.95), digits=3),
            max_W_pct=round(maximum(W_dfC.pctW_m), digits=3),

            # true_inferrred_B_pct=round(C_inferred_conduct_true_B / total_mktC, digits=3),
            # true_inferred_C_pct=round(C_inferred_conduct_true_C / total_mktC, digits=3),
            true_subsidy_pct=round(C_S_optimal_true_subsidy / total_mktC, digits=3),
            true_tax_pct=round(C_S_optimal_true_tax / total_mktC, digits=3),
            # true_none_pct=round(C_S_optimal_true_none / total_mktC, digits=3),
            # true_NaN_pct=round(C_S_true_NaN / total_mktC, digits=3),
            # avg_true_subsidy=round(mean(W_true_dfC.S_optimal_true_m), digits=3),
            # avg_true_tax=round(mean(W_true_dfC.S_optimal_true_m), digits=3),
            avg_true_W_pct=round(mean(W_dfC.pctW_true_m), digits=3)
        )
    end
    return B_policy_table, C_policy_table
end




function policy_table_inferred_true_w_numT_reduced_f(dfB, dfC, total_TB::Int64, total_TC::Int64, M::Int64)
    """
    Create the policy table, check the number of T
    Input:
        - W_dfB: DataFrame containing the optimal conduct for Bertrand
        - W_dfC: DataFrame containing the optimal conduct for Cournot
    Output:
        - B_policy_table: DataFrame containing the policy table for Bertrand
        - C_policy_table: DataFrame containing the policy table for Cournot
   """
    W_dfB0 = deepcopy(dfB)
    W_dfC0 = deepcopy(dfC)
    total_mktB = total_TB * M
    total_mktC = total_TC * M
    B_policy_table = DataFrame()
    C_policy_table = DataFrame()
    if isempty(W_dfB0)
        B_policy_table = DataFrame(
            #T=total_TB, 
            total_nonNaN_obs = missing,
            #NaN_pct = missing,
            #inferred_B_pct= missing,
            #inferred_C_pct= missing,
            subsidy_pct= missing,
            tax_pct= missing,
            #none_pct= missing,
            #avg_subsidy=missing,
            #avg_tax=missing,
            avg_W_pct=missing,
            min_W_pct = missing,
            #pct_10_W_pct = missing,
            pct_25_W_pct = missing,
            pct_75_W_pct = missing,
            #pct_90_W_pct = missing,
            max_W_pct = missing,

            true_subsidy_pct=missing,
            true_tax_pct=missing,
            avg_true_W_pct=missing
        )
    end
    if isempty(W_dfC0)
        C_policy_table = DataFrame(
            #T=total_TC, 
            total_nonNaN_obs = missing,
            #NaN_pct = missing,
            #inferred_B_pct= missing,
            #inferred_C_pct= missing,
            subsidy_pct= missing,
            tax_pct= missing,
            #none_pct= missing,
            #avg_subsidy=missing,
            #avg_tax=missing,
            avg_W_pct=missing,
            min_W_pct = missing,
            #pct_10_W_pct = missing,
            pct_25_W_pct = missing,
            pct_75_W_pct = missing,
            #pct_90_W_pct = missing,
            max_W_pct = missing,

            true_subsidy_pct=missing,
            true_tax_pct=missing,
            avg_true_W_pct=missing
        )
    end
    if !isempty(W_dfB0)

        # drop NaN based on both inferred and true
        W_dfB = @subset(W_dfB0, .!isnan.(:S_optimal_m))
        W_dfB = @subset(W_dfB, .!isnan.(:S_optimal_true_m))
        W_dfB = @subset(W_dfB, .!isnan.(:Woptimal_m)) 
        W_dfB = @subset(W_dfB, .!isnan.(:Woptimal_true_m)) 

        # W_dfB_subset = @subset(W_dfB, .!isnan.(:S_optimal_true_m))
        
        # check the number of dropped observations
        B_S_union_NaN = total_mktB - size(W_dfB, 1)

        B_S_optimal_subsidy = sum(W_dfB.S_optimal_m .> 0)
        B_S_optimal_tax = sum(W_dfB.S_optimal_m .< 0)
        # B_S_optimal_none = total_mktB - B_S_optimal_subsidy - B_S_optimal_tax - B_S_NaN
        B_S_optimal_none = sum(W_dfB.inferred_conduct .== "N")

        B_inferred_conduct_B = sum(W_dfB.inferred_conduct .== "B")
        B_inferred_conduct_C = sum(W_dfB.inferred_conduct .== "C")

        B_S_optimal_true_subsidy = sum(W_dfB.S_optimal_true_m .> 0)
        B_S_optimal_true_tax = sum(W_dfB.S_optimal_true_m .< 0)
        # # B_S_optimal_true_none = total_mktB - B_S_optimal_true_subsidy - B_S_optimal_true_tax - B_S_true_NaN
        # B_S_optimal_true_none = sum(W_true_dfB.S_optimal_true_m .== 0) # should be 0， in other words, the sum of subsidy and tax percentage should be 1

        # B_inferred_conduct_true_B = total_mktB - B_S_optimal_true_none - B_S_true_NaN  
        # B_inferred_conduct_true_C = 0.0

        B_policy_table = DataFrame(
            #T=total_TB, 
            total_nonNaN_obs = size(W_dfB, 1),
            #NaN_pct=round(B_S_union_NaN / total_mktB, digits=3),
            #inferred_B_pct=round(B_inferred_conduct_B / total_mktB, digits=3),
            #inferred_C_pct=round(B_inferred_conduct_C / total_mktB, digits=3),
            subsidy_pct=round(B_S_optimal_subsidy / total_mktB, digits=3),
            tax_pct=round(B_S_optimal_tax / total_mktB, digits=3),
            #none_pct=round(B_S_optimal_none / total_mktB, digits=3),
            #avg_subsidy=round(mean(W_dfB.S_optimal_m[W_dfB.S_optimal_m.>0]), digits=3),
            #avg_tax=round(mean(W_dfB.S_optimal_m[W_dfB.S_optimal_m.<0]), digits=3),
            avg_W_pct=round(mean(W_dfB.pctW_m), digits=3),
            min_W_pct=round(minimum(W_dfB.pctW_m), digits=3),
            #pct_5_W_pct=round(quantile(W_dfB.pctW_m, 0.05), digits=3),
            #pct_10_W_pct=round(quantile(W_dfB.pctW_m, 0.10), digits=3),
            pct_25_W_pct=round(quantile(W_dfB.pctW_m, 0.25), digits=3),
            pct_75_W_pct=round(quantile(W_dfB.pctW_m, 0.75), digits=3), 
            #pct_90_W_pct=round(quantile(W_dfB.pctW_m, 0.90), digits=3),
            # pct_95_W_pct=round(quantile(W_dfB.pctW_m, 0.95), digits=3),
            max_W_pct=round(maximum(W_dfB.pctW_m), digits=3),   
            #true_inferrred_B_pct=round(B_inferred_conduct_true_B / total_mktB, digits=3),
            #true_inferrred_C_pct=round(B_inferred_conduct_true_C / total_mktB, digits=3),
            true_subsidy_pct=round(B_S_optimal_true_subsidy / total_mktB, digits=3),
            true_tax_pct=round(B_S_optimal_true_tax / total_mktB, digits=3),
            # true_none_pct=round(B_S_optimal_true_none / total_mktB, digits=3),
            # true_NaN_pct=round(B_S_true_NaN / total_mktB, digits=3),
            # avg_true_subsidy=round(mean(W_true_dfB.S_optimal_true_m), digits=3),
            # avg_true_tax=round(mean(W_true_dfB.S_optimal_true_m), digits=3),
            avg_true_W_pct=round(mean(W_dfB.pctW_true_m), digits=3)
        )
    end
    if !isempty(W_dfC0)
        
        #  # drop NaN based on both inferred and true
        W_dfC = @subset(W_dfC0, .!isnan.(:S_optimal_m))
        W_dfC = @subset(W_dfC, .!isnan.(:S_optimal_true_m)) 
        W_dfC = @subset(W_dfC, .!isnan.(:Woptimal_m)) 
        W_dfC = @subset(W_dfC, .!isnan.(:Woptimal_true_m)) 
         # check the number of dropped observations
        C_S_union_NaN = total_mktC - size(W_dfC, 1)
        
        C_S_optimal_subsidy = sum(W_dfC.S_optimal_m .> 0)
        C_S_optimal_tax = sum(W_dfC.S_optimal_m .< 0)
        # C_S_optimal_none = total_mktC - C_S_optimal_subsidy - C_S_optimal_tax - C_S_NaN
        C_S_optimal_none = sum(W_dfC.inferred_conduct .== "N") 
        C_inferred_conduct_B = sum(W_dfC.inferred_conduct .== "B")
        C_inferred_conduct_C = sum(W_dfC.inferred_conduct .== "C")

        C_S_optimal_true_subsidy = sum(W_dfC.S_optimal_true_m .> 0)
        C_S_optimal_true_tax = sum(W_dfC.S_optimal_true_m .< 0)
        # C_S_optimal_true_none = total_mktC - C_S_optimal_true_subsidy - C_S_optimal_true_tax - C_S_true_NaN
        # C_S_optimal_true_none = sum(W_true_dfC.S_optimal_true_m .== 0) # should be 0
        # C_inferred_conduct_true_B = 0.0 # should be 0, in other words, the sum of subsidy and tax percentage should be 1
        # C_inferred_conduct_true_C = total_mktC - C_S_optimal_true_none - C_S_true_NaN
        
        # # calculate the distribution of W_dfC.pctW_m
        # sum(isnan.(W_dfC.pctW_m))
        # W_dfC_nan = @subset(W_dfC, isnan.(:pctW_m))

        C_policy_table = DataFrame(
            #T=total_TC, 
            total_nonNaN_obs = size(W_dfC, 1),
            #NaN_pct=round(C_S_union_NaN / total_mktC, digits=3),
            #inferred_B_pct=round(C_inferred_conduct_B / total_mktC, digits=3),
            #inferred_C_pct=round(C_inferred_conduct_C / total_mktC, digits=3),
            subsidy_pct=round(C_S_optimal_subsidy / total_mktC, digits=3),
            tax_pct=round(C_S_optimal_tax / total_mktC, digits=3),
            #none_pct=round(C_S_optimal_none / total_mktC, digits=3),
            #avg_subsidy=round(mean(W_dfC.S_optimal_m[W_dfC.S_optimal_m.>0]), digits=3),
            #avg_tax=round(mean(W_dfC.S_optimal_m[W_dfC.S_optimal_m.<0]), digits=3),
            avg_W_pct=round(mean(W_dfC.pctW_m), digits=3), 
            min_W_pct=round(minimum(W_dfC.pctW_m), digits=3),
            #pct_5_W_pct=round(quantile(W_dfC.pctW_m, 0.05), digits=3),
            #pct_10_W_pct=round(quantile(W_dfC.pctW_m, 0.10), digits=3),
            pct_25_W_pct=round(quantile(W_dfC.pctW_m, 0.25), digits=3),
            pct_75_W_pct=round(quantile(W_dfC.pctW_m, 0.75), digits=3),
            #pct_90_W_pct=round(quantile(W_dfC.pctW_m, 0.90), digits=3),
            # pct_95_W_pct=round(quantile(W_dfC.pctW_m, 0.95), digits=3),
            max_W_pct=round(maximum(W_dfC.pctW_m), digits=3),

            # true_inferrred_B_pct=round(C_inferred_conduct_true_B / total_mktC, digits=3),
            # true_inferred_C_pct=round(C_inferred_conduct_true_C / total_mktC, digits=3),
            true_subsidy_pct=round(C_S_optimal_true_subsidy / total_mktC, digits=3),
            true_tax_pct=round(C_S_optimal_true_tax / total_mktC, digits=3),
            # true_none_pct=round(C_S_optimal_true_none / total_mktC, digits=3),
            # true_NaN_pct=round(C_S_true_NaN / total_mktC, digits=3),
            # avg_true_subsidy=round(mean(W_true_dfC.S_optimal_true_m), digits=3),
            # avg_true_tax=round(mean(W_true_dfC.S_optimal_true_m), digits=3),
            avg_true_W_pct=round(mean(W_dfC.pctW_true_m), digits=3)
        )
    end
    return B_policy_table, C_policy_table
end









# at this point simply keep a same sample that have non-NaN for both inferred and true conduct
function policy_table_inferred_true_w_numT_pct_f(dfB, dfC, total_TB::Int64, total_TC::Int64, M::Int64)
    """
    Create the policy table, check the number of T, output for figures
    Input:
        - W_dfB: DataFrame containing the optimal conduct for Bertrand
        - W_dfC: DataFrame containing the optimal conduct for Cournot
    Output:
        - B_policy_table: DataFrame containing the policy table for Bertrand
        - C_policy_table: DataFrame containing the policy table for Cournot
   """
    W_dfB0 = deepcopy(dfB)
    W_dfC0 = deepcopy(dfC)
    total_mktB = total_TB * M
    total_mktC = total_TC * M
    B_policy_table = DataFrame()
    C_policy_table = DataFrame()
    if isempty(W_dfB0)
        B_policy_table = DataFrame(
            T=total_TB, 
            total_nonNaN_obs = missing,
            NaN_pct = missing,
            inferred_B_pct= missing,
            inferred_C_pct= missing,
            subsidy_pct= missing,
            tax_pct= missing,
            none_pct= missing,
            avg_subsidy=missing,
            avg_tax=missing,
            avg_inferred_W_pct_true = missing,

            true_subsidy_pct=missing,
            true_tax_pct=missing
        )
    end
    if isempty(W_dfC0)
        C_policy_table = DataFrame(
            T=total_TC, 
            total_nonNaN_obs = missing,
            NaN_pct = missing,
            inferred_B_pct= missing,
            inferred_C_pct= missing,
            subsidy_pct= missing,
            tax_pct= missing,
            none_pct= missing,
            avg_subsidy=missing,
            avg_tax=missing,
            avg_inferred_W_pct_true = missing,

            true_subsidy_pct=missing,
            true_tax_pct=missing
        )
    end
    if !isempty(W_dfB0)

        # drop NaN based on both inferred and true
        W_dfB = @subset(W_dfB0, .!isnan.(:S_optimal_m))
        W_dfB = @subset(W_dfB, .!isnan.(:S_optimal_true_m))
        W_dfB = @subset(W_dfB, .!isnan.(:Woptimal_m)) 
        W_dfB = @subset(W_dfB, .!isnan.(:Woptimal_true_m)) 

        # W_dfB_subset = @subset(W_dfB, .!isnan.(:S_optimal_true_m))
        
        # check the number of dropped observations
        B_S_union_NaN = total_mktB - size(W_dfB, 1)

        B_S_optimal_subsidy = sum(W_dfB.S_optimal_m .> 0)
        B_S_optimal_tax = sum(W_dfB.S_optimal_m .< 0)
        # B_S_optimal_none = total_mktB - B_S_optimal_subsidy - B_S_optimal_tax - B_S_NaN
        B_S_optimal_none = sum(W_dfB.inferred_conduct .== "N")

        B_inferred_conduct_B = sum(W_dfB.inferred_conduct .== "B")
        B_inferred_conduct_C = sum(W_dfB.inferred_conduct .== "C")

        B_S_optimal_true_subsidy = sum(W_dfB.S_optimal_true_m .> 0)
        B_S_optimal_true_tax = sum(W_dfB.S_optimal_true_m .< 0)
        # # B_S_optimal_true_none = total_mktB - B_S_optimal_true_subsidy - B_S_optimal_true_tax - B_S_true_NaN
        # B_S_optimal_true_none = sum(W_true_dfB.S_optimal_true_m .== 0) # should be 0， in other words, the sum of subsidy and tax percentage should be 1

        # B_inferred_conduct_true_B = total_mktB - B_S_optimal_true_none - B_S_true_NaN  
        # B_inferred_conduct_true_C = 0.0
        # W_dfB_clean = dropmissing(W_dfB, [:pctW_m, :pctW_true_m])
        W_dfB_clean = dropmissing(W_dfB, [:Woptimal_m, :Woptimal_true_m])

        W_dfB_clean_true = @subset(W_dfB0, .!isnan.(:S_optimal_true_m)) 
        W_dfB_clean_true = @subset(W_dfB_clean_true, .!isnan.(:Woptimal_true_m))
        W_dfB_clean_true =  dropmissing(W_dfB_clean_true, [:Woptimal_true_m])
        if size(W_dfB_clean, 1) != 0
           W_dfB_clean[!, :W_pct_ratio] = W_dfB_clean.deltaW_m ./ W_dfB_clean.deltaW_true_m
           # select rows with positive deltaW_true_m
           W_dfB_clean = @subset(W_dfB_clean,  :deltaW_true_m .> 0)
           W_dfB_clean = @subset(W_dfB_clean, .!isnan.(:W_pct_ratio))
           W_dfB_clean = @subset(W_dfB_clean, .!isinf.(:W_pct_ratio))

           W_dfB_clean_true = @subset(W_dfB_clean_true,  :deltaW_true_m .> 0)
        end
        if size(W_dfB_clean) != 0
            B_policy_table = DataFrame(
            T=total_TB, 
            total_nonNaN_obs = size(W_dfB, 1),
            NaN_pct=round(B_S_union_NaN / total_mktB, digits=3),
            inferred_B_pct=round(B_inferred_conduct_B / total_mktB, digits=3),
            inferred_C_pct=round(B_inferred_conduct_C / total_mktB, digits=3),
            subsidy_pct=round(B_S_optimal_subsidy / total_mktB, digits=3),
            tax_pct=round(B_S_optimal_tax / total_mktB, digits=3),
            none_pct=round(B_S_optimal_none / total_mktB, digits=3),
            avg_subsidy=round(mean(W_dfB.S_optimal_m[W_dfB.S_optimal_m.>0]), digits=3),
            avg_tax=round(mean(W_dfB.S_optimal_m[W_dfB.S_optimal_m.<0]), digits=3),

            # avg_inferred_W_pct_true = nrow(W_dfB_clean) == 0 ? missing : round(mean(W_dfB_clean.pctW_m[!W_dfB_clean.pctW_true_m .== 0] ./ W_dfB_clean.pctW_true_m[!W_dfB_clean.pctW_true_m .== 0]), digits=3),
            avg_true_W_pct = round(mean(W_dfB_clean_true.pctW_true_m), digits=3),
            avg_inferred_W_pct_true = round(mean(W_dfB_clean.W_pct_ratio), digits=3),
            median_inferred_W_pct_true = round(median(W_dfB_clean.W_pct_ratio), digits=3),

            true_subsidy_pct=round(B_S_optimal_true_subsidy / total_mktB, digits=3),
            true_tax_pct=round(B_S_optimal_true_tax / total_mktB, digits=3)   
        )

        else
            B_policy_table = DataFrame(
            T=total_TB, 
            total_nonNaN_obs = size(W_dfB, 1),
            NaN_pct=round(B_S_union_NaN / total_mktB, digits=3),
            inferred_B_pct=round(B_inferred_conduct_B / total_mktB, digits=3),
            inferred_C_pct=round(B_inferred_conduct_C / total_mktB, digits=3),
            subsidy_pct=round(B_S_optimal_subsidy / total_mktB, digits=3),
            tax_pct=round(B_S_optimal_tax / total_mktB, digits=3),
            none_pct=round(B_S_optimal_none / total_mktB, digits=3),
            avg_subsidy=round(mean(W_dfB.S_optimal_m[W_dfB.S_optimal_m.>0]), digits=3),
            avg_tax=round(mean(W_dfB.S_optimal_m[W_dfB.S_optimal_m.<0]), digits=3),

            # avg_inferred_W_pct_true = nrow(W_dfB_clean) == 0 ? missing : round(mean(W_dfB_clean.pctW_m[!W_dfB_clean.pctW_true_m .== 0] ./ W_dfB_clean.pctW_true_m[!W_dfB_clean.pctW_true_m .== 0]), digits=3),
            avg_true_W_pct = missing,
            avg_inferred_W_pct_true = missing,
            median_inferred_W_pct_true = missing,

            true_subsidy_pct=round(B_S_optimal_true_subsidy / total_mktB, digits=3),
            true_tax_pct=round(B_S_optimal_true_tax / total_mktB, digits=3)   
        )
        end


    end
    if !isempty(W_dfC0)
        
        #  # drop NaN based on both inferred and true
        W_dfC = @subset(W_dfC0, .!isnan.(:S_optimal_m))
        W_dfC = @subset(W_dfC, .!isnan.(:S_optimal_true_m)) 
        W_dfC = @subset(W_dfC, .!isnan.(:Woptimal_m)) 
        W_dfC = @subset(W_dfC, .!isnan.(:Woptimal_true_m)) 
         # check the number of dropped observations
        C_S_union_NaN = total_mktC - size(W_dfC, 1)
        
        C_S_optimal_subsidy = sum(W_dfC.S_optimal_m .> 0)
        C_S_optimal_tax = sum(W_dfC.S_optimal_m .< 0)
        # C_S_optimal_none = total_mktC - C_S_optimal_subsidy - C_S_optimal_tax - C_S_NaN
        C_S_optimal_none = sum(W_dfC.inferred_conduct .== "N") 
        C_inferred_conduct_B = sum(W_dfC.inferred_conduct .== "B")
        C_inferred_conduct_C = sum(W_dfC.inferred_conduct .== "C")

        C_S_optimal_true_subsidy = sum(W_dfC.S_optimal_true_m .> 0)
        C_S_optimal_true_tax = sum(W_dfC.S_optimal_true_m .< 0)
        # C_S_optimal_true_none = total_mktC - C_S_optimal_true_subsidy - C_S_optimal_true_tax - C_S_true_NaN
        # C_S_optimal_true_none = sum(W_true_dfC.S_optimal_true_m .== 0) # should be 0
        # C_inferred_conduct_true_B = 0.0 # should be 0, in other words, the sum of subsidy and tax percentage should be 1
        # C_inferred_conduct_true_C = total_mktC - C_S_optimal_true_none - C_S_true_NaN
        # W_dfC_clean = dropmissing(W_dfC, [:pctW_m, :pctW_true_m])
        W_dfC_clean = dropmissing(W_dfC, [:Woptimal_m, :Woptimal_true_m])
        
        W_dfC_clean_true = @subset(W_dfC0, .!isnan.(:S_optimal_true_m)) 
        W_dfC_clean_true = @subset(W_dfC_clean_true, .!isnan.(:Woptimal_true_m))
        W_dfC_clean_true =  dropmissing(W_dfC_clean_true, [:Woptimal_true_m])
        if nrow(W_dfC_clean) > 0
           #W_dfC_clean[!, :W_pct_ratio] = W_dfC_clean.pctW_m ./ W_dfC_clean.pctW_true_m
            W_dfC_clean[!, :W_pct_ratio] = W_dfC_clean.deltaW_m ./ W_dfC_clean.deltaW_true_m
            # select rows with positive deltaW_true_m
            W_dfC_clean = @subset(W_dfC_clean,  :deltaW_true_m .> 0)
            # drop NaN
            W_dfC_clean = @subset(W_dfC_clean, .!isnan.(:W_pct_ratio))
            # skip inf and -inf
            W_dfC_clean = @subset(W_dfC_clean, .!isinf.(:W_pct_ratio))


            W_dfB_clean_true = @subset(W_dfB_clean_true,  :deltaW_true_m .> 0)
        end
        if nrow(W_dfC_clean) > 0
        C_policy_table = DataFrame(
            T=total_TC, 
            total_nonNaN_obs = size(W_dfC, 1),
            NaN_pct=round(C_S_union_NaN / total_mktC, digits=3),
            inferred_B_pct=round(C_inferred_conduct_B / total_mktC, digits=3),
            inferred_C_pct=round(C_inferred_conduct_C / total_mktC, digits=3),
            subsidy_pct=round(C_S_optimal_subsidy / total_mktC, digits=3),
            tax_pct=round(C_S_optimal_tax / total_mktC, digits=3),
            none_pct=round(C_S_optimal_none / total_mktC, digits=3),
            avg_subsidy=round(mean(W_dfC.S_optimal_m[W_dfC.S_optimal_m.>0]), digits=3),
            avg_tax=round(mean(W_dfC.S_optimal_m[W_dfC.S_optimal_m.<0]), digits=3),
            # need to skip missing values

            avg_true_W_pct = round(mean(W_dfC_clean_true.pctW_true_m), digits=3),
            avg_inferred_W_pct_true = round(mean(W_dfC_clean.W_pct_ratio), digits=3),
            median_inferred_W_pct_true = round(median(W_dfC_clean.W_pct_ratio), digits=3),


            true_subsidy_pct=round(C_S_optimal_true_subsidy / total_mktC, digits=3),
            true_tax_pct=round(C_S_optimal_true_tax / total_mktC, digits=3)
            
        )
        else
            C_policy_table = DataFrame(
            T=total_TC, 
            total_nonNaN_obs = size(W_dfC, 1),
            NaN_pct=round(C_S_union_NaN / total_mktC, digits=3),
            inferred_B_pct=round(C_inferred_conduct_B / total_mktC, digits=3),
            inferred_C_pct=round(C_inferred_conduct_C / total_mktC, digits=3),
            subsidy_pct=round(C_S_optimal_subsidy / total_mktC, digits=3),
            tax_pct=round(C_S_optimal_tax / total_mktC, digits=3),
            none_pct=round(C_S_optimal_none / total_mktC, digits=3),
            avg_subsidy=round(mean(W_dfC.S_optimal_m[W_dfC.S_optimal_m.>0]), digits=3),
            avg_tax=round(mean(W_dfC.S_optimal_m[W_dfC.S_optimal_m.<0]), digits=3),
            # need to skip missing values
            avg_true_W_pct = missing,
            avg_inferred_W_pct_true = missing,
            median_inferred_W_pct_true = missing,
            true_subsidy_pct=round(C_S_optimal_true_subsidy / total_mktC, digits=3),
            true_tax_pct=round(C_S_optimal_true_tax / total_mktC, digits=3)
            
        )
    end
end
    return B_policy_table, C_policy_table
end





function policy_table_summary_stats(dfB, dfC, total_TB::Int64, total_TC::Int64, M::Int64)
    """
    Check the summary statistics behind the policy table
    Input:
        - W_dfB: DataFrame containing the optimal conduct for Bertrand
        - W_dfC: DataFrame containing the optimal conduct for Cournot
    Output:
        - B_stats_table
        - C_stats_table
   """
    W_dfB0 = deepcopy(dfB)
    W_dfC0 = deepcopy(dfC)
    total_mktB = total_TB * M
    total_mktC = total_TC * M
    B_stats_table = DataFrame()
    C_stats_table = DataFrame()
    if isempty(W_dfB0)
        B_stats_table = DataFrame(
            T = total_TB, 
            total_nonNaN_obs = missing,
            avg_W_pct=missing,
            min_W_pct = missing,
            st_dev_W_pct = missing,
            total_positive_N_pct = missing,
            total_zero_N_pct = missing,
            total_negative_N_pct = missing
        )
    end
    if isempty(W_dfC0)
        C_stats_table = DataFrame(
            T = total_TC, 
            total_nonNaN_obs = missing,
            avg_W_pct=missing,
            min_W_pct = missing,
            st_dev_W_pct = missing,
            total_positive_N_pct = missing,
            total_zero_N_pct = missing,
            total_negative_N_pct = missing
        )
    end
    if !isempty(W_dfB0)
        # drop NaN based on both inferred and true
        W_dfB = @subset(W_dfB0, .!isnan.(:S_optimal_m))
        W_dfB = @subset(W_dfB, .!isnan.(:S_optimal_true_m))
        W_dfB = @subset(W_dfB, .!isnan.(:Woptimal_m)) 
        W_dfB = @subset(W_dfB, .!isnan.(:Woptimal_true_m)) 

        # check the number of dropped observations
        # B_S_union_NaN = total_mktB - size(W_dfB, 1)
        if size(W_dfB, 1) > 0
            B_positive_N_shr = sum(W_dfB.pctW_m .> 0)/size(W_dfB, 1)
            B_zero_N_shr = sum(W_dfB.pctW_m .== 0)/size(W_dfB, 1)
            B_negative_N_shr = sum(W_dfB.pctW_m .< 0)/size(W_dfB, 1)
        else
            B_positive_N_shr = missing
            B_zero_N_shr = missing
            B_negative_N_shr = missing
        end

        B_stats_table = DataFrame(
            T = total_TB, 
            total_nonNaN_obs = size(W_dfB, 1),
            avg_W_pct=round(mean(W_dfB.pctW_m), digits=3),
            min_W_pct=round(minimum(W_dfB.pctW_m), digits=3),
            st_dev_W_pct=round(std(W_dfB.pctW_m), digits=3),
            total_positive_N_pct=round(B_positive_N_shr, digits=3),
            total_zero_N_pct=round(B_zero_N_shr, digits=3),
            total_negative_N_pct=round(B_negative_N_shr, digits=3)
        )
    end
    if !isempty(W_dfC0)
        #  # drop NaN based on both inferred and true
        W_dfC = @subset(W_dfC0, .!isnan.(:S_optimal_m))
        W_dfC = @subset(W_dfC, .!isnan.(:S_optimal_true_m)) 
        W_dfC = @subset(W_dfC, .!isnan.(:Woptimal_m)) 
        W_dfC = @subset(W_dfC, .!isnan.(:Woptimal_true_m)) 
         # check the number of dropped observations
        # C_S_union_NaN = total_mktC - size(W_dfC, 1)
        if size(W_dfC, 1) > 0
            C_positive_N_shr = sum(W_dfC.pctW_m .> 0)/size(W_dfC, 1)
            C_zero_N_shr = sum(W_dfC.pctW_m .== 0)/size(W_dfC, 1)
            C_negative_N_shr = sum(W_dfC.pctW_m .< 0)/size(W_dfC, 1)
        else
            C_positive_N_shr = missing
            C_zero_N_shr = missing
            C_negative_N_shr = missing
        end
        C_stats_table = DataFrame(
            T = total_TC, 
            total_nonNaN_obs = size(W_dfC, 1),
            avg_W_pct=round(mean(W_dfC.pctW_m), digits=3),
            min_W_pct=round(minimum(W_dfC.pctW_m), digits=3),
            st_dev_W_pct=round(std(W_dfC.pctW_m), digits=3),
            total_positive_N_pct=round(C_positive_N_shr, digits=3),
            total_zero_N_pct=round(C_zero_N_shr, digits=3),
            total_negative_N_pct=round(C_negative_N_shr, digits=3)
        )
    end
        
    return B_stats_table, C_stats_table
end
