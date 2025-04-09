using Parameters, DataFrames, RCall, Revise, DataFramesMeta, Statistics
using JLD2, LaTeXStrings
using Plots, StatsPlots
includet("../CES_simDT_structs.jl")
includet("../CES_policy_table_functions.jl")
includet("../texout_f.jl")

############################## parameter 3 ################################
seed = 1224
C = 2
N = 1
M = 1000
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
σ_true = 4.0  # modified 
w = 0.95 # dampening parameter
MAXIT = 5000 # maximum number of iterations
TOL = 1e-8 # tolerance level
MAX_RETRIES = 500 # maximum number of retries for each DGP

T = 50
B = 250


# read in all JLD2 files starting with ...
directory = readdir("Data/Out/CES_MC/main_est/policy_table_ad_valorem/C2N1_varyM/")
strings_to_test = ["Muomega01.5", "Vomega00.1", "Muxi0-3.0", "Vpsi2.0", "parameter3", "retry.jld2"]
WdfB_WdfC_files = filter(x -> all(s -> occursin(s, x), strings_to_test), directory)

# only keep M > 200
WdfB_WdfC_files = filter!(file -> parse(Int, match(r"M(\d+)_", file).captures[1]) > 200, WdfB_WdfC_files)



directory = readdir("Data/Out/CES_MC/main_est/policy_table_ad_valorem/C2N1_varyM/")
strings_to_test = ["Muomega01.5", "Vomega00.1", "Muxi0-3.0", "Vpsi2.0", "parameter3", "retry_new.jld2"]
WdfB_WdfC_files_new = filter(x -> all(s -> occursin(s, x), strings_to_test), directory)
# only keep M <= 200
WdfB_WdfC_files_new = filter!(file -> parse(Int, match(r"M(\d+)_", file).captures[1]) <= 200, WdfB_WdfC_files_new)

# merge the two lists
WdfB_WdfC_files = vcat(WdfB_WdfC_files, WdfB_WdfC_files_new)

# select files with "closest"
WdfB_WdfC_closest_files = filter(x -> occursin(r"closest", x), WdfB_WdfC_files)
# select files with "critical_value"
WdfB_WdfC_hausman_files = filter(x -> occursin(r"critical_value", x), WdfB_WdfC_files)


M_list = [50, 100, 200, 500, 800, 1000]

# file = WdfB_WdfC_hausman_files[1]
# W_dfB, W_dfC = load("Data/Out/CES_MC/main_est/policy_table/alternative_revised/C2N1_varyM/$file", "W_dfB", "W_dfC")
# names(W_dfB)
# "true_conduct"
# "inferred_conduct"
# "t"
# "c"
# "m"
# "S_optimal_m"
# "W0_m"
# "Woptimal_m"
# "deltaW_m"
# "pctW_m"
# "S_optimal_true_m"
# "W0_true_m"
# "Woptimal_true_m"
# "deltaW_true_m"
# "pctW_true_m"
B_stats_table = DataFrame()
C_stats_table = DataFrame()

for file in WdfB_WdfC_hausman_files
    M = parse(Int, match(r"M(\d+)_", file).captures[1])

    InputEstFile = "Data/Out/CES_MC/main_est/C2N1_varyM/C2N1_varyM_est/new/dfB_dfC_GP_C$(C)_N$(N)_M$(M)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_parameter3_retry.jld2"
    dfB, dfC, GP = load(InputEstFile, "dfB", "dfC", "GP");
    t_list_B = unique(dfB.t)
    t_list_C = unique(dfC.t)
    T_length_B = length(t_list_B)
    T_length_C = length(t_list_C)

    critical_value = parse(Float64, match(r"critical_value(\d+\.\d+)", file).captures[1])

    W_dfB, W_dfC = load("Data/Out/CES_MC/main_est/policy_table_ad_valorem/C2N1_varyM/$file", "W_dfB", "W_dfC")

    B_stats_table_temp, C_stats_table_temp = policy_table_summary_stats(W_dfB, W_dfC, T_length_B, T_length_C, M)
    rename!(B_stats_table_temp, names(B_stats_table_temp) .=> "Bertrand_" .* string.(names(B_stats_table_temp)))
    rename!(C_stats_table_temp, names(C_stats_table_temp) .=> "Cournot_" .* string.(names(C_stats_table_temp)))
    B_stats_table_temp[!, :C] .= C
    B_stats_table_temp[!, :N] .= N
    B_stats_table_temp[!, :M] .= M
    B_stats_table_temp[!, :critical_value] .= critical_value
    C_stats_table_temp[!, :C] .= C
    C_stats_table_temp[!, :N] .= N
    C_stats_table_temp[!, :M] .= M
    C_stats_table_temp[!, :critical_value] .= critical_value
    B_stats_table_temp = B_stats_table_temp[:, vcat(["C", "N", "M", "Bertrand_T","critical_value"], names(B_stats_table_temp)[2:end-4])]
    C_stats_table_temp = C_stats_table_temp[:, vcat(["C", "N", "M", "Cournot_T","critical_value"], names(C_stats_table_temp)[2:end-4])]
    for col in names(B_stats_table_temp)
        B_stats_table_temp[!, col] = convert(Vector{Union{eltype(B_stats_table_temp[!, col]), Missing}}, B_stats_table_temp[!, col])
    end
    for col in names(C_stats_table_temp)
        C_stats_table_temp[!, col] = convert(Vector{Union{eltype(C_stats_table_temp[!, col]), Missing}}, C_stats_table_temp[!, col])
    end
    append!(B_stats_table, B_stats_table_temp)
    append!(C_stats_table, C_stats_table_temp)

end

for file in WdfB_WdfC_closest_files
    M = parse(Int, match(r"M(\d+)_", file).captures[1])

    InputEstFile = "Data/Out/CES_MC/main_est/C2N1_varyM/C2N1_varyM_est/new/dfB_dfC_GP_C$(C)_N$(N)_M$(M)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_parameter3_retry.jld2"
    dfB, dfC, GP = load(InputEstFile, "dfB", "dfC", "GP");
    t_list_B = unique(dfB.t)
    t_list_C = unique(dfC.t)
    T_length_B = length(t_list_B)
    T_length_C = length(t_list_C)

    W_dfB, W_dfC = load("Data/Out/CES_MC/main_est/policy_table_ad_valorem/C2N1_varyM/$file", "W_dfB", "W_dfC")

    B_stats_table_temp, C_stats_table_temp = policy_table_summary_stats(W_dfB, W_dfC, T_length_B, T_length_C, M)
    rename!(B_stats_table_temp, names(B_stats_table_temp) .=> "Bertrand_" .* string.(names(B_stats_table_temp)))
    rename!(C_stats_table_temp, names(C_stats_table_temp) .=> "Cournot_" .* string.(names(C_stats_table_temp)))
    B_stats_table_temp[!, :C] .= C
    B_stats_table_temp[!, :N] .= N
    B_stats_table_temp[!, :M] .= M
    B_stats_table_temp[!, :critical_value] .= missing
    C_stats_table_temp[!, :C] .= C
    C_stats_table_temp[!, :N] .= N
    C_stats_table_temp[!, :M] .= M
    C_stats_table_temp[!, :critical_value] .= missing
    B_stats_table_temp = B_stats_table_temp[:, vcat(["C", "N", "M", "Bertrand_T","critical_value"], names(B_stats_table_temp)[2:end-4])]
    C_stats_table_temp = C_stats_table_temp[:, vcat(["C", "N", "M", "Cournot_T","critical_value"], names(C_stats_table_temp)[2:end-4])]
    for col in names(B_stats_table_temp)
        B_stats_table_temp[!, col] = convert(Vector{Union{eltype(B_stats_table_temp[!, col]), Missing}}, B_stats_table_temp[!, col])
    end
    for col in names(C_stats_table_temp)
        C_stats_table_temp[!, col] = convert(Vector{Union{eltype(C_stats_table_temp[!, col]), Missing}}, C_stats_table_temp[!, col])
    end
    append!(B_stats_table, B_stats_table_temp)
    append!(C_stats_table, C_stats_table_temp)
end

# keep subset of B_stats_table with M >= 50
B_stats_table = @subset(B_stats_table, :M .>= 50)
C_stats_table = @subset(C_stats_table, :M .>= 50)
# fill in critical_value = "closest_est" for missing values
B_stats_table[!, :critical_value] = coalesce.(B_stats_table[!, :critical_value], "nearest")
C_stats_table[!, :critical_value] = coalesce.(C_stats_table[!, :critical_value], "nearest")
# sort B_stats_table and C_policy_table by M and critical_value
B_stats_table = sort(B_stats_table, [:M])
C_stats_table = sort(C_stats_table, [:M])


texout_new(B_stats_table, "Tables/CES_MC/main_est_varyCN/policy_table_ad_valorem/B_policy_table_stats_separate_draw_C2N1_varyM_inferred_true_reduced_filter_extremes_parameter3_retry_new.tex")
texout_new(C_stats_table, "Tables/CES_MC/main_est_varyCN/policy_table_ad_valorem/C_policy_table_stats_separate_draw_C2N1_varyM_inferred_true_reduced_filter_extremes_parameter3_retry_new.tex")


# extract critical_value = "nearest"
B_stats_table_nearest = @subset(B_stats_table, :critical_value .== "nearest")
names(B_stats_table_nearest)
C_stats_table_nearest = @subset(C_stats_table, :critical_value .== "nearest")
names(C_stats_table_nearest)


# plot Bertrand_total_negative_N_pct against M
p = plot(B_stats_table_nearest.M, B_stats_table_nearest.Bertrand_total_negative_N_pct*100, seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5, label = "nearest neighbor", ylabel = "Share of negative welfare cases (%)", xlabel = "M", title = "Bertrand DGP")
savefig(p, "Figures/CES_MC/policy_analysis_ad_valorem/extended_3rd_mkt/varyM/Bertrand_total_negative_N_pct_C2N1_varyM_inferred_true_reduced_filter_extremes_parameter3_retry_new.pdf")


# plot Bertrand_st_dev_W_pct against M
p = plot()
for critical_value in unique(B_stats_table.critical_value)[1:3]
    df = @subset(B_stats_table, :critical_value .== critical_value)
    plot!(p, df.M, df.Bertrand_st_dev_W_pct*100, seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5, label = "Hausman w critical value " * string(critical_value))
end
plot!(p, B_stats_table_nearest.M, B_stats_table_nearest.Bertrand_st_dev_W_pct*100, seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5, label = "nearest neighbor", ylabel = "Standard deviation of welfare change (%)", xlabel = "Number of third markets (M)", title = "Bertrand DGP")
savefig(p, "Figures/CES_MC/policy_analysis_ad_valorem/extended_3rd_mkt/varyM/Bertrand_st_dev_W_pct_C2N1_varyM_inferred_true_reduced_filter_extremes_parameter3_retry_new.pdf")


p = plot()
for critical_value in unique(C_stats_table.critical_value)[1:3]
    df = @subset(C_stats_table, :critical_value .== critical_value)
    plot!(p, df.M, df.Cournot_st_dev_W_pct*100, seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5, label = "Hausman w critical value " * string(critical_value))
end
plot!(p, C_stats_table_nearest.M, C_stats_table_nearest.Cournot_st_dev_W_pct*100, seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5, label = "nearest neighbor", ylabel = "Standard deviation of welfare change (%)", xlabel = "Number of third markets (M)", title = "Cournot DGP")
savefig(p, "Figures/CES_MC/policy_analysis_ad_valorem/extended_3rd_mkt/varyM/Cournot_st_dev_W_pct_C2N1_varyM_inferred_true_reduced_filter_extremes_parameter3_retry_new.pdf")




# check the minimum and maximum values of the RS0
M0 = 1000
InputDTFile = "Data/Out/CES_MC/main_est/simDT_output/eqbaB_eqbaC_GP_seed$(seed)_C$(C)_N$(N)_M$(M0)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_filter_extremes.jld2"
eqbm_output_B, eqbm_output_C, GP = load(InputDTFile, "eqbm_output_B", "eqbm_output_C", "GP");
# retrieve each RS0_m from eqbm_output_B and eqbm_output_C, through eqbm_t_list
RS0_B = [eqbm_output_B.eqbm_t_list[i].eqbm_m_list[j].RS0_m for i in 1:length(eqbm_output_B.eqbm_t_list), j in 1:length(eqbm_output_B.eqbm_t_list[1].eqbm_m_list)]
#minimum and maximum values of RS0_B
minimum(RS0_B), mean(RS0_B), maximum(RS0_B)
RS0_C = [eqbm_output_C.eqbm_t_list[i].eqbm_m_list[j].RS0_m for i in 1:length(eqbm_output_C.eqbm_t_list), j in 1:length(eqbm_output_C.eqbm_t_list[1].eqbm_m_list)]
# minimum and maximum values of RS0_C
minimum(RS0_C), mean(RS0_C), maximum(RS0_C)


InputEstFile = "Data/Out/CES_MC/main_est/C2N1_varyM/C2N1_varyM_est/new/dfB_dfC_GP_C$(C)_N$(N)_M$(M0)_MuY$(μ_Y)_VY$(V_Y)_Muxi$(μ_ξ)_Muomega$(μ_ω)_rho$(ρ)_Vomega$(V_ω)_Vxi$(V_ξ)_Mupsi$(μ_psi)_Vpsi$(V_psi)_beta$(β)_sigma$(σ_true)_Muxi0$(μ_ξ0)_Vxi0$(V_ξ0)_Muomega0$(μ_ω0)_Vomega0$(V_ω0)_T$(T)_B$(B)_parameter3_retry.jld2"
dfB, dfC, GP = load(InputEstFile, "dfB", "dfC", "GP");
minimum(dfB.mean_RS0), mean(dfB.mean_RS0), maximum(dfB.mean_RS0), median(dfB.mean_RS0)
minimum(dfC.mean_RS0), mean(dfC.mean_RS0), maximum(dfC.mean_RS0), median(dfC.mean_RS0)