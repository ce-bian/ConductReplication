using Parameters, DataFrames, RCall, Revise, DataFramesMeta, Statistics
using JLD2, LaTeXStrings
using Plots, StatsPlots
includet("../CES_simDT_structs.jl")
includet("../CES_policy_table_functions.jl")
# includet("../texout_f.jl")

# input main parameters
includet("../../../main_parameters.jl")

C = 2
N = 1
M = 1000

w = 0.7 # dampening parameter
MAXIT = 2000 # maximum number of iterations
TOL = 1e-8 # tolerance level
MAX_RETRIES = 500 # maximum number of retries for each DGP


# read in all JLD2 files starting with ...
directory = readdir("Data/")
strings_to_test = ["WdfB_WdfC","C2_N1", ".jld2"]
WdfB_WdfC_files = filter(x -> all(s -> occursin(s, x), strings_to_test), directory)


# select files with "closest"
WdfB_WdfC_closest_files = filter(x -> occursin(r"closest", x), WdfB_WdfC_files)
# select files with "critical_value"
WdfB_WdfC_hausman_files = filter(x -> occursin(r"critical_value", x), WdfB_WdfC_files)


B_policy_table = DataFrame()
C_policy_table = DataFrame()
B_policy_table_plot = DataFrame()
C_policy_table_plot = DataFrame()
# add
B_policy_table_reduced = DataFrame()
C_policy_table_reduced = DataFrame()



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

    B_policy_table_temp, C_policy_table_temp = policy_table_inferred_true_w_numT_f(W_dfB, W_dfC, T_length_B, T_length_C, M)
    rename!(B_policy_table_temp, names(B_policy_table_temp) .=> "Bertrand_" .* string.(names(B_policy_table_temp)))
    rename!(C_policy_table_temp, names(C_policy_table_temp) .=> "Cournot_" .* string.(names(C_policy_table_temp)))
    B_policy_table_temp[!, :C] .= C
    B_policy_table_temp[!, :N] .= N
    B_policy_table_temp[!, :M] .= M
    B_policy_table_temp[!, :critical_value] .= critical_value
    C_policy_table_temp[!, :C] .= C
    C_policy_table_temp[!, :N] .= N
    C_policy_table_temp[!, :M] .= M
    C_policy_table_temp[!, :critical_value] .= critical_value
    B_policy_table_temp = B_policy_table_temp[:, vcat(["C", "N", "M", "Bertrand_T","critical_value"], names(B_policy_table_temp)[2:end-4])]
    C_policy_table_temp = C_policy_table_temp[:, vcat(["C", "N", "M", "Cournot_T","critical_value"], names(C_policy_table_temp)[2:end-4])]
    for col in names(B_policy_table_temp)
        B_policy_table_temp[!, col] = convert(Vector{Union{eltype(B_policy_table_temp[!, col]), Missing}}, B_policy_table_temp[!, col])
    end
    for col in names(C_policy_table_temp)
        C_policy_table_temp[!, col] = convert(Vector{Union{eltype(C_policy_table_temp[!, col]), Missing}}, C_policy_table_temp[!, col])
    end
    append!(B_policy_table, B_policy_table_temp)
    append!(C_policy_table, C_policy_table_temp)

    B_policy_table_reduced_temp, C_policy_table_reduced_temp = policy_table_inferred_true_w_numT_reduced_f(W_dfB, W_dfC, T_length_B, T_length_C, M)
    rename!(B_policy_table_reduced_temp, names(B_policy_table_reduced_temp) .=> "Bertrand_" .* string.(names(B_policy_table_reduced_temp)))
    rename!(C_policy_table_reduced_temp, names(C_policy_table_reduced_temp) .=> "Cournot_" .* string.(names(C_policy_table_reduced_temp)))
    B_policy_table_reduced_temp[!, :C] .= C
    B_policy_table_reduced_temp[!, :N] .= N
    B_policy_table_reduced_temp[!, :M] .= M
    B_policy_table_reduced_temp[!, :critical_value] .= critical_value
    C_policy_table_reduced_temp[!, :C] .= C
    C_policy_table_reduced_temp[!, :N] .= N
    C_policy_table_reduced_temp[!, :M] .= M
    C_policy_table_reduced_temp[!, :critical_value] .= critical_value
    B_policy_table_reduced_temp = B_policy_table_reduced_temp[:, vcat(["C", "N", "M","critical_value"], names(B_policy_table_reduced_temp)[1:end-4])]
    C_policy_table_reduced_temp = C_policy_table_reduced_temp[:, vcat(["C", "N", "M", "critical_value"], names(C_policy_table_reduced_temp)[1:end-4])]
    for col in names(B_policy_table_reduced_temp)
        B_policy_table_reduced_temp[!, col] = convert(Vector{Union{eltype(B_policy_table_reduced_temp[!, col]), Missing}}, B_policy_table_reduced_temp[!, col])
    end
    for col in names(C_policy_table_reduced_temp)
        C_policy_table_reduced_temp[!, col] = convert(Vector{Union{eltype(C_policy_table_reduced_temp[!, col]), Missing}}, C_policy_table_reduced_temp[!, col])
    end
    append!(B_policy_table_reduced, B_policy_table_reduced_temp)
    append!(C_policy_table_reduced, C_policy_table_reduced_temp)

    # B_policy_table_temp, C_policy_table_temp = policy_table_inferred_true_w_numT_reduced_f(W_dfB, W_dfC, T_length_B, T_length_C, M)
    B_policy_table_temp, C_policy_table_temp = policy_table_inferred_true_w_numT_f(W_dfB, W_dfC, T_length_B, T_length_C, M)
    rename!(B_policy_table_temp, names(B_policy_table_temp) .=> "Bertrand_" .* string.(names(B_policy_table_temp)))
    rename!(C_policy_table_temp, names(C_policy_table_temp) .=> "Cournot_" .* string.(names(C_policy_table_temp)))
    B_policy_table_temp[!, :C] .= C
    B_policy_table_temp[!, :N] .= N
    B_policy_table_temp[!, :M] .= M
    B_policy_table_temp[!, :critical_value] .= critical_value
    C_policy_table_temp[!, :C] .= C
    C_policy_table_temp[!, :N] .= N
    C_policy_table_temp[!, :M] .= M
    C_policy_table_temp[!, :critical_value] .= critical_value

    B_policy_table_temp = B_policy_table_temp[:, vcat(["C", "N", "M", "Bertrand_T","critical_value"], names(B_policy_table_temp)[2:end-4])]
    C_policy_table_temp = C_policy_table_temp[:, vcat(["C", "N", "M", "Cournot_T","critical_value"], names(C_policy_table_temp)[2:end-4])]
    # B_policy_table_temp = B_policy_table_temp[:, vcat(["C", "N", "M","critical_value"], names(B_policy_table_temp)[1:end-4])]
    # C_policy_table_temp = C_policy_table_temp[:, vcat(["C", "N", "M","critical_value"], names(C_policy_table_temp)[1:end-4])]
    for col in names(B_policy_table_temp)
        B_policy_table_temp[!, col] = convert(Vector{Union{eltype(B_policy_table_temp[!, col]), Missing}}, B_policy_table_temp[!, col])
    end
    for col in names(C_policy_table_temp)
        C_policy_table_temp[!, col] = convert(Vector{Union{eltype(C_policy_table_temp[!, col]), Missing}}, C_policy_table_temp[!, col])
    end
    append!(B_policy_table_plot, B_policy_table_temp)
    append!(C_policy_table_plot, C_policy_table_temp)
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

    B_policy_table_temp, C_policy_table_temp = policy_table_inferred_true_w_numT_f(W_dfB, W_dfC, T_length_B, T_length_C, M)
    rename!(B_policy_table_temp, names(B_policy_table_temp) .=> "Bertrand_" .* string.(names(B_policy_table_temp)))
    rename!(C_policy_table_temp, names(C_policy_table_temp) .=> "Cournot_" .* string.(names(C_policy_table_temp)))
    B_policy_table_temp[!, :C] .= C
    B_policy_table_temp[!, :N] .= N
    B_policy_table_temp[!, :M] .= M
    B_policy_table_temp[!, :critical_value] .= missing
    C_policy_table_temp[!, :C] .= C
    C_policy_table_temp[!, :N] .= N
    C_policy_table_temp[!, :M] .= M
    C_policy_table_temp[!, :critical_value] .= missing
    B_policy_table_temp = B_policy_table_temp[:, vcat(["C", "N", "M", "Bertrand_T","critical_value"], names(B_policy_table_temp)[2:end-4])]
    C_policy_table_temp = C_policy_table_temp[:, vcat(["C", "N", "M", "Cournot_T","critical_value"], names(C_policy_table_temp)[2:end-4])]
    for col in names(B_policy_table_temp)
        B_policy_table_temp[!, col] = convert(Vector{Union{eltype(B_policy_table_temp[!, col]), Missing}}, B_policy_table_temp[!, col])
    end
    for col in names(C_policy_table_temp)
        C_policy_table_temp[!, col] = convert(Vector{Union{eltype(C_policy_table_temp[!, col]), Missing}}, C_policy_table_temp[!, col])
    end
    append!(B_policy_table, B_policy_table_temp)
    append!(C_policy_table, C_policy_table_temp)

    B_policy_table_reduced_temp, C_policy_table_reduced_temp = policy_table_inferred_true_w_numT_reduced_f(W_dfB, W_dfC, T_length_B, T_length_C, M)
    rename!(B_policy_table_reduced_temp, names(B_policy_table_reduced_temp) .=> "Bertrand_" .* string.(names(B_policy_table_reduced_temp)))
    rename!(C_policy_table_reduced_temp, names(C_policy_table_reduced_temp) .=> "Cournot_" .* string.(names(C_policy_table_reduced_temp)))
    B_policy_table_reduced_temp[!, :C] .= C
    B_policy_table_reduced_temp[!, :N] .= N
    B_policy_table_reduced_temp[!, :M] .= M
    B_policy_table_reduced_temp[!, :critical_value] .= missing
    C_policy_table_reduced_temp[!, :C] .= C
    C_policy_table_reduced_temp[!, :N] .= N
    C_policy_table_reduced_temp[!, :M] .= M
    C_policy_table_reduced_temp[!, :critical_value] .= missing
    B_policy_table_reduced_temp = B_policy_table_reduced_temp[:, vcat(["C", "N", "M","critical_value"], names(B_policy_table_reduced_temp)[1:end-4])]
    C_policy_table_reduced_temp = C_policy_table_reduced_temp[:, vcat(["C", "N", "M","critical_value"], names(C_policy_table_reduced_temp)[1:end-4])]
    for col in names(B_policy_table_reduced_temp)
        B_policy_table_reduced_temp[!, col] = convert(Vector{Union{eltype(B_policy_table_reduced_temp[!, col]), Missing}}, B_policy_table_reduced_temp[!, col])
    end
    for col in names(C_policy_table_reduced_temp)
        C_policy_table_reduced_temp[!, col] = convert(Vector{Union{eltype(C_policy_table_reduced_temp[!, col]), Missing}}, C_policy_table_reduced_temp[!, col])
    end
    append!(B_policy_table_reduced, B_policy_table_reduced_temp)
    append!(C_policy_table_reduced, C_policy_table_reduced_temp)


    # B_policy_table_temp, C_policy_table_temp = policy_table_inferred_true_w_numT_reduced_f(W_dfB, W_dfC, T_length_B, T_length_C, M)
    B_policy_table_temp, C_policy_table_temp = policy_table_inferred_true_w_numT_f(W_dfB, W_dfC, T_length_B, T_length_C, M)
    
    rename!(B_policy_table_temp, names(B_policy_table_temp) .=> "Bertrand_" .* string.(names(B_policy_table_temp)))
    rename!(C_policy_table_temp, names(C_policy_table_temp) .=> "Cournot_" .* string.(names(C_policy_table_temp)))
    B_policy_table_temp[!, :C] .= C
    B_policy_table_temp[!, :N] .= N
    B_policy_table_temp[!, :M] .= M
    B_policy_table_temp[!, :critical_value] .= missing
    C_policy_table_temp[!, :C] .= C
    C_policy_table_temp[!, :N] .= N
    C_policy_table_temp[!, :M] .= M
    C_policy_table_temp[!, :critical_value] .= missing

    B_policy_table_temp = B_policy_table_temp[:, vcat(["C", "N", "M", "Bertrand_T","critical_value"], names(B_policy_table_temp)[2:end-4])]
    C_policy_table_temp = C_policy_table_temp[:, vcat(["C", "N", "M", "Cournot_T","critical_value"], names(C_policy_table_temp)[2:end-4])]
    # B_policy_table_temp = B_policy_table_temp[:, vcat(["C", "N", "M", "critical_value"], names(B_policy_table_temp)[1:end-4])]
    # C_policy_table_temp = C_policy_table_temp[:, vcat(["C", "N", "M", "critical_value"], names(C_policy_table_temp)[1:end-4])]
    for col in names(B_policy_table_temp)
        B_policy_table_temp[!, col] = convert(Vector{Union{eltype(B_policy_table_temp[!, col]), Missing}}, B_policy_table_temp[!, col])
    end
    for col in names(C_policy_table_temp)
        C_policy_table_temp[!, col] = convert(Vector{Union{eltype(C_policy_table_temp[!, col]), Missing}}, C_policy_table_temp[!, col])
    end
    append!(B_policy_table_plot, B_policy_table_temp)
    append!(C_policy_table_plot, C_policy_table_temp)
end

# sort B_policy_table and C_policy_table by M and critical_value
B_policy_table = sort(B_policy_table, [:M, :critical_value])
C_policy_table = sort(C_policy_table, [:M, :critical_value])

B_policy_table_plot = sort(B_policy_table_plot, [:M, :critical_value])
C_policy_table_plot = sort(C_policy_table_plot, [:M, :critical_value])

B_policy_table_reduced = sort(B_policy_table_reduced, [:M, :critical_value])
C_policy_table_reduced = sort(C_policy_table_reduced, [:M, :critical_value])

# keep subset of B_policy_table with M >= 50
B_policy_table_sub = @subset(B_policy_table, :M .>= 50)
C_policy_table_sub = @subset(C_policy_table, :M .>= 50)
# fill in critical_value = "closest_est" for missing values
B_policy_table_sub[!, :critical_value] = coalesce.(B_policy_table_sub[!, :critical_value], "nearest")
C_policy_table_sub[!, :critical_value] = coalesce.(C_policy_table_sub[!, :critical_value], "nearest")
# sort B_policy_table and C_policy_table by M and critical_value
B_policy_table_sub = sort(B_policy_table_sub, [:M])
C_policy_table_sub = sort(C_policy_table_sub, [:M])
# print(names(B_policy_table))


B_policy_table_reduced_sub = @subset(B_policy_table_reduced, :M .>= 50)
C_policy_table_reduced_sub = @subset(C_policy_table_reduced, :M .>= 50)
# fill in critical_value = "closest_est" for missing values
B_policy_table_reduced_sub[!, :critical_value] = coalesce.(B_policy_table_reduced_sub[!, :critical_value], "nearest")
C_policy_table_reduced_sub[!, :critical_value] = coalesce.(C_policy_table_reduced_sub[!, :critical_value], "nearest")
# sort B_policy_table and C_policy_table by M and critical_value
B_policy_table_reduced_sub = sort(B_policy_table_reduced_sub, [:M])
C_policy_table_reduced_sub = sort(C_policy_table_reduced_sub, [:M])

# full version, need to modify the code to use _reduced functions to output correct tex files
# # output
texout_new(B_policy_table_reduced_sub, "Tables/CES_MC/main_est_varyCN/policy_table_ad_valorem/B_policy_table_separate_draw_C2N1_varyM_inferred_true_reduced_filter_extremes_parameter3_retry_new.tex")
texout_new(C_policy_table_reduced_sub, "Tables/CES_MC/main_est_varyCN/policy_table_ad_valorem/C_policy_table_separate_draw_C2N1_varyM_inferred_true_reduced_filter_extremes_parameter3_retry_new.tex")

# texout_new(B_policy_table_plot, "Tables/CES_MC/main_est_varyCN/policy_table/B_policy_table_separate_draw_C2N1_varyM_inferred_true_revised_plot_filter_extremes_parameter2.tex")
# texout_new(C_policy_table_plot, "Tables/CES_MC/main_est_varyCN/policy_table/C_policy_table_separate_draw_C2N1_varyM_inferred_true_revised_plot_filter_extremes_parameter2.tex")



# simple plots 
# plot Bertrand_inferred_B_pct against M for each critical_value for B_policy_table
# keep subset of B_policy_table with M >= 50
B_policy_table_sub = @subset(B_policy_table, :M .>= 50)
p = plot()
for critical_value in unique(B_policy_table_sub.critical_value)[1:3]
    df = @subset(B_policy_table_sub, :critical_value .== critical_value)
    plot!(p, df.M, df.Bertrand_inferred_B_pct*100, label = "Hausman w critical value " * string(critical_value), seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5)
    #plot!(p, df.M, df.Bertrand_inferred_B_pct)
end
df = @subset(B_policy_table_sub, ismissing.(:critical_value))
plot!(p, df.M, df.Bertrand_inferred_B_pct*100, label = "nearest neighbor", seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5)
#plot!(p, df.M, df.Bertrand_inferred_B_pct)
plot!(p, xlabel = "Number of third markets (M)", ylabel = "% of correct inferred conduct", legend=:bottomright, xticks = [50, 100, 200, 500, 800, 1000])
savefig(p, "Figures/CES_MC/policy_analysis_ad_valorem/extended_3rd_mkt/varyM/Bertrand_inferred_B_pct_vs_M_C2N1_varyM2_filter_extremes_parameter3_retry_new.pdf")


C_policy_table_sub = @subset(C_policy_table, :M .>= 50)
p = plot()
for critical_value in unique(C_policy_table_sub.critical_value)[1:3]
    df = @subset(C_policy_table_sub, :critical_value .== critical_value)
    plot!(p, df.M, df.Cournot_inferred_C_pct*100, label = "Hausman w critical value " * string(critical_value), seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5)
    #plot!(p, df.M, df.Cournot_inferred_C_pct)
end
df = @subset(C_policy_table_sub, ismissing.(:critical_value))
plot!(p, df.M, df.Cournot_inferred_C_pct*100, label = "nearest neighbor", seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5)
#plot!(p, df.M, df.Cournot_inferred_C_pct)
plot!(p, xlabel = "Number of third markets (M)", ylabel = "% of correct inferred conduct", legend=:bottomright, xticks = [50, 100, 200, 500, 800, 1000])
savefig(p, "Figures/CES_MC/policy_analysis_ad_valorem/extended_3rd_mkt/varyM/Cournot_inferred_C_pct_vs_M_C2N1_varyM2_filter_extremes_parameter3_retry_new.pdf")






B_policy_table_sub = @subset(B_policy_table, :M .>= 50)
B_policy_table_plot_sub = @subset(B_policy_table_plot, :M .>= 50)
p = plot()
for critical_value in unique(B_policy_table_sub.critical_value)[1:3]
    df = @subset(B_policy_table_sub, :critical_value .== critical_value)
    plot!(p, df.M, df.Bertrand_avg_W_pct*100, label = "Hausman w critical value " * string(critical_value), seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5)
    #plot!(p, df.M, df.Bertrand_inferred_B_pct)
end
df = @subset(B_policy_table_sub, ismissing.(:critical_value))
plot!(p, df.M, df.Bertrand_avg_W_pct*100, label = "nearest neighbor", seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5)

# get unique observations for M and Bertrand_avg_true_W_pct
df = unique(B_policy_table_plot_sub[:, [:M, :Bertrand_avg_true_W_pct]])
plot!(p, df.M, df.Bertrand_avg_true_W_pct*100, linewidth = 2, label = "optimal value based on full sample")
plot!(p, xlabel = "Number of third markets (M)", ylabel = L"avg $\frac{W-W_0}{W_0} \% $", legend=:bottomright, xticks = [50, 100, 200, 500, 800, 1000])
savefig(p, "Figures/CES_MC/policy_analysis_ad_valorem/extended_3rd_mkt/varyM/Bertrand_avg_delta_W_pct_vs_M_C2N1_varyM_full2_filter_extremes_parameter3_retry_new.pdf")


C_policy_table_sub = @subset(C_policy_table, :M .>= 50)
C_policy_table_plot_sub = @subset(C_policy_table_plot, :M .>= 50)
p = plot()
for critical_value in unique(C_policy_table_sub.critical_value)[1:3]
    df = @subset(C_policy_table_sub, :critical_value .== critical_value)
    plot!(p, df.M, df.Cournot_avg_W_pct*100, label = "Hausman w critical value " * string(critical_value), seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5)
    #plot!(p, df.M, df.Cournot_inferred_C_pct)
end
df = @subset(C_policy_table_sub, ismissing.(:critical_value))
plot!(p, df.M, df.Cournot_avg_W_pct*100, label = "nearest neighbor", seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5)

df = unique(C_policy_table_plot_sub[:, [:M, :Cournot_avg_true_W_pct]])
plot!(p, df.M, df.Cournot_avg_true_W_pct*100, linewidth = 2, label = "optimal value based on full sample")
plot!(p, xlabel = "Number of third markets (M)", ylabel = L"avg $\frac{W-W_0}{W_0} \% $ ", legend=:bottomright, xticks = [50, 100, 200, 500, 800, 1000], ylim = (-0.1, 2.7),
yticks = ([0, 0.5, 1.0, 1.5, 2.0, 2.5], ["0", "0.5", "1.0", "1.5", "2.0", "2.5"]))
savefig(p, "Figures/CES_MC/policy_analysis_ad_valorem/extended_3rd_mkt/varyM/Cournot_avg_delta_W_pct_vs_M_C2N1_varyM_full2_filter_extremes_parameter3_retry_new.pdf")



# # dots with min
# print(names(B_policy_table_sub))
# B_policy_table_sub = @subset(B_policy_table, :M .>= 50)
# B_policy_table_plot_sub = @subset(B_policy_table_plot, :M .>= 50)
# p = plot()
# for critical_value in unique(B_policy_table_sub.critical_value)[1:3]
#     df = @subset(B_policy_table_sub, :critical_value .== critical_value)
#     plot!(p, df.M, df.Bertrand_min_W_pct*100, label = "Hausman w critical value " * string(critical_value), seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5)
#     #plot!(p, df.M, df.Bertrand_inferred_B_pct)
# end
# df = @subset(B_policy_table_sub, ismissing.(:critical_value))
# plot!(p, df.M, df.Bertrand_min_W_pct*100, label = "nearest neighbor", seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5)

# # get unique observations for M and Bertrand_min_W_pct
# #df = unique(B_policy_table_plot_sub[:, [:M, :Bertrand_min_W_pct]])
# #plot!(p, df.M, df.Bertrand_min_W_pct*100, linewidth = 2, label = "optimal value based on full sample")
# plot!(p, xlabel = "Number of third markets", ylabel = L"min $\frac{W-W_0}{W_0} \% $", legend=:bottomright, xticks = [50, 100, 200, 500, 800, 1000])
# savefig(p, "Figures/CES_MC/policy_analysis_ad_valorem/extended_3rd_mkt/varyM/Bertrand_min_delta_W_pct_vs_M_C2N1_varyM_full2_filter_extremes_parameter3_new.pdf")


# # dots with min
# print(names(C_policy_table_sub))
# C_policy_table_sub = @subset(C_policy_table, :M .>= 50)
# C_policy_table_plot_sub = @subset(C_policy_table_plot, :M .>= 50)
# p = plot()
# for critical_value in unique(C_policy_table_sub.critical_value)[1:3]
#     df = @subset(C_policy_table_sub, :critical_value .== critical_value)
#     plot!(p, df.M, df.Cournot_min_W_pct*100, label = "Hausman w critical value " * string(critical_value), seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5)
#     #plot!(p, df.M, df.Bertrand_inferred_B_pct)
# end
# df = @subset(C_policy_table_sub, ismissing.(:critical_value))
# plot!(p, df.M, df.Cournot_min_W_pct*100, label = "nearest neighbor", seriestype = :scatterpath, markershape=:circle, line = (:solid, 1), alpha = 0.5)

# # get unique observations for M and Bertrand_min_W_pct
# #df = unique(B_policy_table_plot_sub[:, [:M, :Bertrand_min_W_pct]])
# #plot!(p, df.M, df.Bertrand_min_W_pct*100, linewidth = 2, label = "optimal value based on full sample")
# plot!(p, xlabel = "Number of third markets", ylabel = L"min $\frac{W-W_0}{W_0} \% $", legend=:right, xticks = [50, 100, 200, 500, 800, 1000])
# savefig(p, "Figures/CES_MC/policy_analysis/policy_analysis_ad_valorem/varyM/Cournot_min_delta_W_pct_vs_M_C2N1_varyM_full2_filter_extremes_parameter3_new.pdf")