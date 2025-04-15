using Roots, Revise, Plots, LaTeXStrings
using Plots.PlotMeasures 
includet("../CES_3rd_mkt_reaction_function_plot_functions.jl")

############ intergated version of the figures #########
# Cournot - strategic substitutes
σ = 5.0
ξ1 = 1.0
ξ2 = 1.0
MC1 = 0.25
MC2 = 0.25
s1 = 0.0
s2 = 0.0
s1_subsidy = 0.35

# Q_1 series
Q_1_series = 0.0001:0.0001:1.2
Q_2_series = 0.0001:0.0001:1.2

# get R1(Q_2) series
R1_S = R1_series(Q_2_series, ξ1, ξ2, MC1, s1, σ)
R2_S = R2_series(Q_1_series, ξ1, ξ2, MC2, s2, σ)
R1_S_s = R1_series(Q_2_series, ξ1, ξ2, MC1, s1_subsidy, σ)

# get the equilibrium quantities
Q1_eq, Q2_eq = Cournot_solve_eqba_q(ξ1, ξ2, MC1, MC2, s1, s2, σ)
Q1_eq_s, Q2_eq_s = Cournot_solve_eqba_q(ξ1, ξ2, MC1, MC2, s1_subsidy, s2, σ)
# plot
# initial plot
p = plot( xlabel=L"Q_1", ylabel=L"Q_2", lw=2, xlims=(0,120), ylims=(0,120), legend = false, size = (600, 600), yguidefontrotation=-90, margin=3mm)
p = plot!([0, 100*Q1_eq], [100*Q2_eq, 100*Q2_eq], lw=1, ls=:dash, color=:black)
p = plot!([100*Q1_eq, 100*Q1_eq], [0, 100*Q2_eq], lw=1, ls=:dash, color=:black)
p = plot!([0, 100*Q1_eq_s], [100*Q2_eq_s, 100*Q2_eq_s], lw=1, ls=:dash, color=:black)
p = plot!([100*Q1_eq_s, 100*Q1_eq_s], [0, 100*Q2_eq_s], lw=1, ls=:dash, color=:black)
# add 45 degree dotted line
p = plot!([0,120], [0,120], lw=1, ls=:dot, color=:gray)
# specify the range of x and y axis from 0 to 70
p = plot!(100*Q_1_series, 100*R2_S, lw=2, color=:blue)
p = plot!(100*R1_S, 100*Q_2_series, lw=2, color=:red)
p = plot!(100*R1_S_s, 100*Q_2_series, lw=2, color=:orange)
p = scatter!([100*Q1_eq], [100*Q2_eq], color=:black)
p = scatter!([100*Q1_eq_s], [100*Q2_eq_s], color=:black)
# add labels next to the curve
p = annotate!(63,90 , text(L"R_1(Q_2)", :left, color=:red))
p = annotate!(5,30 , text(L"R_2(Q_1)", :left, color=:blue))
p = annotate!(66,8 , text(L"R_1(Q_2, s_1>0)", :left, color=:orange))
savefig(p, "Figures/Figure1a.pdf")



# cournot - strategic complements
σ = 5.0
ξ1 = 1.0
ξ2 = 1.0
MC1 = 0.25
MC2 = 0.25
MC1_s = 1.0
s1 = 0.0
s2 = 0.0
s1_tax = -0.35

# Q_1 series
Q_1_series = 0.0001:0.0001:1.0
Q_2_series = 0.0001:0.0001:1.0

# get R1(Q_2) series
R1_S = R1_series(Q_2_series, ξ1, ξ2, MC1, s1, σ)
R2_S = R2_series(Q_1_series, ξ1, ξ2, MC2, s2, σ)
R1_S_s = R1_series(Q_2_series, ξ1, ξ2, MC1, s1_tax, σ)

# get the equilibrium quantities
Q1_eq, Q2_eq = Cournot_solve_eqba_q(ξ1, ξ2, MC1, MC2, s1, s2, σ)
Q1_eq_s, Q2_eq_s = Cournot_solve_eqba_q(ξ1, ξ2, MC1, MC2, s1_tax, s2, σ)

# plot
# initial plot
p =plot( xlabel=L"Q_1", ylabel=L"Q_2", lw=2, xlims=(0,100), ylims=(0,100), legend = false, size = (600, 600), yguidefontrotation=-90, margin=3mm)
p = plot!([0, 100*Q1_eq], [100*Q2_eq, 100*Q2_eq], lw=1, ls=:dash, color=:black)
p = plot!([100*Q1_eq, 100*Q1_eq], [0, 100*Q2_eq], lw=1, ls=:dash, color=:black)
p = plot!([0, 100*Q1_eq_s], [100*Q2_eq_s, 100*Q2_eq_s], lw=1, ls=:dash, color=:black)
p = plot!([100*Q1_eq_s, 100*Q1_eq_s], [0, 100*Q2_eq_s], lw=1, ls=:dash, color=:black)
# add 45 degree dotted line
p = plot!([0,100], [0,100], lw=1, ls=:dot, color=:gray)
# specify the range of x and y axis from 0 to 70
p = plot!(100*Q_1_series, 100*R2_S, lw=2, color=:blue)
p = plot!(100*R1_S, 100*Q_2_series, lw=2, color=:red)
p = plot!(100*R1_S_s, 100*Q_2_series, lw=2, color=:orange)
p = scatter!([100*Q1_eq], [100*Q2_eq], color=:black)
p = scatter!([100*Q1_eq_s], [100*Q2_eq_s], color=:black)
# add labels next to the curve
p = annotate!(65,90 , text(L"R_1(Q_2)", :left, color=:red))
p = annotate!(5,30 , text(L"R_2(Q_1)", :left, color=:blue))
p = annotate!(20,90 , text(L"R_1(Q_2, s_1< 0)", :left, color=:orange))
savefig(p, "Figures/Figures1b.pdf")






# Bertrand 
σ = 5.0
ξ1 = 1.0
ξ2 = 1.0
MC1 = 0.25
MC2 = 0.25
s1 = 0.0
s2 = 0.0
s1_tax = -0.35

P1_eq, P2_eq = Bertrand_solve_eqba_p(ξ1, ξ2, MC1, MC2, s1, s2, σ)
P1_eq_s, P2_eq_s = Bertrand_solve_eqba_p(ξ1, ξ2, MC1, MC2, s1_tax, s2, σ)

# P_1 series
P_1_series = 0.1001:0.0001:0.6
P_2_series = 0.1001:0.0001:0.6

# get R1(P_2) series
R1_S = R1_series_p2(P_2_series, ξ1, ξ2, MC1, s1, σ)
R2_S = R2_series_p1(P_1_series, ξ1, ξ2, MC2, s2, σ)
R1_S_s = R1_series_p2(P_2_series, ξ1, ξ2, MC1, s1_tax, σ)

# plot
p = plot(xlabel=L"P_1", ylabel=L"P_2", xlims=(10,60), ylims=(10,60), legend = false, size = (600, 600), yguidefontrotation=-90, margin=3mm)
p = plot!([10, 100*P1_eq], [100*P2_eq, 100*P2_eq], lw=1, ls=:dash, color=:black)
p = plot!([100*P1_eq, 100*P1_eq], [10, 100*P2_eq], lw=1, ls=:dash, color=:black)
p = plot!([10, 100*P1_eq_s], [100*P2_eq_s, 100*P2_eq_s], lw=1, ls=:dash, color=:black)
p = plot!([100*P1_eq_s, 100*P1_eq_s], [10, 100*P2_eq_s], lw=1, ls=:dash, color=:black)
# add 45 degree dashed line
p = plot!([10,60], [10,60], lw=1, ls=:dot, color=:gray)

p = plot!(100*P_1_series, 100*R2_S, lw=2, color=:blue)
p = plot!(100*R1_S, 100*P_2_series, lw=2, color=:red)
p = plot!(100*R1_S_s, 100*P_2_series, lw=2, color=:orange)
# add the equilibrium point as a black dot, add two dashed lines to the x and y axis
p = scatter!([100*P1_eq], [100*P2_eq], color=:black)
p = scatter!([100*P1_eq_s], [100*P2_eq_s], color=:black)
# add labels next to the curve
p = annotate!(46,55 , text(L"R_1(P_2)", :left, color=:red))
p = annotate!(17,30 , text(L"R_2(P_1)", :left, color=:blue))
p = annotate!(36.5,30 , text(L"R_1(P_2, s_1<0)", :left, color=:orange))
savefig(p, "Figures/Figure1d.pdf")



σ = 5.0
ξ1 = 1.0
ξ2 = 1.0
MC1 = 0.25
MC2 = 0.25
s1_subsidy = 0.35

P1_eq, P2_eq = Bertrand_solve_eqba_p(ξ1, ξ2, MC1, MC2, s1, s2, σ)
P1_eq_s, P2_eq_s = Bertrand_solve_eqba_p(ξ1, ξ2, MC1, MC2, s1_subsidy, s2, σ)

# P_1 series
P_1_series = 0.1001:0.0001:0.6
P_2_series = 0.1001:0.0001:0.6

# get R1(P_2) series
R1_S = R1_series_p2(P_2_series, ξ1, ξ2, MC1, s1, σ)
R2_S = R2_series_p1(P_1_series, ξ1, ξ2, MC2, s2, σ)
R1_S_s = R1_series_p2(P_2_series, ξ1, ξ2, MC1, s1_subsidy, σ)

# plot
p = plot(xlabel=L"P_1", ylabel=L"P_2", xlims=(10,60), ylims=(10,60), legend = false, size = (600, 600), yguidefontrotation=-90, margin=3mm)
p = plot!([10, 100*P1_eq], [100*P2_eq, 100*P2_eq], lw=1, ls=:dash, color=:black)
p = plot!([100*P1_eq, 100*P1_eq], [10, 100*P2_eq], lw=1, ls=:dash, color=:black)
p = plot!([10, 100*P1_eq_s], [100*P2_eq_s, 100*P2_eq_s], lw=1, ls=:dash, color=:black)
p = plot!([100*P1_eq_s, 100*P1_eq_s], [10, 100*P2_eq_s], lw=1, ls=:dash, color=:black)
# add 45 degree dashed line
p = plot!([10,60], [10,60], lw=1, ls=:dot, color=:gray)

p = plot!(100*P_1_series, 100*R2_S, lw=2, color=:blue)
p = plot!(100*R1_S, 100*P_2_series, lw=2, color=:red)
p = plot!(100*R1_S_s, 100*P_2_series, lw=2, color=:orange)
# add the equilibrium point as a black dot, add two dashed lines to the x and y axis
p = scatter!([100*P1_eq], [100*P2_eq], color=:black)
p = scatter!([100*P1_eq_s], [100*P2_eq_s], color=:black)
# add labels next to the curve
p = annotate!(46,55 , text(L"R_1(P_2)", :left, color=:red))
p = annotate!(17,30 , text(L"R_2(P_1)", :left, color=:blue))
p = annotate!(22,45 , text(L"R_1(P_2, s_1>0)", :left, color=:orange))
savefig(p, "Figures/Figure1c.pdf")
