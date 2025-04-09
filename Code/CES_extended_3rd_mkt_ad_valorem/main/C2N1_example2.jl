using Parameters, Optim, ForwardDiff, LinearAlgebra, Distributions, Random
using PrettyTables, DataFrames, NLsolve, NLSolversBase, RCall, ThreadsX, Revise
using CSV, JLD2, Test, Plots, StatsPlots, LaTeXStrings

includet("../CES_extended_3rd_mkt_structs.jl")
includet("../CES_solve_eqba_m_functions.jl")
includet("../CES_eqba_properties_m_functions.jl")
includet("../CES_welfare_m_functions.jl")
includet("../CES_optimal_s_det_m_functions.jl")

C = 2 
N = 1
M = 1
σ_true = 5.0 
w = 0.9 # dampening parameter: key!
MAXIT = 2000 # maximum number of iterations
TOL = 1e-8 # tolerance level
MAX_RETRIES = 300


Y_m = exp(1.0)
exp_ξ_m =reshape([exp(0.8); exp(1.2)], :, 1)
exp_ω_m = reshape([1.0; 1.0], :, 1)
S_m = [0.0,0.0]
τ = [1.0 2.0 1+exp(1.0);
     2.0 1.0 1+exp(1.0)]
exp_ξ0_m = exp(6.0)
exp_ω0_m = 1.0

IDT_m = input_dt_m(C=C, N=N, m = C+M, σ=σ_true, w=w, MAXIT=MAXIT, TOL=TOL, MAX_RETRIES=MAX_RETRIES, Y_m=Y_m, exp_ξ_m=exp_ξ_m, exp_ξ0_m= exp_ξ0_m, exp_ω_m=exp_ω_m,exp_ω0_m=exp_ω0_m, τ=τ, S_m=S_m)

eqB = DGP_m("B", IDT_m)
eqC = DGP_m("C", IDT_m)
eqB.RS_m
eqB.RS0_m
eqC.RS_m
eqC.RS0_m
# combine eqC.RS_m and eqC.RS0_m
C_shr_list = [[eqC.RS_m; eqC.RS0_m]...] 
B_shr_list = [[eqB.RS_m; eqB.RS0_m]...]


C_S2_optimal =Cournot_optimal_s_m_c_f(IDT_m, 2)
B_S2_optimal = Bertrand_optimal_s_m_c_f(IDT_m, 2)
S_list = -0.3:0.005:0.15
S_list = collect(S_list)
S_list = [S_list; C_S2_optimal; B_S2_optimal]
S_list = sort(S_list)

resultC2 = delta_W_s_list("C", IDT_m, S_list, 2)
DTC2 = DataFrame(S2 = S_list, delta_W1 = getindex.(resultC2, 1), delta_W2 = getindex.(resultC2, 2))
resultB2 = delta_W_s_list("B", IDT_m, S_list, 2)
DTB2 = DataFrame(S2 = S_list, delta_W1 = getindex.(resultB2, 1), delta_W2 = getindex.(resultB2, 2))

p =plot(DTC2.S2*100, [DTC2.delta_W2*100, DTB2.delta_W2*100], 
     xlabel="", ylabel="", legend=false, framestyle = :origin, ylim = (-0.15*100, 0.1*100),
     grid = false)

C_max_value = DTC2.S2[argmax(DTC2.delta_W2)]
B_max_value = DTB2.S2[argmax(DTB2.delta_W2)]

# add vertical line at B_max_value and C_max_value
p= vline!([B_max_value]*100, color=:orange)
p = vline!([C_max_value]*100, color=:lightblue)
# add horizontal dashed line at DTB2.delta_W2[argmax(DTC2.delta_W2)] from 0 to C_max_value
x = [0, C_max_value*100] # x coordinates from 0 to C_max_value
y = fill(DTB2.delta_W2[argmax(DTC2.delta_W2)]*100, 2) # y coordinates at the specified height
plot!(p, x, y, color=:gray, linestyle=:dash, label="")
# add horizontal dashed line at DTC2.delta_W2[argmax(DTC2.delta_W2)] from 0 to B_max_value
x = [0, B_max_value*100] # x coordinates from 0 to B_max_value
y = fill(DTC2.delta_W2[argmax(DTB2.delta_W2)]*100, 2) # y coordinates at the specified height
plot!(p, x, y, color=:gray, linestyle=:dash, label="")
annotate!(p, [(-27, 0.05*100, text("Bertrand", :red, 14, :left))])
annotate!(p, [(4, -0.02*100, text("Cournot", :blue, 14, :left))])
annotate!(p, [(14.0, -0.01*100, text(L"S (\%)", 10, :left))])
annotate!(p, [(0.5, -0.15*100, text(L"$\frac{W- W_0}{W_0} \times 100\% $", 10, :left))]) 
C_shr_text = "Cournot mkt shrs w/o S: " * string(round.(C_shr_list, digits=3))
B_shr_text = "Bertrand mkt shrs w/o S: " * string(round.(B_shr_list, digits=3))
annotate!(p, [(-30, -0.08*100, text(C_shr_text, 7, :left, color=:gray))])
annotate!(p, [(-30, -0.09*100, text(B_shr_text, 7, :left, color=:gray))])
savefig(p, "Figures/C2N1_delta_Wbigger_revised2.pdf")


C_S1_optimal =Cournot_optimal_s_m_c_f(IDT_m, 1)
B_S1_optimal = Bertrand_optimal_s_m_c_f(IDT_m, 1)
S_list = -0.3:0.005:0.15
S_list = collect(S_list)
S_list = [S_list; C_S1_optimal; B_S1_optimal]
S_list = sort(S_list)
resultC1 = delta_W_s_list("C", IDT_m, S_list, 1)
DTC1 = DataFrame(S1 = S_list, delta_W1 = getindex.(resultC1, 1), delta_W2 = getindex.(resultC1, 2))
resultB1 = delta_W_s_list("B", IDT_m, S_list, 1)
DTB1 = DataFrame(S1 = S_list, delta_W1 = getindex.(resultB1, 1), delta_W2 = getindex.(resultB1, 2))
p =plot(DTC1.S1*100, [DTC1.delta_W1*100, DTB1.delta_W1*100], 
     xlabel="", ylabel="", legend=false, framestyle = :origin, ylim = (-0.15*100, 0.05*100),
     grid = false, yticks = -0.15*100:0.05*100:0.05*100)
C_max_value = DTC1.S1[argmax(DTC1.delta_W1)]
B_max_value = DTB1.S1[argmax(DTB1.delta_W1)]
p= vline!([B_max_value]*100, color=:orange)
p = vline!([C_max_value]*100, color=:lightblue)

# add horizontal dashed line at DTB1.delta_W1[argmax(DTC1.delta_W1)] from 0 to C_max_value
x = [0, C_max_value*100] # x coordinates from 0 to C_max_value
y = fill(DTB1.delta_W1[argmax(DTC1.delta_W1)]*100, 2) # y coordinates at the specified height
plot!(p, x, y, color=:gray, linestyle=:dash, label="")
# add horizontal dashed line at DTC1.delta_W1[argmax(DTC1.delta_W1)] from 0 to B_max_value
x = [0, B_max_value*100] # x coordinates from 0 to B_max_value
y = fill(DTC1.delta_W1[argmax(DTB1.delta_W1)]*100, 2) # y coordinates at the specified height
plot!(p, x, y, color=:gray, linestyle=:dash, label="")

annotate!(p, [(-20, 0.025*100, text("Bertrand", :red, 14, :left))])
annotate!(p, [(5, -0.025*100, text("Cournot", :blue, 14, :left))])
annotate!(p, [(14.0, -0.01*100, text(L"S (\%)", 10, :left))])
annotate!(p, [(0.5, -0.15*100, text(L"$\frac{W- W_0}{W_0} \times 100\% $", 10, :left))]) 


C_shr_list = [[eqC.RS_m; eqC.RS0_m]...] 
B_shr_list = [[eqB.RS_m; eqB.RS0_m]...]
C_shr_text = "Cournot mkt shrs w/o S: " * string(round.(C_shr_list, digits=3))
B_shr_text = "Bertrand mkt shrs w/o S: " * string(round.(B_shr_list, digits=3))
annotate!(p, [(-30, -0.08*100, text(C_shr_text, 7, :left, color=:gray))])
annotate!(p, [(-30, -0.09*100, text(B_shr_text, 7, :left, color=:gray))])
savefig(p, "Figures/C2N1_delta_Wsmaller_revised2.pdf")
