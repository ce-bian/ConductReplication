using Parameters, Optim, ForwardDiff, LinearAlgebra, Distributions, Random
using PrettyTables, DataFrames, NLsolve, NLSolversBase, RCall, ThreadsX, Revise
using CSV, JLD2, Test, Plots, StatsPlots, LaTeXStrings

includet("../../CES_structs.jl")
includet("../CES_solve_eqba_m_functions.jl")
includet("../CES_eqba_properties_m_functions.jl")
includet("../CES_welfare_m_functions.jl")
includet("../CES_optimal_s_det_m_functions.jl")


C = 2
N = 2
M = 1
σ_true = 5.0 
w = 0.9 # dampening parameter: key!
MAXIT = 2000 # maximum number of iterations
TOL = 1e-8 # tolerance level
MAX_RETRIES = 300


Y_m = 1.0
exp_ξ_m = [1.0; 2.0;; 1.0; 2.0]
exp_ω_m = [1.0; 1.0;; 1.0; 1.0]
S_m = [0.0,0.0]
τ = [1.0 2.0 10.0;
     2.0 1.0 5.0]
exp_ξ0_m = 1.0
exp_ω0_m = 5.0

IDT_m = input_dt_m(C=C, N=N, m = C+M, σ=σ_true, w=w, MAXIT=MAXIT, TOL=TOL, MAX_RETRIES=MAX_RETRIES, Y_m=Y_m, exp_ξ_m=exp_ξ_m, exp_ξ0_m= exp_ξ0_m, exp_ω_m=exp_ω_m,exp_ω0_m=exp_ω0_m, τ=τ, S_m=S_m)

eqB = DGP_m("B", IDT_m)
eqC = DGP_m("C", IDT_m)
eqB.RS_m
eqC.RS_m
eqC.RS0_m
C_shr_list = [[vec(eqC.RS_m'); eqC.RS0_m]...] 
B_shr_list = [[vec(eqB.RS_m'); eqB.RS0_m]...]
C_S2_optimal = Cournot_optimal_s_m_c_f(IDT_m, 2)
B_S2_optimal = Bertrand_optimal_s_m_c_f(IDT_m, 2)
S_list = -0.7:0.01:0.6
# transform S_list to a vector
S_list = collect(S_list)
S_list = [S_list; C_S2_optimal; B_S2_optimal]
S_list = sort(S_list)

resultC2 = delta_W_s_list("C", IDT_m, S_list, 2)
resultC2_sum = [mapslices(sum, resultC2[i], dims=2) for i in 1:length(resultC2)]
DTC2 = DataFrame(S2 = S_list, delta_W1 = getindex.(resultC2_sum, 1), delta_W2 = getindex.(resultC2_sum, 2))
resultB2 = delta_W_s_list("B", IDT_m, S_list, 2)
resultB2_sum = [mapslices(sum, resultB2[i], dims=2) for i in 1:length(resultB2)]
DTB2 = DataFrame(S2 = S_list, delta_W1 = getindex.(resultB2_sum, 1), delta_W2 = getindex.(resultB2_sum, 2))

p =plot(DTC2.S2*100, [DTC2.delta_W2*100, DTB2.delta_W2*100], 
     xlabel="", ylabel="", legend=false, framestyle = :origin, ylims=(-0.4*100, 0.4*100),
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
annotate!(p, [(-45, 0.35*100, text("Bertrand", :red, 14))])
annotate!(p, [(25, -0.07*100, text("Cournot", :blue, 14))])
annotate!(p, [(60, -0.026*100, text(L"S (\%)", 10))])
annotate!(p, [(16, -0.4*100, text(L"$\frac{W- W_0}{W_0} \times 100\% $", 10))])
C_shr_text = "Cournot mkt shrs w/o S: " * string(round.(C_shr_list, digits=3))
B_shr_text = "Bertrand mkt shrs w/o S: " * string(round.(B_shr_list, digits=3))
annotate!(p, [(-40, -0.13*100, text(C_shr_text, 7, color=:gray))])
annotate!(p, [(-40, -0.16*100, text(B_shr_text, 7, color=:gray))])
savefig(p, "Figures/C2N2_delta_Wbigger_example1.pdf")


C_S1_optimal =Cournot_optimal_s_m_c_f(IDT_m, 1)
B_S1_optimal = Bertrand_optimal_s_m_c_f(IDT_m, 1)
S_list = -0.02:0.0001:0.02
S_list = collect(S_list)
S_list = [S_list; C_S1_optimal; B_S1_optimal]
S_list = sort(S_list)
resultC1 = delta_W_s_list("C", IDT_m, S_list, 1)
resultC1_sum = [mapslices(sum, resultC1[i], dims=2) for i in 1:length(resultC1)]
DTC1 = DataFrame(S1 = S_list, delta_W1 = getindex.(resultC1_sum, 1), delta_W2 = getindex.(resultC1_sum, 2))
resultB1 = delta_W_s_list("B", IDT_m, S_list, 1)
resultB1_sum = [mapslices(sum, resultB1[i], dims=2) for i in 1:length(resultB1)]
DTB1 = DataFrame(S1 = S_list, delta_W1 = getindex.(resultB1_sum, 1), delta_W2 = getindex.(resultB1_sum, 2))

p =plot(DTC1.S1*100, [DTC1.delta_W1*100, DTB1.delta_W1*100], 
     xlabel="", ylabel="", legend=false, framestyle = :origin, ylim = (-0.0025*100, 0.0015*100),
     grid = false)   
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
annotate!(p, [(-1.5, 0.001*100, text("Bertrand", :red, 14))])
annotate!(p, [(1.0, -0.0016*100, text("Cournot", :blue, 14))])
annotate!(p, [(2.03, -0.00027*100, text(L"S (\%)", 10))])
annotate!(p, [(0.45, -0.0025*100, text(L"$\frac{W- W_0}{W_0} \times 100\% $", 10))])

C_shr_list = [[vec(eqC.RS_m'); eqC.RS0_m]...] 
B_shr_list = [[vec(eqB.RS_m'); eqB.RS0_m]...]
C_shr_text = "Cournot mkt shrs w/o S: " * string(round.(C_shr_list, digits=3))
B_shr_text = "Bertrand mkt shrs w/o S: " * string(round.(B_shr_list, digits=3))
annotate!(p, [(-1.25, -0.0005*100, text(C_shr_text, 7, color=:gray))])
annotate!(p, [(-1.25, -0.0008*100, text(B_shr_text, 7, color=:gray))])
savefig(p, "Figures/C2N2_delta_Wsmaller_example1.pdf")

