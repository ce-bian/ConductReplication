using Parameters, Optim, ForwardDiff, LinearAlgebra, Distributions, Random
using PrettyTables, DataFrames, NLsolve, NLSolversBase, RCall, ThreadsX, Revise
using CSV, JLD2, Test, Plots, StatsPlots, LaTeXStrings
includet("../../CES_structs.jl")
includet("../CES_3rd_mkt_solve_eqba_functions.jl")
includet("../CES_3rd_mkt_eqba_properties_functions.jl")
includet("../CES_3rd_mkt_optimal_s_functions.jl")
includet("../CES_3rd_mkt_s_welfare_diagram_functions.jl")

C = 2
N = 1
M = 1
σ_true = 5.0 
w = 0.95# dampening parameter: key!
MAXIT = 5000 # maximum number of iterations
TOL = 1e-12 # tolerance level
MAX_RETRIES = 300

# only third market matters
Y = [1.0, 1.0, exp(1.0)]
exp_ξ = [1.0; 1.0;;; 1.0; 1.0;;; exp(0.6); exp(1.4)]
exp_ω = [1.0; 1.0;;; 1.0; 1.0;;; 1.0; 1.0]
S = [0.0; 0.0;;; 0.0; 0.0;;; 0.0; 0.0]
τ = [1.0 2.0 1+exp(1.0);
     2.0 1.0 1+exp(1.0)]

     
IDT_m = input_dt_m_no_OG(C=C, N=N, m = 3, σ=σ_true, w=w, MAXIT=MAXIT, TOL=TOL, MAX_RETRIES=MAX_RETRIES, Y_m=Y[3], exp_ξ_m=exp_ξ[:,:,3], exp_ω_m=exp_ω[:,:,3], τ=τ, S_m=S[:,:,3])

eqB = DGP_m("B", IDT_m)
eqB.RS_m
eqC = DGP_m("C", IDT_m)
eqC.RS_m

S_list = -0.5:0.01:0.5
# transform S_list to a vector
S_list = collect(S_list)
C_S2_optimal = Cournot_optimal_s2_m_f(IDT_m)[2] 
B_S2_optimal = Bertrand_optimal_s2_m_f(IDT_m)[2] 
# add C_S2_optimal and B_S2_optimal to S_list
S_list = [S_list; C_S2_optimal; B_S2_optimal]
# sort S_list
S_list = sort(S_list)


#resultC1 = delta_W_s_list("C", IDT_m, S_list, 1)
#DTC1 = DataFrame(S1 = S_list, delta_W1 = getindex.(resultC1, 1), delta_W2 = getindex.(resultC1, 2))
resultC2 = delta_W_s_list("C", IDT_m, S_list, 2)
DTC2 = DataFrame(S2 = S_list, delta_W1 = getindex.(resultC2, 1), delta_W2 = getindex.(resultC2, 2))

#resultB1 = delta_W_s_list("B", IDT_m, S_list, 1)
#DTB1 = DataFrame(S1 = S_list, delta_W1 = getindex.(resultB1, 1), delta_W2 = getindex.(resultB1, 2))
resultB2 = delta_W_s_list("B", IDT_m, S_list, 2)
DTB2 = DataFrame(S2 = S_list, delta_W1 = getindex.(resultB2, 1), delta_W2 = getindex.(resultB2, 2))

p =plot(DTC2.S2*100, [DTC2.delta_W2*100, DTB2.delta_W2*100], 
     xlabel="", ylabel="", legend=false, framestyle = :origin,
     grid = false, ylim=(-0.20*100, 0.05*100), yticks = -20:5:5)

# get S2 tha corresponds to the maximum delta_W2
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

annotate!(p, [(-40, 0.03*100, text("Bertrand", :red, 14))])
annotate!(p, [(-30, -0.06*100, text("Cournot", :blue, 14))])
annotate!(p, [(47, -0.01*100, text(L"S (\%)", 10, :left))])
annotate!(p, [(1.0, -0.20*100, text(L"$\frac{W- W_0}{W_0} \times 100\% $", 10, :left))])
C_shr_list = [eqC.RS_m...]
C_shr_list = reverse(C_shr_list)
B_shr_list = [eqB.RS_m...]
B_shr_list = reverse(B_shr_list)
C_shr_text = "Cournot mkt shrs w/o S: " * string(round.(C_shr_list, digits=3))
B_shr_text = "Bertrand mkt shrs w/o S: " * string(round.(B_shr_list, digits=3))
annotate!(p, [(-45, -0.10*100, text(C_shr_text, 7, :left, color=:gray))])
annotate!(p, [(-45, -0.11*100, text(B_shr_text, 7, :left, color=:gray))])
savefig(p, "Figures/delta_W_bigger.pdf")







S_list = -0.5:0.01:0.2
# transform S_list to a vector
S_list = collect(S_list)
C_S1_optimal = Cournot_optimal_s1_m_f(IDT_m)[1]
B_S1_optimal = Bertrand_optimal_s1_m_f(IDT_m)[1]
# add C_S1_optimal and B_S1_optimal to S_list
S_list = [S_list; C_S1_optimal; B_S1_optimal]
# sort S_list
S_list = sort(S_list)

resultC1 = delta_W_s_list("C", IDT_m, S_list, 1)
DTC1 = DataFrame(S1 = S_list, delta_W1 = getindex.(resultC1, 1), delta_W2 = getindex.(resultC1, 2))
resultB1 = delta_W_s_list("B", IDT_m, S_list, 1)
DTB1 = DataFrame(S1 = S_list, delta_W1 = getindex.(resultB1, 1), delta_W2 = getindex.(resultB1, 2))


p =plot(DTC1.S1*100, [DTC1.delta_W1*100, DTB1.delta_W1*100], 
     xlabel="", ylabel="", legend=false, framestyle = :origin,grid=false,ylim = (-0.20*100, 0.12*100), yticks = -20:5:10)
# get S2 tha corresponds to the maximum delta_W2
C_max_value = DTC1.S1[argmax(DTC1.delta_W1)]
B_max_value = DTB1.S1[argmax(DTB1.delta_W1)]


# add vertical line at B_max_value and C_max_value
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

annotate!(p, [(-50, 0.09*100, text("Bertrand", :red, 14, :left))])
annotate!(p, [(-43, -0.06*100, text("Cournot", :blue, 14, :left))])
annotate!(p, [(19, -0.023*100, text(L"S (\%)", 10, :left))])
annotate!(p, [(1, -0.20*100, text(L"$\frac{W- W_0}{W_0} \times 100\% $", 10, :left))])
C_shr_list = [eqC.RS_m...]
C_shr_list = reverse(C_shr_list)
B_shr_list = [eqB.RS_m...]
B_shr_list = reverse(B_shr_list)
C_shr_text = "Cournot mkt shrs w/o S: " * string(round.(C_shr_list, digits=3))
B_shr_text = "Bertrand mkt shrs w/o S: " * string(round.(B_shr_list, digits=3))
annotate!(p, [(-45, -0.10*100, text(C_shr_text, 7, :left, color=:gray))])
annotate!(p, [(-45, -0.115*100, text(B_shr_text, 7, :left, color=:gray))])
savefig(p, "Figures/delta_W_smaller.pdf")