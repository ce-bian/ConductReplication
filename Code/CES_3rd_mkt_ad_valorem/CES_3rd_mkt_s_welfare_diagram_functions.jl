function W_s(R::String, IDT::input_dt_m_no_OG, S_m :: Array{Float64,2})
    IDT_m = deepcopy(IDT)
    IDT_m.S_m .= S_m
    eq = DGP_m(R, IDT_m)
    # if R == "B"
    #     Bertrand_assertion(eq)
    # elseif R == "C"
    #     Cournot_assertion(eq)
    # end
    P_m = eq.P_m
    Q_m = eq.Q_m
    MC_m = eq.MC_m
    W = P_m .* Q_m - MC_m .* Q_m
    return W # return a matrix of welfare
end


function delta_W_s(R::String, IDT::input_dt_m_no_OG, S_m :: Array{Float64,2})

    W_initial = W_s(R, IDT, IDT.S_m)
    W = W_s(R, IDT, S_m)

    return (W .- W_initial)
end


function delta_W_s_list(R::String, IDT::input_dt_m_no_OG, S_list::Vector{Float64}, firmindex::Int)
    IDT_0 = deepcopy(IDT) # initial IDT, fixed
    W_initial = W_s(R, IDT_0, IDT_0.S_m)
    delta_W_s_list = []
    for S_temp in S_list
        S = deepcopy(IDT_0.S_m) 
        S[firmindex, :] .= S_temp # change the row of S
        #delta = delta_W_s(R, IDT_0, S)
        W = W_s(R, IDT_0, S)
        # push!(delta_W_s_list, delta)
        push!(delta_W_s_list, (W .- W_initial)./W_initial)  # change: ratio
    end
    return delta_W_s_list
end




function Bertrand_assertion(eqB::eqbm_m_no_OG)
    IDT_m = eqB.IDT
    B_S1_optimal = Bertrand_optimal_s1_m_f(IDT_m)[1]   # < 0
    B_S2_optimal = Bertrand_optimal_s2_m_f(IDT_m)[2] # < 0
    B_FOC_mat = Bertrand_dπ_dp_m_single_f(eqB)
    B_SOC_mat = Bertrand_d2π_dp2_single_f(eqB)
    B_d2πx_dxdy = Bertrand_d2π_icm_dp_icm_dp_jkm_single_f(eqB, 1, 1, 2, 1)
    B_d2πy_dxdy = Bertrand_d2π_icm_dp_icm_dp_jkm_single_f(eqB, 2, 1, 1, 1)
    
    @assert all(x -> isapprox(x, 0, atol=1e-5), B_FOC_mat) # FOC, Bertrand
    @assert all(x -> x < 0, B_SOC_mat) # SOC, Bertrand
    @assert B_d2πx_dxdy > 0 # Strategic complements, Bertrand
    @assert B_d2πy_dxdy > 0 # Strategic complements, Bertrand
    @assert B_S1_optimal < 0
    @assert B_S2_optimal < 0
end


function Cournot_assertion(eqC::eqbm_m_no_OG)
    IDT_m = eqC.IDT
    C_S1_optimal = Cournot_optimal_s1_m_f(IDT_m)[1]   
    C_S2_optimal = Cournot_optimal_s2_m_f(IDT_m)[2] 
    C_FOC_mat = Cournot_dπ_dq_m_single_f(eqC)
    C_SOC_mat = Cournot_d2π_dq2_single_f(eqC)
    C_d2πx_dxdy = Cournot_d2π_icm_dq_icm_dq_jkm_single_f(eqC, 1, 1, 2, 1)
    C_d2πy_dxdy = Cournot_d2π_icm_dq_icm_dq_jkm_single_f(eqC, 2, 1, 1, 1)
    # println((eqC.RS_m[1] - 0.5) * C_S1_optimal)
    @assert all(x -> isapprox(x, 0, atol=1e-5), C_FOC_mat) # FOC, Cournot
    @assert all(x -> x < 0, C_SOC_mat) # SOC, Cournot
    @assert (eqC.RS_m[1] - 0.5) * C_d2πx_dxdy > 0 # Strategic complements/substitutes depend on market share, Cournot
    @assert (eqC.RS_m[2] - 0.5) * C_d2πy_dxdy > 0 # Strategic complements/substitutes depend on market share, Cournot
    @assert isapprox(C_d2πx_dxdy + C_d2πy_dxdy, 0, atol=1e-5) # π_xy + π_yx = 0, Cournot
    @assert (eqC.RS_m[1] - 0.5) * C_S1_optimal > 0 
    @assert (eqC.RS_m[2] - 0.5) * C_S2_optimal > 0
end

