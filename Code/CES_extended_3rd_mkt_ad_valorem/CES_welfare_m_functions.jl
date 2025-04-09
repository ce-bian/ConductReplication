function W_s_c(R::String, IDT::input_dt_m, S_cm :: Float64; test::Bool=false)
    """
    Calculate the welfare of production origin c
    """
    IDT_m = deepcopy(IDT)
    IDT_m.S_m[c] = S_cm # update the subsidy/tax vector
    # initiate an eqbm_m object
    C = IDT_m.C
    N = IDT_m.N
    eq = eqbm_m(IDT, R, zeros(C, N), 0.0, zeros(C, N), 0.0, zeros(C, N), 0.0, zeros(C, N), 0.0, 0.0, zeros(C, N), 0.0)
    W_c = NaN
    try
        eq = DGP_m(R, IDT_m)
        if test
            if R == "B"
                Bertrand_assertion(eq)
            elseif R == "C"
                Cournot_assertion(eq)
            end
        end
        # focus on the welfare of foreign firms in the destination market
        P_cm = eq.P_m[c, :]  # C*N matrix
        Q_cm = eq.Q_m[c, :]  # C*N matrix
        MC_cm = eq.MC_m[c, :] # C*N matrix, including τ and ω
        # expand S_m to a matrix by repeating the vector
        S_cm = fill(S_cm, IDT_m.N) # C*N matrix, every element in each row is the same

        # W_c = sum(P_cm .* Q_cm - MC_cm .* Q_cm - S_cm .* Q_cm)
        W_c = sum(P_cm .* Q_cm - MC_cm .* Q_cm)


    catch
        W_c = NaN # return NaN if the equilibrium does not exist, no convergence when solving equilibrium prices
    end

    return W_c # return a scalar of welfare of production origin c
end


function W_s_c_given_eq(R::String, eq_m::eqbm_m; test::Bool=false)
    """
    Calculate the welfare of production origin c given the equilibrium eq
    """
    eq = deepcopy(eq_m)
    if test
        if R == "B"
            Bertrand_assertion(eq)
        elseif R == "C"
            Cournot_assertion(eq)
        end
    end
    # focus on the welfare of foreign firms in the destination market
    P_cm = eq.P_m[c,:]  # C*N matrix
    Q_cm = eq.Q_m[c,:]  # C*N matrix
    MC_cm = eq.MC_m[c,:] # C*N matrix, including τ and ω
    # expand S_m to a matrix by repeating the vector
    # S_cm = fill(S_cm, IDT_m.N) # C*N matrix, every element in each row is the same
    S_cm = fill(eq.IDT.S_m[c], eq.IDT.N) # C*N matrix, every element in each row is the same
    # W_c = sum(P_cm .* Q_cm - MC_cm .* Q_cm - S_cm .* Q_cm)
    W_c = sum(P_cm .* Q_cm - MC_cm .* Q_cm )
    return W_c # return a scalar of welfare of production origin c
end




function W_s(R::String, IDT::input_dt_m, S_m :: Vector{Float64}; test::Bool=false)
    IDT_m = deepcopy(IDT)
    IDT_m.S_m .= S_m # update the subsidy/tax vector
    eq = DGP_m(R, IDT_m)
    if test
        if R == "B"
            Bertrand_assertion(eq)
        elseif R == "C"
            Cournot_assertion(eq)
        end
    end
    # focus on the welfare of foreign firms in the destination market
    P_m = eq.P_m  # C*N matrix
    Q_m = eq.Q_m  # C*N matrix
    MC_m = eq.MC_m # C*N matrix, including τ and ω
    # expand S_m to a matrix by repeating the vector
    S_m = repeat(S_m, 1, IDT_m.N) # C*N matrix, every element in each row is the same

    W = P_m .* Q_m - MC_m .* Q_m
    return W # return a matrix of welfare (C*N)
end


function delta_W_s_list(R::String, IDT::input_dt_m, S_list::Vector{Float64}, Cindex::Int; test::Bool=false)
    IDT_0 = deepcopy(IDT) # initial IDT, fixed
    W_initial = W_s(R, IDT_0, IDT_0.S_m, test=test)
    delta_W_s_list = []
    for S_temp in S_list
        S = deepcopy(IDT_0.S_m)  # a vector
        S[Cindex] = S_temp # change one element of S
        #delta = delta_W_s(R, IDT_0, S)
        W = W_s(R, IDT_0, S)
        # push!(delta_W_s_list, delta)
        push!(delta_W_s_list, (W .- W_initial)./W_initial)  # change: ratio
    end
    return delta_W_s_list
end



function Bertrand_assertion(eqB::eqbm_m)
    B_FOC_mat = Bertrand_dπ_dp_m_single_f(eqB)
    B_FOC_0 = Bertrand_dπ0_dp0_m_single_f(eqB)
    B_SOC_mat = Bertrand_d2π_dp2_single_f(eqB)
    B_SOC_0 = Bertrand_d2π0_dp02_m_single_f(eqB)
    B_d2π_icm_dp_icm_dp_jkm = [Bertrand_d2π_icm_dp_icm_dp_jkm_single_f(eqB, c, i, k, j) for c in 1:eqB.IDT.C, i in 1:eqB.IDT.N, k in 1:eqB.IDT.C, j in 1:eqB.IDT.N]
    B_dπ_icm_dp_0mm_single_f = [Bertrand_dπ_icm_dp_0mm_single_f(eqB, c, i) for c in 1:eqB.IDT.C, i in 1:eqB.IDT.N]
    B_dπ_0mm_dp_jkm_single_f = [Bertrand_dπ_0mm_dp_jkm_single_f(eqB, k, j) for k in 1:eqB.IDT.C, j in 1:eqB.IDT.N]
    
    @assert all(x -> isapprox(x, 0, atol=1e-5), B_FOC_mat) # FOC, Bertrand, inside goods
    @assert isapprox(B_FOC_0, 0, atol=1e-5) # FOC, Bertrand, outside good
    @assert all(x -> x < 0, B_SOC_mat) # SOC, Bertrand, inside goods
    @assert B_SOC_0 < 0 # SOC, Bertrand, outside good
    @assert all(x -> x > 0, B_d2π_icm_dp_icm_dp_jkm) # Strategic complements, Bertrand
    @assert all(x -> x > 0, B_dπ_icm_dp_0mm_single_f) # Strategic complements, Bertrand
    @assert all(x -> x > 0, B_dπ_0mm_dp_jkm_single_f) # Strategic complements, Bertrand      
end

function Cournot_assertion(eqC::eqbm_m)
    C_FOC_mat = Cournot_dπ_dq_m_single_f(eqC)
    C_FOC_0 = Cournot_dπ0_dq0_m_single_f(eqC)
    C_SOC_mat = Cournot_d2π_dq2_m_single_f(eqC)
    C_SOC_0 = Cournot_d2π0_dq02_m_single_f(eqC)
    # println((eqC.RS_m[1] - 0.5) * C_S1_optimal)
    @assert all(x -> isapprox(x, 0, atol=1e-5), C_FOC_mat) # FOC, Cournot
    @assert isapprox(C_FOC_0, 0, atol=1e-5) # FOC, Cournot, outside good
    @assert all(x -> x < 0, C_SOC_mat) # SOC, Cournot
    @assert C_SOC_0 < 0 # SOC, Cournot, outside good
end