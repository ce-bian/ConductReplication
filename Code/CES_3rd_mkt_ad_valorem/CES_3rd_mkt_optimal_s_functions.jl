
function π_icm_single_f(eq::eqbm_m_no_OG, c::Int, i::Int)
    """
    profit of firm i in market m, rely on eqbm_m only
    """
    P_m = eq.P_m
    Q_m = eq.Q_m
    MC_m = eq.MC_m
    return P_m[c,i]*Q_m[c,i] - MC_m[c,i]*Q_m[c,i]
end

# From now on assume only 3rd market model =================================
function Cournot_optimal_s1_m_f(IDT::input_dt_m_no_OG)
    """
    Optimal tax/subsidy for firm 1. Solve by fixed point iteration?
    """
    IDT_m = deepcopy(IDT)
    TOL = IDT_m.TOL
    MAXIT = IDT_m.MAXIT
    k = 1
    diff = 1.0

    S_old = [0.0, 0.0]
    while diff > TOL && k<MAXIT
        IDT_m.S_m .= S_old
        eq = DGP_m("C", IDT_m)
        P_m = eq.P_m
        Q_m = eq.Q_m
        S_new = [-Cournot_dπ_icm_dq_jkm_single_f(eq, 1, 1, 2, 1)*Cournot_d2π_icm_dq_icm_dq_jkm_single_f(eq, 2, 1, 1, 1)/(P_m[1,1]*Cournot_d2π_dq2_single_f(eq)[2,1] + Q_m[1,1]*Cournot_dp_dq_single_f(eq)[1,1]*Cournot_d2π_dq2_single_f(eq)[2,1] - Q_m[1,1]*Cournot_dp_icm_dq_jkm_single_f(eq, 1, 1, 2, 1)*Cournot_d2π_icm_dq_icm_dq_jkm_single_f(eq, 2, 1, 1, 1)), 0.0]
        diff = norm(S_new - S_old)
        S_old = w.*S_old + (1-w).*S_new
        k += 1
    end

    @assert k < MAXIT "No convergence of S after $(MAXIT) iterations"
    return S_old    
end



function Cournot_optimal_s2_m_f(IDT::input_dt_m_no_OG)
    """
    Optimal tax/subsidy for firm 2. Solve by fixed point iteration?
    """
    IDT_m = deepcopy(IDT)
    TOL = IDT_m.TOL
    MAXIT = IDT_m.MAXIT
    k = 1
    diff = 1.0

    S_old = [0.0, 0.0]
    while diff > TOL && k<MAXIT
        IDT_m.S_m .= S_old
        eq = DGP_m("C", IDT_m)
        P_m = eq.P_m
        Q_m = eq.Q_m
        S_new = [0.0, -Cournot_dπ_icm_dq_jkm_single_f(eq, 2, 1, 1, 1)*Cournot_d2π_icm_dq_icm_dq_jkm_single_f(eq, 1, 1, 2, 1)/(P_m[2,1]*Cournot_d2π_dq2_single_f(eq)[1,1] + Q_m[2,1]*Cournot_dp_dq_single_f(eq)[2,1]*Cournot_d2π_dq2_single_f(eq)[1,1] - Q_m[2,1]*Cournot_dp_icm_dq_jkm_single_f(eq, 2, 1, 1, 1)*Cournot_d2π_icm_dq_icm_dq_jkm_single_f(eq, 1, 1, 2, 1))]
        diff = norm(S_new - S_old)
        S_old = w.*S_old + (1-w).*S_new
        k += 1
    end

    @assert k < MAXIT "No convergence of S after $(MAXIT) iterations"

    return S_old    
end





function Bertrand_optimal_s1_m_f(IDT::input_dt_m_no_OG)
    """
    Optimal tax/subsidy for firm 1. Solve by fixed point iteration?
    """
    IDT_m = deepcopy(IDT)
    TOL = IDT_m.TOL
    MAXIT = IDT_m.MAXIT
    k = 1
    diff = 1.0

    S_old = [0.0, 0.0]
    while diff > TOL && k<MAXIT
        IDT_m.S_m .= S_old
        eq = DGP_m("B", IDT_m)
        # S_new = [-Bertrand_dπ_icm_dp_jkm_single_f(eq, 1, 1, 2, 1)*Bertrand_d2π_icm_dp_icm_dp_jkm_single_f(eq, 2, 1, 1, 1)/(Bertrand_dq_dp_single_f(eq)[1,1]*Bertrand_d2π_dp2_single_f(eq)[2,1])
        # , 0.0]
        P_m = eq.P_m
        Q_m = eq.Q_m
        S_new = [-Bertrand_dπ_icm_dp_jkm_single_f(eq, 1, 1, 2, 1)*Bertrand_d2π_icm_dp_icm_dp_jkm_single_f(eq, 2, 1, 1, 1)/(P_m[1,1]*Bertrand_dq_dp_single_f(eq)[1,1]*Bertrand_d2π_dp2_single_f(eq)[2,1] + Q_m[1,1]*Bertrand_d2π_dp2_single_f(eq)[2,1] - P_m[1,1]*Bertrand_dq_icm_dp_jkm_single_f(eq, 1,1,2,1)*Bertrand_d2π_icm_dp_icm_dp_jkm_single_f(eq, 2, 1, 1, 1))
        , 0.0]
        diff = norm(S_new - S_old)
        S_old = w.*S_old + (1-w).*S_new
        k += 1
    end

    @assert k < MAXIT "No convergence of S after $(MAXIT) iterations"
    return S_old    
end




function Bertrand_optimal_s2_m_f(IDT::input_dt_m_no_OG)
    """
    Optimal tax/subsidy for firm 2. Solve by fixed point iteration?
    """
    IDT_m = deepcopy(IDT)
    TOL = IDT_m.TOL
    MAXIT = IDT_m.MAXIT
    k = 1
    diff = 1.0

    S_old = [0.0, 0.0]
    while diff > TOL && k<MAXIT
        IDT_m.S_m .= S_old
        eq = DGP_m("B", IDT_m)
        P_m = eq.P_m
        Q_m = eq.Q_m
        # S_new = [0.0, -Bertrand_dπ_icm_dp_jkm_single_f(eq, 2, 1, 1, 1)*Bertrand_d2π_icm_dp_icm_dp_jkm_single_f(eq, 1, 1, 2, 1)/(Bertrand_dq_dp_single_f(eq)[2,1]*Bertrand_d2π_dp2_single_f(eq)[1,1])]
        S_new = [0.0, -Bertrand_dπ_icm_dp_jkm_single_f(eq, 2, 1, 1, 1)*Bertrand_d2π_icm_dp_icm_dp_jkm_single_f(eq, 1, 1, 2, 1)/(P_m[2,1]*Bertrand_dq_dp_single_f(eq)[2,1]*Bertrand_d2π_dp2_single_f(eq)[1,1] + Q_m[2,1]*Bertrand_d2π_dp2_single_f(eq)[1,1] - P_m[2,1]*Bertrand_dq_icm_dp_jkm_single_f(eq, 2,1,1,1)*Bertrand_d2π_icm_dp_icm_dp_jkm_single_f(eq, 1, 1, 2, 1))]
        diff = norm(S_new - S_old)
        S_old = w.*S_old + (1-w).*S_new
        k += 1
    end

    @assert k < MAXIT "No convergence of S after $(MAXIT) iterations"
    return S_old    
end


RS_m_vec(P_m, exp_ξ_m, σ) = Σ((P_m.^(1-σ)).*exp_ξ_m)  # for ForwardDiff
π_icm_p_vec(i, P_m, MC_m, exp_ξ_m, Y_m, σ) = (P_m[i]-MC_m[i])*Y_m*P_m[i]^(-σ)*exp_ξ_m[i]/RS_m_vec(P_m,exp_ξ_m,σ)


function π_p_jacobian(P_m, MC_m, exp_ξ_m, Y_m, σ)
    # Flatten the input matrices to vectors
    p_vec = vec(P_m)
    exp_ξ_vec = vec(exp_ξ_m)
    mc_vec = vec(MC_m)
    jacobian_matrix = zeros(length(p_vec), length(p_vec))
    for i in 1:length(p_vec)
        jacobian_matrix[i,:] =  ForwardDiff.gradient(p -> π_icm_p_vec(i, p, mc_vec, exp_ξ_vec, Y_m, σ), p_vec)
    end
   

    return jacobian_matrix
end

function π_p_hessian(P_m, MC_m, exp_ξ_m, Y_m, σ::Real)
    # Flatten the input matrices to vectors
    p_vec = vec(P_m)
    exp_ξ_vec = vec(exp_ξ_m)
    mc_vec = vec(MC_m)

    # Wrapper function that computes the combined values as a scalar sum
    # Compute the Hessian matrix using ForwardDiff
    hessian_matrix = zeros(size(P_m, 1),size(P_m, 1), length(p_vec))
    for i in 1:length(p_vec)
        hessian_matrix[:,:, i] = ForwardDiff.hessian(p -> π_icm_p_vec(i, p, mc_vec, exp_ξ_vec, Y_m, σ), p_vec)
    end
    return hessian_matrix
end
