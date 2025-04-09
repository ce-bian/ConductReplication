# functions for 3rd mkt model only
# market size is normalized as 1
P_1(Q1::Float64, Q2::Float64, ξ1::Float64, ξ2::Float64, σ::Float64) = Q1^(-1/σ)*ξ1^((σ-1)/σ)/(Q1^((σ-1)/σ)*ξ1^((σ-1)/σ) + Q2^((σ-1)/σ)*ξ2^((σ-1)/σ))
P_2(Q1::Float64, Q2::Float64, ξ1::Float64, ξ2::Float64, σ::Float64) = Q2^(-1/σ)*ξ2^((σ-1)/σ)/(Q1^((σ-1)/σ)*ξ1^((σ-1)/σ) + Q2^((σ-1)/σ)*ξ2^((σ-1)/σ))

RS_1(Q1::Float64, Q2::Float64, ξ1::Float64, ξ2::Float64, σ::Float64) = P_1(Q1, Q2, ξ1, ξ2, σ)*Q1/(P_1(Q1, Q2, ξ1, ξ2, σ)*Q1 + P_2(Q1, Q2, ξ1, ξ2, σ)*Q2)
RS_2(Q1::Float64, Q2::Float64, ξ1::Float64, ξ2::Float64, σ::Float64) = P_2(Q1, Q2, ξ1, ξ2, σ)*Q2/(P_1(Q1, Q2, ξ1, ξ2, σ)*Q1 + P_2(Q1, Q2, ξ1, ξ2, σ)*Q2)

Q_1(P1::Float64, P2::Float64, ξ1::Float64, ξ2::Float64, σ::Float64) = P1^(-σ)*ξ1^(σ-1)/(P1^(1-σ)*ξ1^(σ-1) + P2^(1-σ)*ξ2^(σ-1))
Q_2(P1::Float64, P2::Float64, ξ1::Float64, ξ2::Float64, σ::Float64) = P2^(-σ)*ξ2^(σ-1)/(P1^(1-σ)*ξ1^(σ-1) + P2^(1-σ)*ξ2^(σ-1))

RSP_1(P1::Float64, P2::Float64, ξ1::Float64, ξ2::Float64, σ::Float64) = Q_1(P1, P2, ξ1, ξ2, σ)*P1/(Q_1(P1, P2, ξ1, ξ2, σ)*P1 + Q_2(P1, P2, ξ1, ξ2, σ)*P2)
RSP_2(P1::Float64, P2::Float64, ξ1::Float64, ξ2::Float64, σ::Float64) = Q_2(P1, P2, ξ1, ξ2, σ)*P2/(Q_1(P1, P2, ξ1, ξ2, σ)*P1 + Q_2(P1, P2, ξ1, ξ2, σ)*P2)



function Cournot_dπ1_dq1(q1::Float64, Q2::Float64, ξ1::Float64, ξ2::Float64, MC1::Float64, s1::Float64, σ::Float64)
    # (σ-1)/σ*P[c,i,m]*(1-RS[c,i,m]) - MC[c,i,m] 
    Q1 = exp(q1)
    return (1+s1)*(σ-1)/σ*P_1(Q1, Q2, ξ1, ξ2, σ)*(1-RS_1(Q1, Q2, ξ1, ξ2, σ)) - MC1
end

function Cournot_dπ2_dq2(Q1::Float64, q2::Float64, ξ1::Float64, ξ2::Float64, MC2::Float64, s2::Float64, σ::Float64)
    # (σ-1)/σ*P[c,i,m]*(1-RS[c,i,m]) - MC[c,i,m] 
    Q2 = exp(q2)
    return (1+s2)*(σ-1)/σ*P_2(Q1, Q2, ξ1, ξ2, σ)*(1-RS_2(Q1, Q2, ξ1, ξ2, σ)) - MC2
end

# R1(Q_2)
function R_f1(Q2::Float64, ξ1::Float64, ξ2::Float64, MC1::Float64, s1::Float64, σ::Float64)
    # find the root of Cournot_dπ1_dq1 to get Q_1
    return exp(find_zero(Q_1 -> Cournot_dπ1_dq1(Q_1, Q2, ξ1, ξ2, MC1, s1, σ), Q2))
end

# R2(Q_1)
function R_f2(Q1::Float64, ξ1::Float64, ξ2::Float64, MC2::Float64, s2::Float64, σ::Float64)
    # find the root of Cournot_dπ2_dq2 to get Q_2
    return exp(find_zero(Q2 -> Cournot_dπ2_dq2(Q1, Q2, ξ1, ξ2, MC2, s2, σ), Q1))
end

# get a series of R1(Q_2) for a series of Q_2
function R1_series(Q_2_series, ξ1::Float64, ξ2::Float64, MC1::Float64, s1::Float64, σ::Float64)
    return [R_f1(Q2, ξ1, ξ2, MC1, s1, σ) for Q2 in Q_2_series]
end

# get a series of R2(Q_1) for a series of Q_1
function R2_series(Q_1_series, ξ1::Float64, ξ2::Float64, MC2::Float64, s2::Float64, σ::Float64)
    return [R_f2(Q1, ξ1, ξ2, MC2, s2, σ) for Q1 in Q_1_series]
end

function Cournot_solve_eqba_q(ξ1::Float64, ξ2::Float64, MC1::Float64, MC2::Float64, s1::Float64, s2::Float64, σ::Float64)
    # find the equilibrium quantities Q1 and Q2
    # initial guess
    Q1 = 0.5
    Q2 = 0.5
    # find the equilibrium quantities, where the Cournot_dπ1_dq1 and Cournot_dπ2_dq2 are both zero
    while abs(Cournot_dπ1_dq1(log(Q1), Q2, ξ1, ξ2, MC1, s1, σ)) > 1e-6 || abs(Cournot_dπ2_dq2(Q1, log(Q2), ξ1, ξ2, MC2, s2, σ)) > 1e-6
        # update Q1 and Q2
        Q1 = exp(find_zero(q1 -> Cournot_dπ1_dq1(q1, Q2, ξ1, ξ2, MC1, s1, σ), log(Q1)))
        Q2 = exp(find_zero(q2 -> Cournot_dπ2_dq2(Q1, q2, ξ1, ξ2, MC2, s2, σ), log(Q2)))
    end
    return Q1, Q2
end


################ Bertrand functions ################

function Bertrand_dπ1_dp1(p1::Float64, p2::Float64, ξ1::Float64, ξ2::Float64, MC1::Float64, s1::Float64, σ::Float64)  # partial term, sufficient for determining the solution
    # Q_m[c,i]*(1-(P_m[c,i]-MC_m[c,i])/P_m[c,i]*(σ+(1-σ)*RS_m[c,i])) 
    return (1-((1+s1)*p1-MC1)/((1+s1)*p1)*(σ+(1-σ)*RSP_1(p1, p2, ξ1, ξ2, σ)))
end

function Bertrand_dπ2_dp2(p1::Float64, p2::Float64, ξ1::Float64, ξ2::Float64, MC2::Float64, s2::Float64, σ::Float64)
    # Q_m[c,i]/P_m[c,i]*(1-(P_m[c,i]-MC_m[c,i])/P_m[c,i]*(σ+(1-σ)*RS_m[c,i])) 
    return (1-((1+s2)*p2-MC2)/((1+s2)*p2)*(σ+(1-σ)*RSP_2(p1, p2, ξ1, ξ2, σ)))
end

# R1(P_2)
function R_f1_p2(P2::Float64, ξ1::Float64, ξ2::Float64, MC1::Float64, s1::Float64, σ::Float64)
    # find the root of Bertrand_dπ1_dp1 to get P_1
    return find_zero(P1 -> Bertrand_dπ1_dp1(P1, P2, ξ1, ξ2, MC1, s1, σ), P2)
end

# R2(P_1)
function R_f2_p1(P1::Float64, ξ1::Float64, ξ2::Float64, MC2::Float64, s2::Float64, σ::Float64)
    # find the root of Bertrand_dπ2_dp2 to get P_2
    return find_zero(P2 -> Bertrand_dπ2_dp2(P1, P2, ξ1, ξ2, MC2, s2, σ), P1)
end

# get a series of R1(P_2) for a series of P_2
function R1_series_p2(P_2_series, ξ1::Float64, ξ2::Float64, MC1::Float64, s1::Float64, σ::Float64)
    return [R_f1_p2(P2, ξ1, ξ2, MC1, s1, σ) for P2 in P_2_series]
end

# get a series of R2(P_1) for a series of P_1
function R2_series_p1(P_1_series, ξ1::Float64, ξ2::Float64, MC2::Float64, s2::Float64, σ::Float64)
    return [R_f2_p1(P1, ξ1, ξ2, MC2, s2, σ) for P1 in P_1_series]
end


function Bertrand_solve_eqba_p(ξ1::Float64, ξ2::Float64, MC1::Float64, MC2::Float64, s1::Float64, s2::Float64, σ::Float64)
    # find the equilibrium prices P1 and P2
    # initial guess
    P1 = 1.0
    P2 = 1.0
    # find the equilibrium prices, where the Bertrand_dπ1_dp1 and Bertrand_dπ2_dp2 are both zero
    while abs(Bertrand_dπ1_dp1(P1, P2, ξ1, ξ2, MC1, s1, σ)) > 1e-6 || abs(Bertrand_dπ2_dp2(P1, P2, ξ1, ξ2, MC2, s2, σ)) > 1e-6
        # update P1 and P2
        P1 = find_zero(p1 -> Bertrand_dπ1_dp1(p1, P2, ξ1, ξ2, MC1, s1, σ), P1)
        P2 = find_zero(p2 -> Bertrand_dπ2_dp2(P1, p2, ξ1, ξ2, MC2, s2, σ), P2)
    end
    return P1, P2
end