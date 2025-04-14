function Cournot_dπ_dsdq_m_f(eq::eqbm_m, c::Int)
    C = eq.IDT.C
    N = eq.IDT.N
    # initialize an empty vector with size C*N+1
    dπ_dsdq = Array{Float64}(undef, C*N + 1)
    for j in 1:C
        for i in 1:N
            #dπ_dsdq[(j-1)*N+i] = if(S_m[c]>0) 1.0 else 0.0 end  # wrong
            if j == c
                dπ_dsdq[(j-1)*N+i] = Cournot_d2π_icm_dq_icm_ds_cm(eq, c, i)
            else
                dπ_dsdq[(j-1)*N+i] = 0.0
            end
        end
    end
    dπ_dsdq[end] = 0.0
    return dπ_dsdq
end


function Cournot_d2π_dq1q2_mat_m_f(eq::eqbm_m)
    C = eq.IDT.C
    N = eq.IDT.N
    
    d2π_dq1q2_mat = Array{Float64}(undef, C*N+1, C*N+1)
    for c in 1:C
        for i in 1:N
            for k in 1:C
                for j in 1:N
                    if c == k && i == j
                        d2π_dq1q2_mat[(c-1)*N+i, (k-1)*N+j] = Cournot_d2π_icm_dq2_icm_single_f(eq, c, i)
                    else
                        d2π_dq1q2_mat[(c-1)*N+i, (k-1)*N+j] = Cournot_d2π_icm_dq_icm_dq_jkm_single_f(eq, c, i, k, j)
                    end
                end
            end
        end
    end
    # add last row and last column, except the last element
    for c in 1:C
        for i in 1:N
            # last column
            d2π_dq1q2_mat[(c-1)*N+i, C*N+1] = Cournot_d2π_icm_dq_icm_dq_0mm_single_f(eq, c, i)
            # last row
            d2π_dq1q2_mat[C*N+1, (c-1)*N+i] = Cournot_d2π_0mm_dq_0mm_dq_jkm_single_f(eq, c, i)
        end
    end
    d2π_dq1q2_mat[C*N+1, C*N+1] = Cournot_d2π0_dq02_m_single_f(eq)
    
    return d2π_dq1q2_mat
end


function Cournot_dq_jkm_ds_cm_f(eq::eqbm_m, j::Int, k::Int, M_mat::Array{Float64, 2}, det_M_mat::Float64, dπ_dsdq_vec::Array{Float64, 1})
    N = eq.IDT.N
    # M_mat = Cournot_d2π_dq1q2_mat_m_f(eq)
    # det_M_mat = det(M_mat)
    # dπ_dsdq_vec = Cournot_dπ_dsdq_m_f(eq)

    # replace the (k-1)*N+j column with -dπ_dsdq_vec
    M_jk_mat = deepcopy(M_mat)
    M_jk_mat[:, (k-1)*N+j] = -dπ_dsdq_vec
    det_M_jk_mat = det(M_jk_mat)
    
    return det_M_jk_mat/det_M_mat
end

function Cournot_dq_0mm_ds_cm_f(eq::eqbm_m, M_mat::Array{Float64, 2}, det_M_mat::Float64, dπ_dsdq_vec::Array{Float64, 1})
    N = eq.IDT.N
    C = eq.IDT.C
    # M_mat = Cournot_d2π_dq1q2_mat_m_f(eq)
    # det_M_mat = det(M_mat)
    # dπ_dsdq_vec = Cournot_dπ_dsdq_m_f(eq)

    # replace the last column with -dπ_dsdq_vec
    M_00_mat = deepcopy(M_mat)
    M_00_mat[:, C*N+1] = -dπ_dsdq_vec
    det_M_00_mat = det(M_00_mat)
    
    return det_M_00_mat/det_M_mat
end



function Cournot_optimal_s_m_c_f(IDT::input_dt_m, c::Int) # return a scalar, optimal s for production origin c
    IDT_m = deepcopy(IDT)
    N = IDT.N
    C = IDT.C
   # C_list = [i for i in 1:C if i != c] 
    TOL = IDT.TOL
    # MAXIT = IDT.MAXIT
    MAXIT = 5000
    k = 1
    diff = 1.0

    S_old = 0.0
    S_temp = deepcopy(IDT.S_m)
    
    # initiate an eqbm_m object
    # eq = DGP_m("C", IDT)
    eq = eqbm_m(IDT, "C", zeros(C,N), 0.0, zeros(C,N), 0.0, zeros(C,N), 0.0, zeros(C,N), 0.0, 0.0, zeros(C,N), 0.0)

    while diff > TOL && k<MAXIT
        S_temp[c] = S_old
        IDT_m.S_m .= S_temp
        try 
          eq = DGP_m("C", IDT_m)
        catch
            return NaN
        end 
        P_m = eq.P_m
        Q_m = eq.Q_m
        M_mat = Cournot_d2π_dq1q2_mat_m_f(eq)
        det_M_mat = det(M_mat)
        
        dπ_dsdq_vec = Cournot_dπ_dsdq_m_f(eq, c) #  C*N+1 vector, including OG
        dπ_icm_dq_0mm = [Cournot_dπ_icm_dq_0mm_single_f(eq, c, i) for i in 1:N]
        # dπ_icm_dq_jkm_3 = [Cournot_dπ_icm_dq_jkm_single_f(eq, c, i, k, j) for i in 1:N, k in C_list, j in 1:N]
        dπ_icm_dq_jkm_3 = [Cournot_dπ_icm_dq_jkm_single_f(eq, c, i, k, j) for i in 1:N, k in 1:C, j in 1:N] # 3d array
        dq_0mm_ds_cm = Cournot_dq_0mm_ds_cm_f(eq, M_mat, det_M_mat, dπ_dsdq_vec) # scalar
        # dq_jkm_ds_cm_2 = [Cournot_dq_jkm_ds_cm_f(eq, j, k, M_mat, det_M_mat, dπ_dsdq_vec) for k in C_list, j in 1:N]
        dq_jkm_ds_cm_2 = [Cournot_dq_jkm_ds_cm_f(eq, j, k, M_mat, det_M_mat, dπ_dsdq_vec) for k in 1:C, j in 1:N] # matrix
        dq_icm_ds_cm_1 = [Cournot_dq_jkm_ds_cm_f(eq, i, c, M_mat, det_M_mat, dπ_dsdq_vec) for i in 1:N] # vector
        
        # add 
        dp_icm_dq_jkm_3 = [Cournot_dp_icm_dq_jkm_single_f(eq, c, i, k, j) for i in 1:N, k in 1:C, j in 1:N] # 3d array
        dp_icm_dq_icm_1 = [Cournot_dp_icm_dq_icm_single_f(eq, c, i) for i in 1:N] # vector
        dp_icm_dq_0mm = [Cournot_dp_icm_dq_0mm_single_f(eq, c, i) for i in 1:N] # vector

        S_new = (
            sum([dπ_icm_dq_0mm[i] * dq_0mm_ds_cm for i in 1:N]) +
            sum([dπ_icm_dq_jkm_3[i, k, j] * dq_jkm_ds_cm_2[k, j] for i in 1:N, k in 1:C, j in 1:N]) -
            sum([dπ_icm_dq_jkm_3[i, c, i] * dq_jkm_ds_cm_2[c, i] for i in 1:N])
        ) / (
            sum([dq_icm_ds_cm_1[i] * P_m[c, i] for i in 1:N]) +
            sum([
                Q_m[c, index] * (
                    sum([dp_icm_dq_jkm_3[index, k, j] * dq_jkm_ds_cm_2[k, j] for k in 1:C, j in 1:N]) -
                    dp_icm_dq_jkm_3[index, c, index] * dq_jkm_ds_cm_2[c, index] +
                    dp_icm_dq_icm_1[index] * dq_icm_ds_cm_1[index]  +
                   dp_icm_dq_0mm[index] * dq_0mm_ds_cm 
                ) for index in 1:N
            ])
        )

        diff = abs(S_new - S_old)
        S_old = w*S_old + (1-w)*S_new
        k += 1
    end
   #  @assert k < MAXIT "No convergence of S after $(MAXIT) iterations"
    if k >= MAXIT
        # println("No convergence of S after $(MAXIT) iterations")
        return NaN
    end

    return S_old
end



# for Cournot check purposes, no subsidy/tax
QS_m_vec(Q_m, exp_ξ_m, σ) =  Σ((Q_m.^((σ-1)/σ)).*(exp_ξ_m).^(1/σ))  # for ForwardDiff
π_icm_q_vec(i, Q_m, MC_m, exp_ξ_m, Y_m, Q0_m, exp_ξ0_m, σ) = (Y_m*Q_m[i]^(-1/σ)*(exp_ξ_m[i])^(1/σ)/(Q0_m^((σ-1)/σ)*exp_ξ0_m^(1/σ) + QS_m_vec(Q_m,exp_ξ_m,σ)) - MC_m[i])*Q_m[i]


function π_q_jacobian(Q_m, MC_m, exp_ξ_m, Y_m, Q0_m, exp_ξ0_m, σ)
    # Flatten the input matrices to vectors
    q_vec = vec(Q_m)
    exp_ξ_vec = vec(exp_ξ_m)
    mc_vec = vec(MC_m)
    jacobian_matrix = zeros(length(q_vec), length(q_vec))
    for i in 1:length(q_vec)
        jacobian_matrix[i,:] =  ForwardDiff.gradient(q -> π_icm_q_vec(i, q, mc_vec, exp_ξ_vec, Y_m, Q0_m, exp_ξ0_m, σ), q_vec)
    end

    return jacobian_matrix
end


function π_q_hessian(Q_m, MC_m, exp_ξ_m, Q0_m, exp_ξ0_m, Y_m, σ::Real)
    # Flatten the input matrices to vectors
    q_vec = vec(Q_m)
    exp_ξ_vec = vec(exp_ξ_m)
    mc_vec = vec(MC_m)

    # Wrapper function that computes the combined values as a scalar sum
    # Compute the Hessian matrix using ForwardDiff
    hessian_matrix = zeros(size(Q_m, 1),size(Q_m, 1), length(q_vec))
    for i in 1:length(q_vec)
        hessian_matrix[:,:, i] = ForwardDiff.hessian(q -> π_icm_q_vec(i, q, mc_vec, exp_ξ_vec, Y_m, Q0_m, exp_ξ0_m, σ), q_vec)
    end
    return hessian_matrix
end


# # function for C=2 N=1 only (haven't modified yet) ===============
# function Cournot_optimal_s1_m_f(IDT::input_dt_m)
#     """
#     Optimal tax/subsidy for firm 1. Solve by fixed point iteration?
#     """
#     IDT_m = deepcopy(IDT)
#     TOL = IDT_m.TOL
#     MAXIT = IDT_m.MAXIT
#     k = 1
#     diff = 1.0

#     S_old = [0.0, 0.0]
#     while diff > TOL && k<MAXIT
#         IDT_m.S_m .= S_old
#         eq = DGP_m("C", IDT_m)
#         S_new = [-(Cournot_dπ_icm_dq_jkm_single_f(eq, 1, 1, 2, 1)*(Cournot_d2π_icm_dq_icm_dq_jkm_single_f(eq, 2, 1, 1, 1)*Cournot_d2π0_dq02_m_single_f(eq)-Cournot_d2π_0mm_dq_0mm_dq_jkm_single_f(eq, 1, 1)*Cournot_d2π_icm_dq_icm_dq_0mm_single_f(eq,2,1))
#         + Cournot_dπ_icm_dq_0mm_single_f(eq,1,1)*(Cournot_d2π_0mm_dq_0mm_dq_jkm_single_f(eq, 1, 1)*Cournot_d2π_icm_dq2_icm_single_f(eq,2,1)-Cournot_d2π_icm_dq_icm_dq_jkm_single_f(eq, 2, 1, 1, 1)*Cournot_d2π_0mm_dq_0mm_dq_jkm_single_f(eq, 2, 1)))/(Cournot_d2π_icm_dq2_icm_single_f(eq,2,1)*Cournot_d2π0_dq02_m_single_f(eq)-Cournot_d2π_0mm_dq_0mm_dq_jkm_single_f(eq, 2, 1)*Cournot_d2π_icm_dq_icm_dq_0mm_single_f(eq,2,1)), 0.0]
#         diff = norm(S_new - S_old)
#         S_old = w.*S_old + (1-w).*S_new
#         k += 1
#     end

#     # @assert k < MAXIT "No convergence of S after $(MAXIT) iterations"
#     if k >= MAXIT
#         return NaN
#     end

#     return S_old    
# end


# function Cournot_optimal_s2_m_f(IDT::input_dt_m)
#     """
#     Optimal tax/subsidy for firm 1. Solve by fixed point iteration?
#     """
#     IDT_m = deepcopy(IDT)
#     TOL = IDT_m.TOL
#     MAXIT = IDT_m.MAXIT
#     k = 1
#     diff = 1.0

#     S_old = [0.0, 0.0]
#     while diff > TOL && k<MAXIT
#         IDT_m.S_m .= S_old
#         eq = DGP_m("C", IDT_m)
#         S_new = [0.0, -(Cournot_dπ_icm_dq_jkm_single_f(eq, 2, 1, 1, 1)*(Cournot_d2π_icm_dq_icm_dq_jkm_single_f(eq, 1, 1, 2, 1)*Cournot_d2π0_dq02_m_single_f(eq)-Cournot_d2π_0mm_dq_0mm_dq_jkm_single_f(eq, 2, 1)*Cournot_d2π_icm_dq_icm_dq_0mm_single_f(eq,1,1))
#         + Cournot_dπ_icm_dq_0mm_single_f(eq,2,1)*(Cournot_d2π_0mm_dq_0mm_dq_jkm_single_f(eq, 2, 1)*Cournot_d2π_icm_dq2_icm_single_f(eq,1,1)-Cournot_d2π_icm_dq_icm_dq_jkm_single_f(eq, 1, 1, 2, 1)*Cournot_d2π_0mm_dq_0mm_dq_jkm_single_f(eq, 1, 1)))/(Cournot_d2π_icm_dq2_icm_single_f(eq,1,1)*Cournot_d2π0_dq02_m_single_f(eq)-Cournot_d2π_0mm_dq_0mm_dq_jkm_single_f(eq, 1, 1)*Cournot_d2π_icm_dq_icm_dq_0mm_single_f(eq,1,1))]
#         diff = norm(S_new - S_old)
#         S_old = w.*S_old + (1-w).*S_new
#         k += 1
#     end

#     # @assert k < MAXIT "No convergence of S after $(MAXIT) iterations"
#     if k >= MAXIT
#         return NaN
#     end
#     return S_old    
# end


function Bertrand_dπ_dsdp_m_f(eq::eqbm_m, c::Int)
    C = eq.IDT.C
    N = eq.IDT.N
    # initialize an empty vector with size C*N+1
    dπ_dsdp = Array{Float64}(undef, C*N + 1)
    for j in 1:C
        for i in 1:N
            #dπ_dsdq[(j-1)*N+i] = if(S_m[c]>0) 1.0 else 0.0 end  # wrong
            if j == c
                dπ_dsdp[(j-1)*N+i] = Bertrand_d2π_icm_dp_icm_ds_cm(eq, c, i)
            else
                dπ_dsdp[(j-1)*N+i] = 0.0
            end
        end
    end
    dπ_dsdp[end] = 0.0
    return dπ_dsdp
end




function Bertrand_d2π_dp1p2_mat_m_f(eq::eqbm_m)
    C = eq.IDT.C
    N = eq.IDT.N
    
    d2π_dp1p2_mat = Array{Float64}(undef, C*N+1, C*N+1)
    for c in 1:C
        for i in 1:N
            for k in 1:C
                for j in 1:N
                    if c == k && i == j
                        d2π_dp1p2_mat[(c-1)*N+i, (k-1)*N+j] = Bertrand_d2π_icm_dp2_icm_single_f(eq, c, i)
                    else
                        d2π_dp1p2_mat[(c-1)*N+i, (k-1)*N+j] = Bertrand_d2π_icm_dp_icm_dp_jkm_single_f(eq, c, i, k, j)
                    end
                end
            end
        end
    end
    # add last row and last column, except the last element
    for c in 1:C
        for i in 1:N
            # last column
            d2π_dp1p2_mat[(c-1)*N+i, C*N+1] = Bertrand_d2π_icm_dp_icm_dp_0mm_single_f(eq, c, i)
            # last row
            d2π_dp1p2_mat[C*N+1, (c-1)*N+i] = Bertrand_d2π_0mm_dp_0mm_dp_jkm_single_f(eq, c, i)
        end
    end
    d2π_dp1p2_mat[C*N+1, C*N+1] = Bertrand_d2π0_dp02_m_single_f(eq)
    
    return d2π_dp1p2_mat
end




function Bertrand_dp_jkm_ds_cm_f(eq::eqbm_m, j::Int, k::Int, M_mat::Array{Float64, 2}, det_M_mat::Float64, dπ_dsdp_vec::Array{Float64, 1})
    N = eq.IDT.N

    # replace the (k-1)*N+j column with -dπ_dsdp_vec
    M_jk_mat = deepcopy(M_mat)
    M_jk_mat[:, (k-1)*N+j] = -dπ_dsdp_vec
    det_M_jk_mat = det(M_jk_mat)
    
    return det_M_jk_mat/det_M_mat
end



function Bertrand_dp_0mm_ds_cm_f(eq::eqbm_m, M_mat::Array{Float64, 2}, det_M_mat::Float64, dπ_dsdp_vec::Array{Float64, 1})
    N = eq.IDT.N
    C = eq.IDT.C
  
    # replace the last column with -dπ_dsdp_vec
    M_00_mat = deepcopy(M_mat)
    M_00_mat[:, C*N+1] = -dπ_dsdp_vec
    det_M_00_mat = det(M_00_mat)
    
    return det_M_00_mat/det_M_mat
end




function Bertrand_optimal_s_m_c_f(IDT::input_dt_m, c::Int) # return a scalar, optimal s for production origin c
    IDT_m = deepcopy(IDT)
    N = IDT.N
    C = IDT.C
   # C_list = [i for i in 1:C if i != c] 
    TOL = IDT.TOL
    # MAXIT = IDT.MAXIT
    MAXIT = 5000
    k = 1
    diff = 1.0

    S_old = 0.0
    S_temp = deepcopy(IDT.S_m)

    # initiate an eqbm_m object
    # eq = DGP_m("B", IDT) 
    eq = eqbm_m(IDT, "B", zeros(C,N), 0.0, zeros(C,N), 0.0, zeros(C,N), 0.0, zeros(C,N), 0.0, 0.0, zeros(C,N), 0.0)

    while diff > TOL && k<MAXIT
        S_temp[c] = S_old
        IDT_m.S_m .= S_temp
        try 
            eq = DGP_m("B", IDT_m)
        catch
            return NaN
        end 
        P_m = eq.P_m
        Q_m = eq.Q_m
        M_mat = Bertrand_d2π_dp1p2_mat_m_f(eq)
        det_M_mat = det(M_mat)

        dπ_dsdp_vec = Bertrand_dπ_dsdp_m_f(eq, c) #  C*N+1 vector, including OG
        dπ_icm_dp_0mm = [Bertrand_dπ_icm_dp_0mm_single_f(eq, c, i) for i in 1:N]
        dπ_icm_dp_jkm_3 = [Bertrand_dπ_icm_dp_jkm_single_f(eq, c, i, k, j) for i in 1:N, k in 1:C, j in 1:N] # 3d array

        dp_0mm_ds_cm = Bertrand_dp_0mm_ds_cm_f(eq, M_mat, det_M_mat, dπ_dsdp_vec) # scalar

        dp_jkm_ds_cm_2 = [Bertrand_dp_jkm_ds_cm_f(eq, j, k, M_mat, det_M_mat, dπ_dsdp_vec) for k in 1:C, j in 1:N] # matrix

        dp_icm_ds_cm_1 = [Bertrand_dp_jkm_ds_cm_f(eq, i, c, M_mat, det_M_mat, dπ_dsdp_vec) for i in 1:N] # vector

        dq_icm_dp_0mm = [Bertrand_dq_icm_dp_0mm_single_f(eq, c, i) for i in 1:N] # vector

        dq_icm_dp_jkm_3 = [Bertrand_dq_icm_dp_jkm_single_f(eq, c, i, k, j) for i in 1:N, k in 1:C, j in 1:N] # 3d array

        dq_icm_dp_icm_1 = [Bertrand_dq_icm_dp_icm_single_f(eq, c, i) for i in 1:N] # vector

        # S_new = (sum([dπ_icm_dp_0mm[i]*dp_0mm_ds_cm for i in 1:N]) 
        # + sum([dπ_icm_dp_jkm_3[i,k,j]*dp_jkm_ds_cm_2[k,j] for i in 1:N, k in 1:C, j in 1:N]) -sum([dπ_icm_dp_jkm_3[i,c,i]*dp_jkm_ds_cm_2[c,i] for i in 1:N]))/ (sum([dp_icm_ds_cm_1[i]*Q_m[c,i] for i in 1:N]) + sum([P_m[c,i]* (sum([dq_icm_dp_jkm_3[i,k,j]*dp_jkm_ds_cm_2[k,j] for i in 1:N, k in 1:C, j in 1:N]) 
        # - sum([dq_icm_dp_jkm_3[i,c,i]*dp_jkm_ds_cm_2[c,i] for i in 1:N])
        # + sum([dq_icm_dp_icm_1[i]*dp_icm_ds_cm_1[i] for i in 1:N])
        # + sum([dq_icm_dp_0mm[i]*dp_0mm_ds_cm for i in 1:N]))] for i in 1:N)
        # )
        S_new = (
    sum([dπ_icm_dp_0mm[i] * dp_0mm_ds_cm for i in 1:N]) +
    sum([dπ_icm_dp_jkm_3[i, k, j] * dp_jkm_ds_cm_2[k, j] for i in 1:N, k in 1:C, j in 1:N]) -
    sum([dπ_icm_dp_jkm_3[i, c, i] * dp_jkm_ds_cm_2[c, i] for i in 1:N])
) / (
    sum([dp_icm_ds_cm_1[i] * Q_m[c, i] for i in 1:N]) +
    sum([
        P_m[c, index] * (
            sum([dq_icm_dp_jkm_3[index, k, j] * dp_jkm_ds_cm_2[k, j] for k in 1:C, j in 1:N]) -
            dq_icm_dp_jkm_3[index, c, index] * dp_jkm_ds_cm_2[c, index]  +
            dq_icm_dp_icm_1[index] * dp_icm_ds_cm_1[index]+
            dq_icm_dp_0mm[index] * dp_0mm_ds_cm 
        ) for index in 1:N
    ])
)

        diff = abs(S_new - S_old)
        S_old = w*S_old + (1-w)*S_new
        k += 1
    end
    # @assert k < MAXIT "No convergence of S after $(MAXIT) iterations"
    if k >= MAXIT
        return NaN
    end

    return S_old
end



# for Bertrand check purposes  =======

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