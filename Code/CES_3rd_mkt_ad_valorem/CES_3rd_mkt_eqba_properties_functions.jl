# Cournot FOC
function Cournot_dπ_dq_f(eq::Union{eqbm_no_OG, eqbm_full_no_OG})
    C = eq.IDT.C
    N = eq.IDT.N
    M = eq.IDT.M
    RS = eq.RS
    σ = eq.IDT.σ
    P = eq.P
    MC = eq.MC
    S = eq.S

    FOC_vec = zeros(C, N, C+M)
    for c in 1:C
        for i in 1:N
            for m in 1:C+M
                FOC_vec[c,i,m] = (1+S[c,i,m])*(σ-1)/σ*P[c,i,m]*(1-RS[c,i,m]) - MC[c,i,m] 
            end
        end
    end
    
    return FOC_vec
end

# Cournot FOC, market m 
function Cournot_dπ_dq_m_f(eq::Union{eqbm_no_OG, eqbm_full_no_OG}, m::Int)
    C = eq.IDT.C
    N = eq.IDT.N
    RS = eq.RS
    σ = eq.IDT.σ
    P = eq.P
    MC = eq.MC
    S = eq.S

    FOC_mat = zeros(C, N)
    for c in 1:C
        for i in 1:N
                FOC_mat[c,i] = (1+S[c,i,m])*(σ-1)/σ*P[c,i,m]*(1-RS[c,i,m]) - MC[c,i,m] 
        end
    end 
    return FOC_mat
end


function Cournot_dπ_dq_m_single_f(eq::eqbm_m_no_OG)
    C = eq.IDT.C
    N = eq.IDT.N
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    MC_m = eq.MC_m
    S_m = eq.S_m

    FOC_mat = zeros(C, N)
    for c in 1:C
        for i in 1:N
                FOC_mat[c,i] = (1+S_m[c,i])*(σ-1)/σ*P_m[c,i]*(1-RS_m[c,i]) - MC_m[c,i] 
        end
    end 
    return FOC_mat
end


# Cournot SOC
function Cournot_d2π_dq2_f(eq::Union{eqbm_no_OG, eqbm_full_no_OG})
    C = eq.IDT.C
    N = eq.IDT.N
    M = eq.IDT.M
    RS = eq.RS
    σ = eq.IDT.σ
    # P = eq.P
    Y = eq.IDT.Y
    Q = eq.Q
    S = eq.S

    SOC_vec = zeros(C, N, C+M)
    for c in 1:C
        for i in 1:N
            for m in 1:C+M
                # q_icm = RS[c,i,m]*Y[m]/P[c,i,m] # quantity
                SOC_vec[c,i,m] = (1+S[c,i,m])*Y[m]/Q[c,i,m]^2*((1-σ)/σ^2*RS[c,i,m] + 
                    (3-2σ)*(σ-1)/σ^2*RS[c,i,m]^2+ 2*(σ-1)^2/σ^2*RS[c,i,m]^3)
            end
        end
    end
    
    return SOC_vec
end

function Cournot_d2π_dq2_single_f(eq::eqbm_m_no_OG)
    C = eq.IDT.C
    N = eq.IDT.N
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    # P = eq.P
    Y_m = eq.IDT.Y_m
    Q_m = eq.Q_m
    S_m = eq.S_m

    SOC_vec = zeros(C, N)
    for c in 1:C
        for i in 1:N
            # q_icm = RS[c,i,m]*Y[m]/P[c,i,m] # quantity
            SOC_vec[c,i] = (1+S_m[c,i])*Y_m/Q_m[c,i]^2*((1-σ)/σ^2*RS_m[c,i] + (3-2σ)*(σ-1)/σ^2*RS_m[c,i]^2+ 2*(σ-1)^2/σ^2*RS_m[c,i]^3)
        end
    end
    
    return SOC_vec
end


# Cournot SOC, market m
function Cournot_d2π_dq2_m_f(eq::Union{eqbm_no_OG, eqbm_full_no_OG}, m::Int)
    C = eq.IDT.C
    N = eq.IDT.N
    RS = eq.RS
    σ = eq.IDT.σ
    # P = eq.P
    Y = eq.IDT.Y
    Q = eq.Q
    S = eq.S

    SOC_mat = zeros(C, N)
    for c in 1:C
        for i in 1:N
                # q_icm = RS[c,i,m]*Y[m]/P[c,i,m] # quantity
                SOC_mat[c,i] = (1+S[c,i,m])*Y[m]/Q[c,i,m]^2*((1-σ)/σ^2*RS[c,i,m] + 
                    (3-2σ)*(σ-1)/σ^2*RS[c,i,m]^2+ 2*(σ-1)^2/σ^2*RS[c,i,m]^3)
        end
    end
    return SOC_mat
end


function Cournot_d2π_dq2_m_single_f(eq::eqbm_m_no_OG)
    C = eq.IDT.C
    N = eq.IDT.N
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    # P = eq.P
    Y_m = eq.IDT.Y_m
    Q_m = eq.Q_m
    S_m = eq.S_m

    SOC_mat = zeros(C, N)
    for c in 1:C
        for i in 1:N
                # q_icm = RS[c,i,m]*Y[m]/P[c,i,m] # quantity
                SOC_mat[c,i] = (1+S_m[c,i])*Y_m/Q_m[c,i]^2*((1-σ)/σ^2*RS_m[c,i] + 
                    (3-2σ)*(σ-1)/σ^2*RS_m[c,i]^2+ 2*(σ-1)^2/σ^2*RS_m[c,i]^3)
        end
    end
    return SOC_mat
end


function Cournot_dπ_icm_dq_jkm_single_f(eq::eqbm_m_no_OG, c::Int, i::Int, k::Int, j::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    S_m = eq.S_m

    return -(1+S_m[c,i])*(σ-1)/σ*RS_m[c,i]*P_m[k,j]
end




# Cournot strategic complementarity or substitutes, i and j from different production origins
function Cournot_d2π_icm_dq_icm_dq_jkm_f(eq::Union{eqbm_no_OG, eqbm_full_no_OG}, m::Int, c::Int, i::Int, k::Int, j::Int)
    RS = eq.RS
    σ = eq.IDT.σ
    P = eq.P
    Y = eq.IDT.Y
    S = eq.S

    return -(1+S[c,i,m])/Y[m]*(σ-1)^2/σ^2*P[c,i,m]*P[k,j,m]*(1-2*RS[c,i,m])
end


function Cournot_d2π_icm_dq_icm_dq_jkm_single_f(eq::eqbm_m_no_OG, c::Int, i::Int, k::Int, j::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    S_m = eq.S_m
    return -(1+S_m[c,i])/Y_m*(σ-1)^2/σ^2*P_m[c,i]*P_m[k,j]*(1-2*RS_m[c,i])
end


function Cournot_dp_dq_single_f(eq::eqbm_m_no_OG)
    C = eq.IDT.C
    N = eq.IDT.N
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    Q_m = eq.Q_m
    Y_m = eq.IDT.Y_m
    dp_dq_vec = zeros(C, N)
    for c in 1:C
        for i in 1:N
            dp_dq_vec[c,i] = - 1/σ*Y_m/Q_m[c,i]^2*(RS_m[c,i] + (σ-1)*RS_m[c,i]^2)
        end
    end
    
    return  dp_dq_vec
end
 


function Cournot_dp_icm_dq_jkm_single_f(eq::eqbm_m_no_OG, c::Int, i::Int, k::Int, j::Int)
    #RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    
    
    # return (σ-1)*Y_m*RS_m[c,i]*RS_m[k,j]/(P_m[c,i]*P_m[k,j])
    return -(σ-1)/σ * P_m[c,i]*P_m[k,j]/Y_m
end



# Bertrand FOC
function Bertrand_dπ_dp_f(eq::Union{eqbm_no_OG, eqbm_full_no_OG})
    C = eq.IDT.C
    N = eq.IDT.N
    M = eq.IDT.M
    RS = eq.RS
    σ = eq.IDT.σ
    P = eq.P
    MC = eq.MC
    Y = eq.IDT.Y
    S = eq.S

    FOC_vec = zeros(C, N, C+M)
    for c in 1:C
        for i in 1:N
            for m in 1:C+M
                FOC_vec[c,i,m] = (1+S[c,i,m])*Y[m]*RS[c,i,m]/P[c,i,m]*(1-((1+S[c,i,m])*P[c,i,m]-MC[c,i,m])/((1+S[c,i,m])*P[c,i,m]) * (σ + (1-σ)*RS[c,i,m]))
            end
        end
    end
    
    return FOC_vec
end



# Bertrand FOC, market m
function Bertrand_dπ_dp_m_f(eq::Union{eqbm_no_OG, eqbm_full_no_OG}, m::Int)
    C = eq.IDT.C
    N = eq.IDT.N
    RS = eq.RS
    σ = eq.IDT.σ
    P = eq.P
    MC = eq.MC
    Y = eq.IDT.Y
    S = eq.S

    FOC_mat = zeros(C, N)
    for c in 1:C
        for i in 1:N
                FOC_mat[c,i] = (1+S[c,i,m])*Y[m]*RS[c,i,m]/P[c,i,m]*(1-((1+S[c,i,m])*P[c,i,m]-MC[c,i,m])/((1+S[c,i,m])*P[c,i,m]) * (σ + (1-σ)*RS[c,i,m]))
        end
    end
    
    return FOC_mat
end


function Bertrand_dπ_dp_m_single_f(eq::eqbm_m_no_OG)
    C = eq.IDT.C
    N = eq.IDT.N
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    MC_m = eq.MC_m
    Y_m = eq.IDT.Y_m
    S_m = eq.S_m

    FOC_mat = zeros(C, N)
    for c in 1:C
        for i in 1:N
                FOC_mat[c,i] = (1+S_m[c,i])*Y_m*RS_m[c,i]/P_m[c,i]*(1-((1+S_m[c,i])*P_m[c,i]-MC_m[c,i])/((1+S_m[c,i])*P_m[c,i])*(σ+(1-σ)*RS_m[c,i]))
        end
    end
    
    return FOC_mat
end



# Bertrand SOC
function Bertrand_d2π_dp2_f(eq::Union{eqbm_no_OG, eqbm_full_no_OG})
    C = eq.IDT.C
    N = eq.IDT.N
    M = eq.IDT.M
    RS = eq.RS
    σ = eq.IDT.σ
    P = eq.P
    Y = eq.IDT.Y
    MC = eq.MC
    S = eq.S

    SOC_vec = zeros(C, N, C+M)
    for c in 1:C
        for i in 1:N
            for m in 1:C+M
                SOC_vec[c,i,m] = -2*(1+S[c,i,m])*Y[m]*RS[c,i,m]/P[c,i,m]^2*(σ+(1-σ)*RS[c,i,m]) + ((1+S[c,i,m])*P[c,i,m]-MC[c,i,m])*Y[m]/P[c,i,m]^3*(σ*(σ+1)*RS[c,i,m]+3*σ*(1-σ)*RS[c,i,m]^2+2*(1-σ)^2*RS[c,i,m]^3)
            end
        end
    end
    
    return SOC_vec
end


# Bertrand SOC, market m 
function Bertrand_d2π_dp2_m_f(eq::Union{eqbm_no_OG, eqbm_full_no_OG}, m::Int)
    C = eq.IDT.C
    N = eq.IDT.N
    M = eq.IDT.M
    RS = eq.RS
    σ = eq.IDT.σ
    P = eq.P
    Y = eq.IDT.Y
    MC = eq.MC
    S = eq.S

    SOC_mat = zeros(C, N)
    for c in 1:C
        for i in 1:N
                SOC_mat[c,i] = -2*(1+S[c,i,m])*Y[m]*RS[c,i,m]/P[c,i,m]^2*(σ+(1-σ)*RS[c,i,m]) + ((1+S[c,i,m])*P[c,i,m]-MC[c,i,m])*Y[m]/P[c,i,m]^3*(σ*(σ+1)*RS[c,i,m]+3*σ*(1-σ)*RS[c,i,m]^2+2*(1-σ)^2*RS[c,i,m]^3)
        end
    end
    
    return SOC_mat
end


# Bertrand SOC
function Bertrand_d2π_dp2_single_f(eq::eqbm_m_no_OG)
    C = eq.IDT.C
    N = eq.IDT.N
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    MC_m = eq.MC_m
    S_m = eq.S_m

    SOC_vec = zeros(C, N)
    for c in 1:C
        for i in 1:N
            SOC_vec[c,i] = -2*(1+S_m[c,i])*Y_m*RS_m[c,i]/P_m[c,i]^2*(σ+(1-σ)*RS_m[c,i]) + ((1+S_m[c,i])*P_m[c,i]-MC_m[c,i])*Y_m/P_m[c,i]^3*(σ*(σ+1)*RS_m[c,i]+3*σ*(1-σ)*RS_m[c,i]^2+2*(1-σ)^2*RS_m[c,i]^3)
        end
    end
    
    return SOC_vec
end


function Bertrand_dπ_icm_dp_jkm_single_f(eq::eqbm_m_no_OG, c::Int, i::Int, k::Int, j::Int)
    Y_m = eq.Y_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Q_m = eq.Q_m
    MC_m = eq.MC_m
    S_m = eq.S_m

    return (σ-1)*((1+S_m[c,i])*P_m[c,i]-MC_m[c,i])/Y_m*Q_m[c,i]*Q_m[k,j]
end


# Bertrand strategic complementarity or substitutes, i and j from different production origins
function Bertrand_d2π_icm_dp_icm_dp_jkm_f(eq::Union{eqbm_no_OG, eqbm_full_no_OG}, m::Int, c::Int, i::Int, k::Int, j::Int)
    RS = eq.RS
    σ = eq.IDT.σ
    P = eq.P
    Y = eq.IDT.Y
    MC = eq.MC
    S = eq.S

    return (1-σ)^2*Y[m]*((1+S[c,i,m])*P[c,i,m]-MC[c,i,m])/P[c,i,m]^2*RS[c,i,m]^2*RS[k,j,m]/P[k,j,m]
end



function Bertrand_d2π_icm_dp_icm_dp_jkm_single_f(eq::eqbm_m_no_OG, c::Int, i::Int, k::Int, j::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    MC_m = eq.MC_m
    S_m = eq.S_m

   return (1-σ)^2*Y_m*((1+S_m[c,i])*P_m[c,i]-MC_m[c,i])/P_m[c,i]^2*RS_m[c,i]^2*RS_m[k,j]/P_m[k,j]
   #return (σ-1)*Y_m/(P_m[c,i]*P_m[k,j])*RS_m[c,i]*RS_m[k,j] +(P_m[c,i]-MC_m[c,i])*(-σ*(σ-1)*Y_m*RS_m[c,i]*RS_m[k,j]/(P_m[c,i]^2*P_m[k,j])+ 2*(σ-1)^2*Y_m*RS_m[c,i]^2*RS_m[k,j]/(P_m[c,i]^2*P_m[k,j]))
end



function Bertrand_dq_dp_single_f(eq::eqbm_m_no_OG)
    C = eq.IDT.C
    N = eq.IDT.N
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    dq_dp_vec = zeros(C, N)
    for c in 1:C
        for i in 1:N
            dq_dp_vec[c,i] = -Y_m*RS_m[c,i]/P_m[c,i]^2*(σ+(1-σ)*RS_m[c,i])
        end
    end
    
    return  dq_dp_vec
end



function Bertrand_dq_icm_dp_jkm_single_f(eq::eqbm_m_no_OG, c::Int, i::Int, k::Int, j::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    
    
    return (σ-1)*Y_m*RS_m[c,i]*RS_m[k,j]/(P_m[c,i]*P_m[k,j])
end