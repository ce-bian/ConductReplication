# slow version... Ideally just input relevant parameters. Not sure if it affects the performance a lot... 
# Cournot FOC (inside goods), market m 
function Cournot_dπ_dq_m_single_f(eq::eqbm_m)
    C = eq.IDT.C
    N = eq.IDT.N
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    MC_m = eq.MC_m
    S_m = eq.IDT.S_m

    FOC_mat = zeros(C, N)
    for c in 1:C
        for i in 1:N
                FOC_mat[c,i] = (1+S_m[c])*(σ-1)/σ*P_m[c,i]*(1-RS_m[c,i]) - MC_m[c,i] 
        end
    end 
    return FOC_mat
end


# Cournot FOC (outside good), market m 
function Cournot_dπ0_dq0_m_single_f(eq::eqbm_m)
    RS0_m = eq.RS0_m
    σ = eq.IDT.σ
    P0_m = eq.P0_m
    MC0_m = eq.MC0_m

    return (σ-1)/σ*P0_m*(1-RS0_m) - MC0_m 
end

# Cournot SOC (inside goods), market m 
function Cournot_d2π_dq2_m_single_f(eq::eqbm_m)
    C = eq.IDT.C
    N = eq.IDT.N
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    # P = eq.P
    Y_m = eq.IDT.Y_m
    Q_m = eq.Q_m
    S_m = eq.IDT.S_m

    SOC_mat = zeros(C, N)
    for c in 1:C
        for i in 1:N
                # q_icm = RS[c,i,m]*Y[m]/P[c,i,m] # quantity
                SOC_mat[c,i] = (1+S_m[c])*Y_m/Q_m[c,i]^2*((1-σ)/σ^2*RS_m[c,i] + 
                    (3-2σ)*(σ-1)/σ^2*RS_m[c,i]^2+ 2*(σ-1)^2/σ^2*RS_m[c,i]^3)
        end
    end
    return SOC_mat
end

 
function Cournot_d2π_icm_dq2_icm_single_f(eq::eqbm_m, c::Int, i::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    # P = eq.P
    Y_m = eq.IDT.Y_m
    Q_m = eq.Q_m
    S_m = eq.IDT.S_m
    return (1+S_m[c])*Y_m/Q_m[c,i]^2*((1-σ)/σ^2*RS_m[c,i] + 
    (3-2σ)*(σ-1)/σ^2*RS_m[c,i]^2+ 2*(σ-1)^2/σ^2*RS_m[c,i]^3)
end


# Cournot SOC (outside good), market m 
function Cournot_d2π0_dq02_m_single_f(eq::eqbm_m)
    RS0_m = eq.RS0_m
    σ = eq.IDT.σ
    Y_m = eq.IDT.Y_m
    Q0_m = eq.Q0_m
    
    return Y_m/Q0_m^2*((1-σ)/σ^2*RS0_m + (3-2σ)*(σ-1)/σ^2*RS0_m^2+ 2*(σ-1)^2/σ^2*RS0_m^3)
end

# Cournot dπ_icm_dq_jkm (ic, jk inside good), market m
function Cournot_dπ_icm_dq_jkm_single_f(eq::eqbm_m, c::Int, i::Int, k::Int, j::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    S_m = eq.IDT.S_m

    return -(1+S_m[c])*(σ-1)/σ*RS_m[c,i]*P_m[k,j]
end

# Cournot dπ_icm_dq_0mm (ic-inside good, 0m-outside good), market m
function Cournot_dπ_icm_dq_0mm_single_f(eq::eqbm_m, c::Int, i::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P0_m = eq.P0_m
    S_m = eq.IDT.S_m

    return -(1+S_m[c])*(σ-1)/σ*RS_m[c,i]*P0_m
end


# Cournot dπ_0mm_dq_jkm (0m-outside good, jk outside good), market m
function Cournot_dπ_0mm_dq_jkm_single_f(eq::eqbm_m, k::Int, j::Int)
    RS0_m = eq.RS0_m
    σ = eq.IDT.σ
    P_m = eq.P_m

    return -(σ-1)/σ*RS0_m*P_m[k,j]
end

# Cournot d2π_icm_dq_icm_dq_jkm (ic, jk inside good), market m
function Cournot_d2π_icm_dq_icm_dq_jkm_single_f(eq::eqbm_m, c::Int, i::Int, k::Int, j::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    S_m = eq.IDT.S_m
    return -(1+S_m[c])/Y_m*(σ-1)^2/σ^2*P_m[c,i]*P_m[k,j]*(1-2*RS_m[c,i])
end


# Cournot d2π_icm_dq_icm_dq_0mm (ic-inside good, 0m-outside good), market m
function Cournot_d2π_icm_dq_icm_dq_0mm_single_f(eq::eqbm_m, c::Int, i::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    P0_m = eq.P0_m
    Y_m = eq.IDT.Y_m
    S_m = eq.IDT.S_m
    return -(1+S_m[c])/Y_m*(σ-1)^2/σ^2*P_m[c,i]*P0_m*(1-2*RS_m[c,i])
end


# Cournot d2π_0mm_dq_0mm_dq_jkm (jk-inside good, 0m-outside good), market m
function Cournot_d2π_0mm_dq_0mm_dq_jkm_single_f(eq::eqbm_m, k::Int, j::Int)
    RS0_m = eq.RS0_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    P0_m = eq.P0_m
    Y_m = eq.IDT.Y_m
    return -1/Y_m*(σ-1)^2/σ^2*P0_m*P_m[k,j]*(1-2*RS0_m)
end


function Cournot_d2π_icm_dq_icm_ds_cm(eq::eqbm_m, c::Int, i::Int)  
    Y_m = eq.IDT.Y_m
    RS_m = eq.RS_m
    P_m = eq.P_m
    σ = eq.IDT.σ
    
    return (σ-1)/σ*P_m[c,i]*(1-RS_m[c,i])
end

# Cournot dq_icm_dp_jkm (ic, jk inside good), market m
function Cournot_dp_icm_dq_jkm_single_f(eq::eqbm_m, c::Int, i::Int, k::Int, j::Int)
    # RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    
    return -(σ-1)/σ * P_m[c,i]*P_m[k,j]/Y_m
end

function Cournot_dp_icm_dq_0mm_single_f(eq::eqbm_m, c::Int, i::Int)
    # RS_m = eq.RS_m
    σ = eq.IDT.σ
    P0_m = eq.P0_m
    Y_m = eq.IDT.Y_m
    P_m = eq.P_m

    return -(σ-1)/σ * P_m[c,i]*P0_m/Y_m
end


function Cournot_dp_icm_dq_icm_single_f(eq::eqbm_m, c::Int, i::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    # P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    Q_m = eq.Q_m
    
    return - 1/σ*Y_m/Q_m[c,i]^2*(RS_m[c,i] + (σ-1)*RS_m[c,i]^2)
end


# ======================================== Bertrand ========================================

# Bertrand FOC (inside goods), market m 
function Bertrand_dπ_dp_m_single_f(eq::eqbm_m)
    C = eq.IDT.C
    N = eq.IDT.N
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    MC_m = eq.MC_m
    Y_m = eq.IDT.Y_m
    S_m = eq.IDT.S_m

    FOC_mat = zeros(C, N)
    for c in 1:C
        for i in 1:N
                FOC_mat[c,i] = (1+S_m[c])*Y_m*RS_m[c,i]/P_m[c,i]*(1-((1+S_m[c])*P_m[c,i]-MC_m[c,i])/((1+S_m[c])*P_m[c,i])*(σ+(1-σ)*RS_m[c,i]))
        end
    end
    
    return FOC_mat
end


# Bertrand FOC (outside good), market m 
function Bertrand_dπ0_dp0_m_single_f(eq::eqbm_m)
    RS0_m = eq.RS0_m
    σ = eq.IDT.σ
    P0_m = eq.P0_m
    MC0_m = eq.MC0_m
    Y_m = eq.IDT.Y_m
    
    return Y_m*RS0_m/P0_m*(1-(P0_m-MC0_m)/P0_m*(σ+(1-σ)*RS0_m))
end


# Bertrand SOC, (inside goods), market m 
function Bertrand_d2π_dp2_single_f(eq::eqbm_m)
    C = eq.IDT.C
    N = eq.IDT.N
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    MC_m = eq.MC_m
    S_m = eq.IDT.S_m

    SOC_vec = zeros(C, N)
    for c in 1:C
        for i in 1:N
            SOC_vec[c,i] = -2*(1+S_m[c])*Y_m*RS_m[c,i]/P_m[c,i]^2*(σ+(1-σ)*RS_m[c,i]) + ((1+S_m[c])*P_m[c,i]-MC_m[c,i])*Y_m/P_m[c,i]^3*(σ*(σ+1)*RS_m[c,i]+3*σ*(1-σ)*RS_m[c,i]^2+2*(1-σ)^2*RS_m[c,i]^3)
        end
    end
    
    return SOC_vec
end


function Bertrand_d2π_icm_dp2_icm_single_f(eq::eqbm_m, c::Int, i::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    MC_m = eq.MC_m
    S_m = eq.IDT.S_m

    return -2*(1+S_m[c])*Y_m*RS_m[c,i]/P_m[c,i]^2*(σ+(1-σ)*RS_m[c,i]) + ((1+S_m[c])*P_m[c,i]-MC_m[c,i])*Y_m/P_m[c,i]^3*(σ*(σ+1)*RS_m[c,i]+3*σ*(1-σ)*RS_m[c,i]^2+2*(1-σ)^2*RS_m[c,i]^3)
end


# Bertrand SOC (outside good), market m 
function Bertrand_d2π0_dp02_m_single_f(eq::eqbm_m)
    RS0_m = eq.RS0_m
    σ = eq.IDT.σ
    P0_m = eq.P0_m
    Y_m = eq.IDT.Y_m
    MC0_m = eq.MC0_m
    
    return -2*Y_m*RS0_m/P0_m^2*(σ+(1-σ)*RS0_m) + (P0_m-MC0_m)*Y_m/P0_m^3*(σ*(σ+1)*RS0_m+3*σ*(1-σ)*RS0_m^2+2*(1-σ)^2*RS0_m^3)
end


# Bertrand dπ_icm_dq_jkm (ic, jk inside good), market m
function Bertrand_dπ_icm_dp_jkm_single_f(eq::eqbm_m, c::Int, i::Int, k::Int, j::Int)
    Y_m = eq.IDT.Y_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Q_m = eq.Q_m
    MC_m = eq.MC_m
    S_m = eq.IDT.S_m

    return (σ-1)*((1+S_m[c])*P_m[c,i]-MC_m[c,i])/Y_m*Q_m[c,i]*Q_m[k,j]
end


# Bertrand dπ_icm_dq_0mm (ic-inside good, 0m-outside good), market m
function Bertrand_dπ_icm_dp_0mm_single_f(eq::eqbm_m, c::Int, i::Int)
    Y_m = eq.IDT.Y_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Q_m = eq.Q_m
    MC_m = eq.MC_m
    Q0_m = eq.Q0_m
    S_m = eq.IDT.S_m

    return (σ-1)*((1+S_m[c])*P_m[c,i]-MC_m[c,i])/Y_m*Q_m[c,i]*Q0_m
end


# Bertrand dπ_0mm_dp_jkm (0m-outside good, jk outside good), market m
function Bertrand_dπ_0mm_dp_jkm_single_f(eq::eqbm_m, k::Int, j::Int)
    Y_m = eq.IDT.Y_m
    σ = eq.IDT.σ
    Q_m = eq.Q_m
    P0_m = eq.P0_m
    Q0_m = eq.Q0_m
    MC0_m = eq.MC0_m

    return (σ-1)*(P0_m-MC0_m)/Y_m*Q0_m*Q_m[k,j]
end


# Bertrand d2π_icm_dp_icm_dp_jkm (ic, jk inside good), market m
function Bertrand_d2π_icm_dp_icm_dp_jkm_single_f(eq::eqbm_m, c::Int, i::Int, k::Int, j::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    MC_m = eq.MC_m
    S_m = eq.IDT.S_m

   return (1-σ)^2*Y_m*((1+S_m[c])*P_m[c,i]-MC_m[c,i])/P_m[c,i]^2*RS_m[c,i]^2*RS_m[k,j]/P_m[k,j]
   #return (σ-1)*Y_m/(P_m[c,i]*P_m[k,j])*RS_m[c,i]*RS_m[k,j] +(P_m[c,i]-MC_m[c,i])*(-σ*(σ-1)*Y_m*RS_m[c,i]*RS_m[k,j]/(P_m[c,i]^2*P_m[k,j])+ 2*(σ-1)^2*Y_m*RS_m[c,i]^2*RS_m[k,j]/(P_m[c,i]^2*P_m[k,j]))
end

# Bertrand d2π_icm_dp_icm_dp_0mm (ic-inside good, om-outside good), market m
function Bertrand_d2π_icm_dp_icm_dp_0mm_single_f(eq::eqbm_m, c::Int, i::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    MC_m = eq.MC_m
    P0_m = eq.P0_m
    RS0_m = eq.RS0_m
    S_m = eq.IDT.S_m

    return (1-σ)^2*Y_m*((1+S_m[c])*P_m[c,i]-MC_m[c,i])/P_m[c,i]^2*RS_m[c,i]^2*RS0_m/P0_m
end


# Bertrand d2π_0mm_dp_0mm_dp_jkm (jk-inside good, om-outside good), market m
function Bertrand_d2π_0mm_dp_0mm_dp_jkm_single_f(eq::eqbm_m, k::Int, j::Int)
    RS0_m = eq.RS0_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    P0_m = eq.P0_m
    MC0_m = eq.MC0_m
    RS_m = eq.RS_m

    return (1-σ)^2*Y_m*(P0_m-MC0_m)/P0_m^2*RS0_m^2*RS_m[k,j]/P_m[k,j]
end

# Bertrand dq_dp (inside goods), market m
function Bertrand_dq_dp_single_f(eq::eqbm_m)
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


function Bertrand_dq_icm_dp_icm_single_f(eq::eqbm_m, c::Int, i::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    return -Y_m*RS_m[c,i]/P_m[c,i]^2*(σ+(1-σ)*RS_m[c,i])
end

# Bertrand dq0_dp0 (outside good), market m
function Bertrand_dq0_dp0_single_f(eq::eqbm_m)
    RS0_m = eq.RS0_m
    σ = eq.IDT.σ
    P0_m = eq.P0_m
    Y_m = eq.IDT.Y_m
    dq0_dp0 = -Y_m*RS0_m/P0_m^2*(σ+(1-σ)*RS0_m)
    return dq0_dp0
end


# Bertrand dq_icm_dp_jkm (ic, jk inside good), market m
function Bertrand_dq_icm_dp_jkm_single_f(eq::eqbm_m, c::Int, i::Int, k::Int, j::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    
    return (σ-1)*Y_m*RS_m[c,i]*RS_m[k,j]/(P_m[c,i]*P_m[k,j])
end

# Bertrand dq_icm_dp_0mm (ic-inside good, 0m-outside good), market m
function Bertrand_dq_icm_dp_0mm_single_f(eq::eqbm_m, c::Int, i::Int)
    RS_m = eq.RS_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    P0_m = eq.P0_m
    RS0_m = eq.RS0_m
    
    return (σ-1)*Y_m*RS_m[c,i]*RS0_m/(P_m[c,i]*P0_m)
end 

# Bertrand dq_0mm_dp_jkm (0m-outside good, jk outside good), market m
function Bertrand_dq_0mm_dp_jkm_single_f(eq::eqbm_m, k::Int, j::Int)
    RS0_m = eq.RS0_m
    σ = eq.IDT.σ
    P_m = eq.P_m
    Y_m = eq.IDT.Y_m
    P0_m = eq.P0_m
    RS_m = eq.RS_m  
    
    return (σ-1)*Y_m*RS0_m*RS_m[k,j]/(P0_m*P_m[k,j])
end


function Bertrand_d2π_icm_dp_icm_ds_cm(eq::eqbm_m, c::Int, i::Int)
    Y_m = eq.IDT.Y_m
    RS_m = eq.RS_m
    P_m = eq.P_m
    σ = eq.IDT.σ
    
    return Y_m*RS_m[c,i]/P_m[c,i]*(1-σ)*(1-RS_m[c,i])
end