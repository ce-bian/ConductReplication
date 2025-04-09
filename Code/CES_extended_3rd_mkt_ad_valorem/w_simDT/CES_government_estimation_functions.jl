# functions that estimate tau, normalized \xi, \omega with given observables and estimated \sigma (GMM_DSB or GMM_DSC depending on the inferred conduct)

# markups given market shares
μB_icm_R(c::Int,i::Int,RS::Array{Float64,2},σ) = 1/((σ-1)*(1-RS[c,i])) + 1
μC_icm_R(c::Int,i::Int,RS::Array{Float64,2},σ) = σ/((σ-1)*(1-RS[c,i]))
μB_0mm_R(RS0::Float64,σ) = 1/((σ-1)*(1-RS0)) + 1
μC_0mm_R(RS0::Float64,σ) = σ/((σ-1)*(1-RS0))

function estimated_eqbm_m_f(R_inferred::String, eq_data::eqbm_m, σ_est::Float64)
    C = eq_data.IDT.C
    N = eq_data.IDT.N
    m = eq_data.IDT.m
    # R = eq_data.R  # true R in DGP

    # observables to government
    P_m = eq_data.P_m
    P0_m = eq_data.P0_m
    RS_mat = eq_data.RS_m 
    RS0_mat = eq_data.RS0_m
    τ = eq_data.IDT.τ
    S_m = eq_data.IDT.S_m

    μR_icm = R_inferred == "B" ? μB_icm_R : μC_icm_R
    μR_0mm = R_inferred == "B" ? μB_0mm_R : μC_0mm_R
    
    est_μ_m = [μR_icm(c, i, RS_mat, σ_est) for c in 1:C, i in 1:N]
    est_μ0_m = μR_0mm(RS0_mat, σ_est)
    est_exp_ω_m = [(P_m[c,i]/est_μ_m[c,i]+S_m[c])/τ[c,m] for c in 1:C, i in 1:N]
    est_exp_ω0_m = P0_m/est_μ0_m

    RS_normalized_m = [RS_mat[c,i]/RS0_mat for c in 1:C, i in 1:N]
    P_normalized_m = [P_m[c,i]/P0_m for c in 1:C, i in 1:N]
    
    est_ξ0_m = 0.0
    # use demand equations to get estimated ξ_m (normalized)
    est_ξ_m = [log(RS_normalized_m[c,i]) - (1-σ_est)*log(P_normalized_m[c,i])  for c in 1:C, i in 1:N]

    est_IDT = deepcopy(eq_data.IDT)
    est_IDT.exp_ξ_m = exp.(est_ξ_m)
    est_IDT.exp_ω_m = est_exp_ω_m
    est_IDT.exp_ξ0_m =exp(est_ξ0_m)
    est_IDT.exp_ω0_m = est_exp_ω0_m
    est_IDT.σ = σ_est # update σ to estimated value

    est_eq_m = deepcopy(eq_data)
    est_eq_m.IDT = est_IDT
    
    return est_eq_m
end


function estimated_eqbm_t_f(R_inferred::String, eq_data::eqbm_t, σ_est::Float64)
    t = eq_data.t
    eqbm_m_list = Array{eqbm_m,1}(undef, length(eq_data.eqbm_m_list))
    for i in 1:length(eq_data.eqbm_m_list)
        eqbm_m = eq_data.eqbm_m_list[i]
        eqbm_m_est = estimated_eqbm_m_f(R_inferred, eqbm_m, σ_est)
        eqbm_m_list[i] = eqbm_m_est
    end
    return eqbm_t(eqbm_m_list, t)
end


# function estimated_eqba_output_f(R_inferred::String, eq_data::eqba_output, σ_est::Float64)
#     seed = eq_data.seed
#     eqbm_t_list = Array{eqbm_t,1}(undef, length(eq_data.eqbm_t_list))
#     for i in 1:length(eq_data.eqbm_t_list)
#         eqbm_t = eq_data.eqbm_t_list[i]
#         eqbm_t_est = estimated_eqbm_t_f(R_inferred, eqbm_t, σ_est)
#         eqbm_t_list[i] = eqbm_t_est
#     end
#     return eqbm_output(seed, eqbm_t_list)
# end


