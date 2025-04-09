# functions to sample subset of markets from the full dataset
function sample_mkt_t_f(eq_full::eqbm_t, mkt_num::Int64)
    eqbm_m_list = eq_full.eqbm_m_list
    # keep the first mkt_num markets, later modify to keep random m markets
    C = eq_full.eqbm_m_list[1].IDT.C
    eqbm_m_list = eqbm_m_list[1:C+mkt_num]
    eq_subset = eqbm_t(eqbm_m_list, eq_full.t)
    return eq_subset
end


function sample_mkt_f(eq_full::eqbm_output, mkt_num::Int64)
    eqbm_t_list = eq_full.eqbm_t_list
    eqbm_t_list = [sample_mkt_t_f(eqbm_t, mkt_num) for eqbm_t in eqbm_t_list]
    eq_subset = eqbm_output(eq_full.seed, eqbm_t_list)
    return eq_subset
end



