# output simulated data to JLD2 files
# at this point try to be consistent with CES_solve_eqba_m_functions so that we can directly apply functions there for welfare analysis
# later will need to merge the two set of files to integrate simulations, estimations and policy analysis (to be done)
@with_kw mutable struct input_dt_m #data for a single market
    C::Int64
    N::Int64
    m::Int64
    σ::Float64
    w::Float64 = 0.7
    MAXIT::Int64 = 2000
    TOL::Float64 = 1e-8
    MAX_RETRIES::Int64 = 300
    Y_m::Float64
    exp_ξ_m::Array{Float64,2}
    exp_ω_m::Array{Float64,2}
    exp_ξ0_m::Float64 # for outside good
    exp_ω0_m::Float64 # for outside good
    τ::Array{Float64,2}
    S_m::Vector{Float64} # subsidy/tax; change: make it a vector! size = C
end


struct eqbm_m #equilibrium objects + input data
    IDT::input_dt_m
    R::String
    P_m :: Array{Float64,2}
    P0_m :: Float64
    Q_m :: Array{Float64,2}
    Q0_m :: Float64
    RS_m :: Array{Float64,2}
    RS0_m :: Float64
    MC_m:: Array{Float64,2}
    MC0_m:: Float64
    herfindahl :: Float64
    RV_m :: Array{Float64,2}
    RV0_m :: Float64
end


struct eqbm_t
    eqbm_m_list::Array{eqbm_m,1}
    t :: Int64
end


struct eqbm_output
    seed :: Int64
    eqbm_t_list::Array{eqbm_t,1}
end



function DGP_output_t(t::Int64, R::String, GP::global_param, global_seed::Int64)
    """
    Output the simulated data in eqbm_output_t format 
    Input:
        - R: "B" (Bertrand) or "C" (Cournot)
        - t: task index
        - GP: global parameters
        - global_seed: global seed for random number generator
    Output:
        - eqbm_t
    """
    C = GP.C
    N = GP.N
    M = GP.M
    σ = GP.σ
    MAX_RETRIES = GP.MAX_RETRIES
    w = GP.w
    TOL = GP.TOL
    MAXIT = GP.MAXIT


    unique_seed = global_seed + t # Create a unique seed for this task
    task_rng = MersenneTwister(unique_seed)  # Initialize RNG with the task-specific seed
    # Y, exp_ξ, exp_ω, τ = set_θ(GP, task_rng)
    SDT = set_θ(GP, task_rng)
    exp_ξ = SDT.exp_ξ
    exp_ω = SDT.exp_ω
    exp_ξ0 = SDT.exp_ξ0
    exp_ω0 = SDT.exp_ω0
    τ = SDT.τ
    Y = SDT.Y


    eqba_m_list = Array{eqbm_m,1}(undef, C+M) # include both C and M markets

    m = 1
    num_retry = 0

    while m <= C+M
        retries = 0
        while true
            try
                P_m, P0_m = solve_P_m_R(R, m, C, N, exp_ξ[:, :, m], exp_ω[:, :, m], exp_ξ0[m], exp_ω0[m], τ, σ, w, MAXIT, TOL)

                RS_mat = [RS_icm(c, i, P_m, exp_ξ[:, :, m], P0_m, exp_ξ0[m], σ) for c in 1:C, i in 1:N]
                RS0_m = RS_0mm(P_m, exp_ξ[:, :, m], P0_m, exp_ξ0[m], σ)
                S_m = zeros(C) # when simulating data, we assume no subsidy/tax
                MC_m = [MC_icm(c, i, m, τ, exp_ω[:, :, m]) for c in 1:C, i in 1:N]
                MC0_m = MC_0mm(τ, exp_ω0[m])
                herfindahl = sum(RS_mat .^ 2) + RS0_m^2
                RV_m = [RS_mat[c, i] * Y[m] for c in 1:C, i in 1:N]
                RV0_m = RS0_m * Y[m]
                Q_m = [RS_mat[c, i] * Y[m] / P_m[c, i] for c in 1:C, i in 1:N]
                Q0_m = RS0_m * Y[m] / P0_m
                # construct eqbm_m and input_dt_m
                IDT_m = input_dt_m(C=C, N=N, m=m, σ=σ, w=w, MAXIT=MAXIT, TOL=TOL, MAX_RETRIES=MAX_RETRIES, Y_m=Y[m], exp_ξ_m=exp_ξ[:, :, m], exp_ξ0_m=exp_ξ0[m], exp_ω_m=exp_ω[:, :, m], exp_ω0_m=exp_ω0[m], τ=τ, S_m=S_m)
                eqba_m_list[m] = eqbm_m(IDT_m, R, P_m, P0_m, Q_m, Q0_m, RS_mat, RS0_m, MC_m, MC0_m, herfindahl, RV_m, RV0_m)

                break # Break the loop if the above line executes successfully
            catch e
                num_retry += 1
                if isa(e, AssertionError)
                    retries += 1
                    # println("AssertionError caught. Retry $retries...")
                    if retries > MAX_RETRIES
                        error("Maximum retries reached. Unable to find a valid solution.")
                    end
                    unique_seed_new = global_seed + t + 1000 * num_retry # Create a unique seed for this task
                    task_rng_new = MersenneTwister(unique_seed_new)  # Initialize RNG with the task-specific seed
                    # Y, exp_ξ, exp_ω, τ = set_θ(GP, task_rng)
                    SDT = set_θ_update_m(m, SDT, GP, task_rng_new) # # update original SDT
                    exp_ξ = SDT.exp_ξ
                    exp_ω = SDT.exp_ω
                    exp_ξ0 = SDT.exp_ξ0
                    exp_ω0 = SDT.exp_ω0
                    τ = SDT.τ
                    Y = SDT.Y
                else
                    rethrow()
                end
            end # end for try-catch
        end # end for while true

        m += 1
    end # end for while m <= C+M
    # add
   
    eqbm_t_output = eqbm_t(eqba_m_list, t)
    return eqbm_t_output
end



function DGP_output(R::String, seed::Int64, GP::global_param, t_list::Vector{Int64})
    """
    Output the simulated data in eqbm_output format 
    Input:
        - R: "B" (Bertrand) or "C" (Cournot)
        - seed: seed for random number generator
        - GP: global parameters
        - t_list: list of periods 
    Output:
        - eqbm_output
    """
    T_length = length(t_list)
    eqbm_t_list = Array{eqbm_t,1}(undef, T_length)
    # construct a Dict to match each element in t_list to 1:T_length
    T_dict = Dict(t_list[i] => i for i in 1:T_length)
    for t in t_list
        eqbm_t_list[T_dict[t]] = DGP_output_t(t, R, GP, seed)
    end
    # for t in t_list
    #     eqbm_t_list[t] = DGP_output_t(t, R, GP, seed)
    # end

    eqbm_output_R = eqbm_output(seed, eqbm_t_list)
    return eqbm_output_R
end




function output_JLD2_file(seed::Int, inputFile::String, outputFile::String)
    """
    Save the output to a JLD2 file for a given set of parameters
    """
    dfB, dfC, GP = load(inputFile, "dfB", "dfC", "GP")
    # get the unique t_list 
    t_list_B = unique(dfB.t)
    t_list_C = unique(dfC.t)

    eqbm_output_B = DGP_output("B", seed, GP, t_list_B)
    eqbm_output_C = DGP_output("C", seed, GP, t_list_C)

    save(outputFile, "eqbm_output_B", eqbm_output_B, "eqbm_output_C", eqbm_output_C, "GP", GP, "seed", seed)
end


