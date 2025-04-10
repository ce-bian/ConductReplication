# ============ structs for CES_MC_DGP ============
@with_kw mutable struct global_param # Global parameters
    C::Int64
    N::Int64
    M::Int64
    μ_Y::Float64
    V_Y::Float64
    μ_ξ::Float64
    V_ξ::Float64
    μ_ξ0::Float64
    V_ξ0::Float64
    μ_ω::Float64
    V_ω::Float64
    μ_ω0::Float64 # add
    V_ω0::Float64 # add
    ρ:: Float64 = 0.5
    off1::Float64
    off2::Float64
    off10::Float64
    off20::Float64
    μ_psi::Float64
    V_psi::Float64
    β::Float64
    σ::Float64
    w::Float64 = 0.7
    MAXIT::Int64 = 2000
    TOL::Float64 = 1e-8
    MAX_RETRIES::Int64 = 300
end


struct sim_dt # simulated exogenous parameters + global parameters
    GP::global_param
    Y::Array{Float64,1}
    exp_ξ::Array{Float64,3}
    exp_ω::Array{Float64,3}
    exp_ξ0::Vector{Float64} # for outside good
    exp_ω0::Vector{Float64} # for outside good
    τ::Array{Float64,2}
end

struct eqbm #equilibrium objects + global parameters + simulated exogenous parameters
    GP::global_param
    R::String
    P :: Array{Float64,3}
    P0:: Array{Float64,1}
    RS :: Array{Float64,3}
    RS0:: Array{Float64,1}
    Y :: Array{Float64,1}
    τ :: Array{Float64,2}
    herfindahl :: Float64
    num_retry :: Int
end

# for test purpose only, memory heavy
struct eqbm_full #equilibrium objects + global parameters + simulated exogenous parameters
    GP::global_param
    R::String
    P :: Array{Float64,3}
    P0:: Array{Float64,1}
    RS :: Array{Float64,3}
    RS0:: Array{Float64,1}
    Y :: Array{Float64,1}
    τ :: Array{Float64,2}
    herfindahl :: Float64
    num_retry :: Int
    exp_ξ ::Array{Float64,3}
    exp_ω :: Array{Float64,3}
    exp_ξ0 :: Vector{Float64} # for outside good
    exp_ω0 :: Vector{Float64}
    RV :: Array{Float64,3}
    RV0 :: Array{Float64,1}
end



# ============ structs for CES_extended_3rd_mkt_ad_valorem (with OG)  ============

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
    S_m::Vector{Float64} # ad valorem subsidy/tax; change: make it a vector! size = C
end


# equilibrium object for a single market
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


mutable struct eqbm_t
    eqbm_m_list::Array{eqbm_m,1}
    t :: Int64
end


mutable struct eqbm_output
    seed :: Int64
    eqbm_t_list::Array{eqbm_t,1}
end



# ============ structs for CES_3rd_mkt_ad_valorem (no OG) ============

@with_kw mutable struct input_dt_no_OG #data
    C::Int64
    N::Int64
    M::Int64
    σ::Float64
    w::Float64 = 0.7
    MAXIT::Int64 = 2000
    TOL::Float64 = 1e-8
    MAX_RETRIES::Int64 = 300
    Y::Array{Float64,1}
    exp_ξ::Array{Float64,3}
    exp_ω::Array{Float64,3}
    τ::Array{Float64,2}
    S::Array{Float64,3} # subsidy/tax
end

@with_kw mutable struct input_dt_m_no_OG #data
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
    τ::Array{Float64,2}
    S_m::Array{Float64,2} # subsidy/tax
end


struct eqbm_no_OG #equilibrium objects + input data
    IDT::input_dt_no_OG
    R::String
    P :: Array{Float64,3}
    Q :: Array{Float64,3}
    RS :: Array{Float64,3}
    S :: Array{Float64,3}
    MC:: Array{Float64,3}
    Y :: Array{Float64,1}
    τ :: Array{Float64,2}
    herfindahl :: Float64
end


# memory heavy
struct eqbm_full_no_OG #equilibrium objects + input data
    IDT::input_dt_no_OG
    R::String
    P :: Array{Float64,3}
    Q :: Array{Float64,3}
    RS :: Array{Float64,3}
    MC:: Array{Float64,3}
    Y :: Array{Float64,1}
    τ :: Array{Float64,2}
    herfindahl :: Float64
    exp_ξ ::Array{Float64,3}
    exp_ω :: Array{Float64,3}
    S :: Array{Float64,3}
    RV :: Array{Float64,3}
end

# equilibrium object for a single market
struct eqbm_m_no_OG #equilibrium objects + input data
    IDT::input_dt_m_no_OG
    R::String
    P_m :: Array{Float64,2}
    Q_m :: Array{Float64,2}
    RS_m :: Array{Float64,2}
    MC_m:: Array{Float64,2}
    Y_m :: Float64
    τ :: Array{Float64,2}
    herfindahl :: Float64
    exp_ξ_m ::Array{Float64,2}
    exp_ω_m :: Array{Float64,2}
    S_m :: Array{Float64,2}
    RV_m :: Array{Float64,2}
end

