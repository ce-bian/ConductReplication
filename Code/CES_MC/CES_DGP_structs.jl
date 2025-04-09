@with_kw mutable struct global_param # Global parameters
    C::Int64
    N::Int64
    M::Int64
    μ_Y::Float64
    V_Y::Float64
    μ_ξ::Float64
    V_ξ::Float64
    μ_ω::Float64
    V_ω::Float64
    ρ:: Float64 = 0.5
    off1::Float64
    off2::Float64
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

