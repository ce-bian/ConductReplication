# no need to simulate any data, but input the data directly
# Given parameters, solve for equilibrium prices and quantities in the CES model

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

