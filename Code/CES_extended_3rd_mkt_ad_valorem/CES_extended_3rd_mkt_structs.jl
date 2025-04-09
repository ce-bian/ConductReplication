
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
