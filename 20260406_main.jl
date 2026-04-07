using Distributions
using LsqFit
using Plots
using Parameters: @unpack
using QuantEcon: tauchen, stationary_distributions, MarkovChain, Categorical
using LaTeXStrings
using Measures
using JLD2: @save, @load
using Random
using FreqTables
using BenchmarkTools
using QuadGK
# using GLMakie
# using CairoMakie
using Polyester
using Optim
using StatsPlots

function adda_cooper(N::Integer, ρ::Real, σ::Real; μ::Real=0.0)
    """
    Approximation of an autoregression process with a Markov chain proposed by Adda and Cooper (2003)
    """
    σ_ϵ = σ / sqrt(1.0 - ρ^2.0)
    q = quantile(Normal(), range(0, 1; length=N + 1))
    ϵ = σ_ϵ .* q .+ μ
    z = [N * σ_ϵ * (pdf(Normal(), (ϵ[i] - μ) / σ_ϵ) - pdf(Normal(), (ϵ[i+1] - μ) / σ_ϵ)) + μ for i ∈ 1:N]

    Π = zeros(N, N)
    dist_σ = Normal(μ, σ_ϵ)
    if ρ == 0.0
        Π .= 1.0 / N
    else
        for i ∈ 1:N, j ∈ 1:N
            f(u) = pdf(dist_σ, u) * (cdf(Normal(), (ϵ[j+1] - μ * (1.0 - ρ) - ρ * u) / σ) - cdf(Normal(), (ϵ[j] - μ * (1.0 - ρ) - ρ * u) / σ))
            Π[i, j], _ = quadgk(u -> f(u), ϵ[i], ϵ[i+1])
            Π[i, j] *= N
        end
    end
    Π = Π ./ sum(Π, dims=2)
    return z, Π
end

function tauchen_grid(N::Int, ρ::Float64, σ::Float64; μ::Float64=0.0, m::Float64=3.0)
    σ_z = σ / sqrt(1 - ρ^2)
    z_max = μ + m * σ_z
    z_min = μ - m * σ_z
    z = collect(range(z_min, z_max; length=N))
    return z
end

function tauchen_transition_matrix(z::Vector{Float64}, ρ::Float64, σ::Float64; μ::Float64=0.0)
    N = length(z)
    Δ = z[2] - z[1]  # uniform spacing
    Π = zeros(N, N)

    for i in 1:N
        μ_cond = μ * (1 - ρ) + ρ * z[i]
        for j in 1:N
            if j == 1
                Π[i, j] = cdf(Normal(μ_cond, σ), z[1] + Δ / 2)
            elseif j == N
                Π[i, j] = 1.0 - cdf(Normal(μ_cond, σ), z[N] - Δ / 2)
            else
                Π[i, j] = cdf(Normal(μ_cond, σ), z[j] + Δ / 2) -
                          cdf(Normal(μ_cond, σ), z[j] - Δ / 2)
            end
        end
    end
    return Π
end


function binomial_matrix_function(n_max::Integer, p::Real)
    """
    Construct transition matrix of dependent children (up to n_max)
    	n_max - largest numder of children considerd
    	p     - probability that a child becomes independent
    """
    Π_n = zeros(1 + n_max, 1 + n_max)
    Π_n[1, 1] = 1.0
    for i ∈ 1:n_max
        d = Binomial(i, 1.0 - p)
        for j ∈ 0:i
            Π_n[1+i, 1+j] = pdf(d, j)
        end
    end
    return Π_n
end

function infertility_risk_function(data_age::Array{Int64,1}, data_inf::Array{Float64,1}, age_min::Integer, age_max::Integer, age_inf::Integer)
    """
    Exponential fit of infertility probability, intrapolated on ages up to age_inf
    """
    model(t, ω) = ω[1] * exp.(ω[2] * t)
    ω_int = [0.5, 0.5]
    fit = curve_fit(model, data_age, data_inf, ω_int)
    model_age = collect(age_min:age_max)
    model_inf = fit.param[1] .* exp.(fit.param[2] .* model_age)
    model_inf[findfirst(model_age .== age_inf):end] .= 1.0
    return model_age, model_inf
end

function infertility_risk_low_function(data_age::Array{Int64,1}, data_inf::Array{Float64,1}, age_min::Integer, age_max::Integer)
    """
    Exponential fit of infertility probability, intrapolated on ages up to age_inf
    """
    model(t, ω) = ω[1] * exp.(ω[2] * t)
    ω_int = [0.5, 0.5]
    fit = curve_fit(model, data_age, data_inf, ω_int)
    model_age = collect(age_min:age_max)
    model_inf = fit.param[1] .* exp.(fit.param[2] .* model_age)
    age_inf = model_age[findlast(model_inf .< 1.0)]
    model_inf[findfirst(model_age .== (age_inf + 1)):end] .= 1.0
    return model_age, model_inf, age_inf
end

@inline function u_CRRA_e(c::Float64, e::Float64, one_m_γ::Float64, inv_one_m_γ::Float64, one_m_κ::Float64, ψ_inv_one_m_κ::Float64, scale::Float64; c_floor::Float64=1e-12, V_penalty::Float64=-1.0e16)
    c < c_floor && return V_penalty
    return c^one_m_γ * inv_one_m_γ + scale * e^one_m_κ * ψ_inv_one_m_κ
end

@inline function u_CRRA(c::Float64, one_m_γ::Float64, inv_one_m_γ::Float64; c_floor::Float64=1e-12, V_penalty::Float64=-1.0e16)
    c < c_floor && return V_penalty
    return c^one_m_γ * inv_one_m_γ
end

@inline function u_log(c::Float64; c_floor::Float64=1e-12, V_penalty::Float64=-1.0e16)
    c < c_floor && return V_penalty
    return log(c)
end

function quality_function(x::Real, l::Real, n::Integer, μ::Real, θ::Real, ψ_1::Real, ψ_2::Real)
    """
    child quality function
    	x - time input
    	l - monetary input
    	n - number of kids
    """
    if n > 0
        return (μ * (x / (n^ψ_1))^θ + (1.0 - μ) * (l / (n^ψ_2))^θ)^(1.0 / θ)
    else
        return 0.0
    end
end

function parameters_function(;
    #----------------------#
    # exogenous parameters #
    #----------------------#
    r::Real=0.04,                   # interest rate
    β::Real=1.0 / (1.0 + r),        # discount factor
    γ::Real=1.5,                    # risk aversion
    ρ::Real=0.95,                   # persistance coefficient
    σ_ϵ::Real=0.21,                 # std of persistent shock
    σ_ν::Real=0.17,                 # std of transitory shock
    b::Real=0.40,                   # replacement rate
    #----------------------#
    # estimated parameters #
    #----------------------#
    κ::Real=0.14,                   # preference curvature
    ψ::Real=3.50,                   # preference scale
    μ::Real=0.35,                   # production share
    θ::Real=0.70,                   # elasticity of substitution in production
    q_bar::Real=0.50,               # lower bound on children's consumption
    ψ_1::Real=0.91,                 # HH economies to money input to production
    ψ_2::Real=0.54,                 # HH economies to time input to production
    p::Real=0.02,                   # prob that a child becomes independent

    # numerical solution #

    age_min::Integer=18,            # min age
    age_max::Integer=80,            # max age
    age_edu::Integer=22,            # education age
    age_inf::Integer=45,            # infertile age
    age_ret::Integer=65,            # retirement age
    n_max::Integer=4,               # max number of kids
    ϵ_size::Integer=7,              # number of persistent shock
    ν_size::Integer=5,              # number of transitory shock
    a_min::Real=0.0,                # min of asset holding
    a_max::Real=120,                # max of asset holding
    a_size::Integer=250,            # number of asset
    a_degree::Integer=1,            # curvature of asset gridpoints
    q_x::Real=1.0,                  # price of monetary input $x$

    # case indicators #
    edu_ind::Symbol=:NC,          # with or without education
)
    """
    Contruct an immutable object containg all paramters
    """

    # auxiliary parameters
    R = 1.0 + r
    inv_R = 1.0 / R
    m_inv_γ = -1.0 / γ
    one_m_γ = 1.0 - γ
    inv_one_m_γ = 1.0 / (1.0 - γ)
    βR = β * R
    EGM_fac = βR^m_inv_γ
    one_m_κ = 1.0 - κ
    inv_κ = 1.0 / κ
    inv_one_m_κ = 1.0 / (1.0 - κ)
    ψ_inv_one_m_κ = ψ / (1.0 - κ)
    γ_by_κ = γ / κ

    # infertility parameters: taken from Trussell and Wilson (1985, Population Studies)
    data_inf = [0.07, 0.131, 0.231, 0.345, 0.576, 0.952]
    data_age = [20, 25, 30, 35, 40, 45]
    age_grid, inf_grid = infertility_risk_function(data_age, data_inf, age_min, age_max, age_inf)
    age_size = length(age_grid)
    inf_size = 2

    # transition of child dependence
    n_grid = collect(0.0:n_max)
    n_size = length(n_grid)
    n_Γ = binomial_matrix_function(n_max, p)

    # life-cycle income
    data_h_NC = [
        8.720850, 9.124881, 9.435658, 9.628212, 9.761456,
        9.860271, 9.921792, 9.983205, 10.033841, 10.065866,
        10.090777, 10.122985, 10.126819, 10.138425, 10.177699,
        10.181802, 10.186769, 10.225693, 10.207617, 10.228332,
        10.230438, 10.213399, 10.222961, 10.225275, 10.237479,
        10.246593, 10.237827, 10.257969, 10.243393, 10.239770,
        10.216180, 10.211378, 10.190308, 10.191437, 10.174690,
        10.164399, 10.186016, 10.153089, 10.128997, 10.122162,
        10.081263, 10.121118, 10.059273, 10.015279, 9.950545,
        9.838649, 9.781660, 9.637810,
    ]    

    data_h_C = [
        9.428264, 9.767262, 10.091448, 10.258264, 10.360586,
        10.471955, 10.512105, 10.589578, 10.640750, 10.681159,
        10.707237, 10.730122, 10.774949, 10.827320, 10.870387,
        10.916551, 10.863727, 10.879524, 10.890114, 10.931122,
        10.943369, 10.966267, 10.959625, 10.992401, 10.987035,
        10.958040, 10.978917, 11.066968, 11.017463, 11.017285,
        10.960259, 11.035227, 11.028443, 11.068266, 10.941890,
        10.945880, 10.893451, 10.880835, 10.917661, 10.824608,
        10.808046, 10.612711, 10.541795, 10.430953,
    ]
    h_mean = mean(data_h_NC)

    if edu_ind == :NC
        p_fertile_at_entry = 1.0
        h_grid = copy(data_h_NC) .- h_mean
        h_grid = vcat(h_grid, repeat([h_grid[end]], age_max - age_ret))
        h_size = length(h_grid)
        @assert age_size == h_size
    elseif edu_ind == :C
        age_min = age_min + 4
        age_grid = age_grid[5:end]
        age_size = length(age_grid)
        p_fertile_at_entry = prod(1.0 .- inf_grid[1:4])
        inf_grid = inf_grid[5:end]
        @assert age_size == length(inf_grid)
        h_grid = copy(data_h_C) .- h_mean
        h_grid = vcat(h_grid, repeat([h_grid[end]], age_max - age_ret))
        h_size = length(h_grid)
        @assert age_size == h_size
    end

    # persistent income shock
    ϵ_grid = tauchen_grid(ϵ_size, ρ, σ_ϵ; m=2.0)
    ϵ_Γ = tauchen_transition_matrix(ϵ_grid, ρ, σ_ϵ)
    # ϵ_grid, ϵ_Γ = adda_cooper(ϵ_size, ρ, σ_ϵ)
    ϵ_G = stationary_distributions(MarkovChain(ϵ_Γ, ϵ_grid))[1]

    # transitory income shock
    ν_grid = tauchen_grid(ν_size, 0.0, σ_ν; m=2.0)
    ν_Γ = tauchen_transition_matrix(ν_grid, 0.0, σ_ν)
    # ν_grid, ν_Γ = adda_cooper(ν_size, 0.0, σ_ν)
    ν_Γ = ν_Γ[1, :]
    ν_G = ν_Γ

    Γ_ret = zeros(ϵ_size, ν_size, ϵ_size)
    for ϵ_i in 1:ϵ_size, ν_p_i in 1:ν_size, ϵ_p_i in 1:ϵ_size
        Γ_ret[ϵ_p_i, ν_p_i, ϵ_i] = ν_Γ[ν_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i]
    end

    Γ_inf = zeros(n_size, ϵ_size, ν_size, n_size, ϵ_size)
    for ϵ_i in 1:ϵ_size, n_i in 1:n_size, ν_p_i in 1:ν_size, ϵ_p_i in 1:ϵ_size, n_p_i in 1:n_size
        Γ_inf[n_p_i, ϵ_p_i, ν_p_i, n_i, ϵ_i] = n_Γ[n_i, n_p_i] * ν_Γ[ν_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i]
    end

    Γ = zeros(n_size, ϵ_size, ν_size, inf_size, n_size, ϵ_size, inf_size, age_inf - age_min)
    for h_i in 1:(age_inf-age_min), ϵ_i in 1:ϵ_size, n_i in 1:n_size, ν_p_i in 1:ν_size, ϵ_p_i in 1:ϵ_size, n_p_i in 1:n_size
        Γ[n_p_i, ϵ_p_i, ν_p_i, 1, n_i, ϵ_i, 1, h_i] = (1.0 - inf_grid[h_i]) * n_Γ[n_i, n_p_i] * ν_Γ[ν_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i]
        Γ[n_p_i, ϵ_p_i, ν_p_i, 2, n_i, ϵ_i, 1, h_i] = inf_grid[h_i] * n_Γ[n_i, n_p_i] * ν_Γ[ν_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i]
        Γ[n_p_i, ϵ_p_i, ν_p_i, 2, n_i, ϵ_i, 2, h_i] = n_Γ[n_i, n_p_i] * ν_Γ[ν_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i]
    end

    w_grid = Array{Float64}(undef, ϵ_size, ν_size, h_size)
    for h_i ∈ 1:h_size, ν_i in 1:ν_size, ϵ_i in 1:ϵ_size
        h = h_grid[h_i]
        ν = ν_grid[ν_i]
        ϵ = ϵ_grid[ϵ_i]
        ret_idx = age_ret - age_min + 1
        w_grid[ϵ_i, ν_i, h_i] = h_i < ret_idx ? exp(h + ν + ϵ) : b * exp(h + ν + ϵ)
    end

    P_grid = Array{Float64}(undef, n_size, ϵ_size, ν_size, h_size)
    σ_θ = 1.0 / (1.0 - θ)
    inv1mσ = 1.0 / (1.0 - σ_θ)
    μσ = μ^σ_θ
    oneμσ = (1.0 - μ)^σ_θ
    pow1 = ψ_1 * (1.0 - σ_θ)
    pow2 = ψ_2 * (1.0 - σ_θ)
    for h_i in 1:h_size, ν_i in 1:ν_size, ϵ_i in 1:ϵ_size, n_i in 1:n_size
        n = n_grid[n_i]
        if n == 0.0
            P_grid[n_i, ϵ_i, ν_i, h_i] = 1.0
        else
            w = w_grid[ϵ_i, ν_i, h_i]
            term1 = μσ * n^pow1
            term2 = oneμσ * w^(1.0 - σ_θ) * n^pow2
            P_grid[n_i, ϵ_i, ν_i, h_i] = (term1 + term2)^inv1mσ
        end
    end
    q_bar_P_grid = q_bar .* P_grid

    # asset holding
    a_grid = ((range(0.0, stop=a_size - 1, length=a_size) / (a_size - 1)) .^ a_degree) * a_max
    a_ind_zero = 1
    aR_grid = a_grid * R

    # numerical guards
    c_floor = 1.0e-12
    V_penalty = -1.0e16

    # return values
    return (
        r=r,
        R=R,
        inv_R=inv_R,
        β=β,
        γ=γ,
        m_inv_γ=m_inv_γ,
        one_m_γ=one_m_γ,
        inv_one_m_γ=inv_one_m_γ,
        βR=βR,
        EGM_fac=EGM_fac,
        one_m_κ=one_m_κ,
        inv_κ=inv_κ,
        inv_one_m_κ=inv_one_m_κ,
        ψ_inv_one_m_κ=ψ_inv_one_m_κ,
        γ_by_κ=γ_by_κ,
        ρ=ρ,
        σ_ϵ=σ_ϵ,
        σ_ν=σ_ν,
        b=b,
        κ=κ,
        ψ=ψ,
        μ=μ,
        θ=θ,
        q_bar=q_bar,
        ψ_1=ψ_1,
        ψ_2=ψ_2,
        p=p,
        age_min=age_min,
        age_max=age_max,
        age_edu=age_edu,
        age_inf=age_inf,
        age_ret=age_ret,
        age_size=age_size,
        age_grid=age_grid,
        inf_grid=inf_grid,
        data_age=data_age,
        data_inf=data_inf,
        inf_size=inf_size,
        n_max=n_max,
        n_size=n_size,
        n_grid=n_grid,
        n_Γ=n_Γ,
        h_size=h_size,
        h_grid=h_grid,
        ϵ_size=ϵ_size,
        ϵ_grid=ϵ_grid,
        ϵ_Γ=ϵ_Γ,
        ϵ_G=ϵ_G,
        ν_size=ν_size,
        ν_grid=ν_grid,
        ν_Γ=ν_Γ,
        ν_G=ν_G,
        Γ_ret=Γ_ret,
        Γ_inf=Γ_inf,
        Γ=Γ,
        w_grid=w_grid,
        P_grid=P_grid,
        q_bar_P_grid=q_bar_P_grid,
        a_min=a_min,
        a_max=a_max,
        a_ind_zero=a_ind_zero,
        a_size=a_size,
        a_grid=a_grid,
        aR_grid=aR_grid,
        a_degree=a_degree,
        q_x=q_x,
        c_floor=c_floor,
        V_penalty=V_penalty,
        p_fertile_at_entry=p_fertile_at_entry,
    )
end

mutable struct EGMWorkspace
    V_endo::Vector{Float64}
    a_endo::Vector{Float64}
    ap_endo::Vector{Float64}
    V1::Vector{Float64}
    c1::Vector{Float64}
    ap1::Vector{Float64}
    e1::Vector{Float64}
    EV_next::Array{Float64,4}
    c_next_CE::Array{Float64,4}
end

mutable struct Mutable_Variables
    """
    Construct a type for mutable variables
    """
    V::Array{Float64,6}
    policy_c::Array{Float64,6}
    policy_a_p::Array{Float64,6}
    policy_e::Array{Float64,6}
    policy_K::Array{Float64,6}
    EGM_ws::EGMWorkspace
end

function variables_function(parameters::NamedTuple)
    """
    Construct a mutable object containing endogenous variables
    """
    # unpack parameters
    @unpack inf_size, a_size, n_size, ϵ_size, ν_size, age_size = parameters

    # define value and policy functions: (a,n,ϵ,ν,f,d,age)
    V = zeros(a_size, n_size, ϵ_size, ν_size, inf_size, age_size)
    policy_c = zeros(a_size, n_size, ϵ_size, ν_size, inf_size, age_size)
    policy_a_p = zeros(a_size, n_size, ϵ_size, ν_size, inf_size, age_size)
    policy_e = zeros(a_size, n_size, ϵ_size, ν_size, inf_size, age_size)
    policy_K = zeros(a_size, n_size, ϵ_size, ν_size, inf_size, age_size)

    # EGM workspace
    V_endo = Vector{Float64}(undef, a_size)
    a_endo = Vector{Float64}(undef, a_size)
    ap_endo = Vector{Float64}(undef, a_size)
    V1 = Vector{Float64}(undef, a_size)
    c1 = Vector{Float64}(undef, a_size)
    ap1 = Vector{Float64}(undef, a_size)
    e1 = Vector{Float64}(undef, a_size)
    EV_next = Array{Float64}(undef, a_size, n_size, ϵ_size, inf_size)
    c_next_CE = Array{Float64}(undef, a_size, n_size, ϵ_size, inf_size)
    EGM_ws = EGMWorkspace(V_endo, a_endo, ap_endo, V1, c1, ap1, e1, EV_next, c_next_CE)

    # return outputs
    variables = Mutable_Variables(V, policy_c, policy_a_p, policy_e, policy_K, EGM_ws)
    return variables
end

@inline function lininterp(x::AbstractVector{Float64}, y::AbstractVector{Float64}, xq::Float64)
    @inbounds begin
        if xq <= x[1]
            return y[1]
        elseif xq >= x[end]
            return y[end]
        else
            j = searchsortedlast(x, xq)
            x0 = x[j]
            x1 = x[j+1]
            y0 = y[j]
            y1 = y[j+1]
            w = (xq - x0) / (x1 - x0)
            return muladd(w, (y1 - y0), y0)
        end
    end
end

@inline function Vcont_safe(a_grid::AbstractVector{Float64},
    EV_next_a::AbstractVector{Float64},
    ap::Float64,
    V_penalty::Float64)
    if !isfinite(ap)
        return V_penalty
    end

    if ap <= a_grid[1]
        v = EV_next_a[1]
        return (isfinite(v) && v > V_penalty) ? v : V_penalty
    elseif ap >= a_grid[end]
        v = EV_next_a[end]
        return (isfinite(v) && v > V_penalty) ? v : V_penalty
    end

    j = searchsortedlast(a_grid, ap)
    @inbounds begin
        x0 = a_grid[j]
        x1 = a_grid[j+1]
        y0 = EV_next_a[j]
        y1 = EV_next_a[j+1]
    end

    if !(isfinite(y0) && y0 > V_penalty && isfinite(y1) && y1 > V_penalty)
        return V_penalty
    end

    w = (ap - x0) / (x1 - x0)
    return muladd(w, (y1 - y0), y0)
end

@inline function terminal_step!(
    V_end_a::AbstractVector{Float64},
    c_end_a::AbstractVector{Float64},
    parameters::NamedTuple,
    w_bar::Float64,
)
    @unpack w_grid, aR_grid, one_m_γ, inv_one_m_γ = parameters
    @unpack c_floor, V_penalty = parameters
    @inbounds for i in eachindex(aR_grid)
        c = aR_grid[i] + w_bar
        c_end_a[i] = c
        V_end_a[i] = u_CRRA(c, one_m_γ, inv_one_m_γ; c_floor=c_floor, V_penalty=V_penalty)
    end
    return nothing
end

function opt_ap_binding_retired(
    M::Float64,
    EV_next_a::AbstractVector{Float64},
    parameters::NamedTuple;
    ap_lo::Float64,
    ap_hi::Float64,
)
    @unpack a_grid, β, one_m_γ, inv_one_m_γ, c_floor, V_penalty = parameters

    tol = 10 * eps(max(abs(ap_lo), abs(ap_hi), 1.0))

    if ap_hi < ap_lo - tol
        return (false, NaN, NaN, V_penalty)
    elseif ap_hi <= ap_lo + tol
        ap = ap_lo
        c = M - ap
        if !isfinite(c) || c < c_floor
            return (false, ap, c_floor, V_penalty)
        end
        Vp = Vcont_safe(a_grid, EV_next_a, ap, V_penalty)
        val = u_CRRA(c, one_m_γ, inv_one_m_γ; c_floor=c_floor, V_penalty=V_penalty) + β * Vp
        ok = isfinite(val) && (Vp > V_penalty) && (c >= c_floor)
        return (ok, ap, c, ok ? val : V_penalty)
    end

    PEN = 1e30

    function obj(ap::Float64)
        if ap < ap_lo || ap > ap_hi
            return PEN
        end
        c = M - ap
        if !isfinite(c) || c < c_floor
            return PEN
        end
        Vp = Vcont_safe(a_grid, EV_next_a, ap, V_penalty)
        if !isfinite(Vp) || Vp <= V_penalty
            return PEN
        end
        val = u_CRRA(c, one_m_γ, inv_one_m_γ; c_floor=c_floor, V_penalty=V_penalty) + β * Vp
        return -val
    end

    res = Optim.optimize(obj, ap_lo, ap_hi, Brent())
    ap = Optim.minimizer(res)

    c = M - ap
    Vp = Vcont_safe(a_grid, EV_next_a, ap, V_penalty)
    val = u_CRRA(c, one_m_γ, inv_one_m_γ; c_floor=c_floor, V_penalty=V_penalty) + β * Vp

    ok = isfinite(val) && val > V_penalty &&
         isfinite(Vp) && Vp > V_penalty &&
         isfinite(c) && c >= c_floor

    return (ok, ap, c, ok ? val : V_penalty)
end

function retired_step!(
    V_endo::Vector{Float64},
    a_endo::Vector{Float64},
    ap_endo::Vector{Float64},
    EV_next_a::AbstractVector{Float64},
    c_next_a::AbstractVector{Float64},
    V_current_a::AbstractVector{Float64},
    c_current_a::AbstractVector{Float64},
    a_p_current_a::AbstractVector{Float64},
    parameters::NamedTuple,
    w_bar::Float64,
)
    @unpack w_grid, a_size, a_grid, a_min, aR_grid, inv_R, EGM_fac, β, one_m_γ, inv_one_m_γ = parameters
    @unpack c_floor, V_penalty = parameters

    nv = 0
    @inbounds for ap_i in 1:a_size
        Vp = EV_next_a[ap_i]
        if !isfinite(Vp) || Vp <= V_penalty
            continue
        end
        cnext = c_next_a[ap_i]
        if !isfinite(cnext) || cnext < c_floor
            continue
        end
        c = EGM_fac * cnext
        if !isfinite(c) || c < c_floor
            continue
        end
        ap = a_grid[ap_i]
        m = c + ap
        a_today = (m - w_bar) * inv_R
        if !isfinite(a_today)
            continue
        end
        nv += 1
        V_endo[nv] = Vp
        ap_endo[nv] = ap
        a_endo[nv] = a_today
    end

    if nv == 0
        @inbounds for a_i in 1:a_size
            a_p_current_a[a_i] = a_min
            c_current_a[a_i] = c_floor
            V_current_a[a_i] = V_penalty
        end
        return nothing
    end

    V_endo_v = @view V_endo[1:nv]
    a_endo_v = @view a_endo[1:nv]
    ap_endo_v = @view ap_endo[1:nv]

    @inbounds for j in 2:nv
        if a_endo_v[j] <= a_endo_v[j-1]
            a_endo_v[j] = nextfloat(a_endo_v[j-1])
        end
    end

    ibind = searchsortedlast(a_grid, a_endo_v[1])

    if ibind > 0
        @inbounds for a_i in 1:ibind
            aR = aR_grid[a_i]
            M = aR + w_bar
            ap_hi = min(a_grid[end], M - c_floor)
            ap_lo = a_min
            ok, ap, c, Vnow = opt_ap_binding_retired(M, EV_next_a, parameters; ap_lo=ap_lo, ap_hi=ap_hi)
            if !ok
                a_p_current_a[a_i] = a_min
                c_current_a[a_i] = c_floor
                V_current_a[a_i] = V_penalty
            else
                a_p_current_a[a_i] = ap
                c_current_a[a_i] = c
                V_current_a[a_i] = Vnow
            end
        end
    end

    @inbounds for a_i in (ibind+1):a_size
        a = a_grid[a_i]
        aR = aR_grid[a_i]
        ap = lininterp(a_endo_v, ap_endo_v, a)
        c = aR + w_bar - ap
        if !isfinite(c) || c < c_floor
            a_p_current_a[a_i] = ap
            c_current_a[a_i] = c_floor
            V_current_a[a_i] = V_penalty
        else
            a_p_current_a[a_i] = ap
            c_current_a[a_i] = c
            Vp = Vcont_safe(a_grid, EV_next_a, ap, V_penalty)
            if !isfinite(Vp) || Vp <= V_penalty
                c_current_a[a_i] = c_floor
                V_current_a[a_i] = V_penalty
                continue
            end
            V_current_a[a_i] = u_CRRA(c, one_m_γ, inv_one_m_γ; c_floor=c_floor, V_penalty=V_penalty) + β * Vp
        end
    end

    return nothing
end

@inline g_eval(e::Float64, m::Float64, A::Float64, γ::Float64, κ::Float64) = A * e^(-κ) - (m - e)^(-γ)

@inline function solve_e_bisect(
    m::Float64,
    n::Float64,
    P::Float64,
    ψ::Float64,
    γ::Float64,
    κ::Float64,
    one_m_κ::Float64,
    e_bar::Float64;
    c_floor::Float64=1e-12,
    maxit::Int64=80,
    tol::Float64=1e-12,
)::Float64

    # basic guards
    if !(isfinite(m) && isfinite(n) && isfinite(P) && isfinite(ψ) &&
         isfinite(γ) && isfinite(κ) && isfinite(one_m_κ) &&
         isfinite(e_bar) && isfinite(c_floor))
        return NaN
    end
    if P <= 0.0 || n < 0.0 || e_bar <= 0.0 || c_floor <= 0.0
        return NaN
    end

    # feasibility upper bound for e given c >= c_floor
    e_hi = m - c_floor
    if !isfinite(e_hi)
        return NaN
    end

    # tolerance around the kink e_hi ≈ e_bar
    tol_e = 10 * eps(max(abs(m), abs(e_bar), 1.0))

    # infeasible: cannot satisfy e >= e_bar and c >= c_floor simultaneously
    if e_hi < e_bar - tol_e
        return NaN
    end

    # degenerate: only feasible point is e = e_bar (within tolerance)
    if e_hi <= e_bar + tol_e
        return e_bar
    end

    # auxiliary parameter A = ψ * (n/P)^(1-κ)
    ratio = n / P
    if !(isfinite(ratio) && ratio >= 0.0)
        return NaN
    end
    scale = ratio^one_m_κ
    if !isfinite(scale)
        return NaN
    end
    A = ψ * scale
    if !isfinite(A)
        return NaN
    end
    # if marginal benefit is non-positive, corner at e_bar
    if A <= 0.0
        return e_bar
    end

    # bracket checks (MUST use isfinite, not just isnan)
    gl = g_eval(e_bar, m, A, γ, κ)
    if !isfinite(gl)
        return NaN
    end
    if gl <= 0.0
        return e_bar
    end

    gh = g_eval(e_hi, m, A, γ, κ)
    if !isfinite(gh)
        return NaN
    end
    if gh >= 0.0
        return e_hi
    end

    # unique root in (e_bar, e_hi)
    lo = e_bar
    hi = e_hi
    mid = 0.5 * (lo + hi)

    @inbounds for _ in 1:maxit
        mid = 0.5 * (lo + hi)
        gm = g_eval(mid, m, A, γ, κ)
        if !isfinite(gm)
            return NaN
        end

        if abs(gm) <= tol || (hi - lo) <= tol * (1.0 + mid)
            return mid
        end

        if gm > 0.0
            lo = mid
        else
            hi = mid
        end
    end

    return mid
end

function opt_ap_binding_infertile(
    M::Float64,
    n::Float64,
    P::Float64,
    e_bar::Float64,
    EV_next_a::AbstractVector{Float64},
    parameters::NamedTuple;
    ap_lo::Float64,
    ap_hi::Float64,
)
    @unpack a_grid, β, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ = parameters
    @unpack ψ, γ, κ, c_floor, V_penalty = parameters

    scale = (n / P)^one_m_κ
    tol = 10 * eps(max(abs(ap_lo), abs(ap_hi), 1.0))

    if ap_hi < ap_lo - tol
        return (false, NaN, NaN, NaN, V_penalty)
    elseif ap_hi <= ap_lo + tol
        ap = ap_lo
        m = M - ap
        if !(isfinite(m) && m >= e_bar + c_floor)
            return (false, ap, c_floor, e_bar, V_penalty)
        end

        Vp = Vcont_safe(a_grid, EV_next_a, ap, V_penalty)
        if !(isfinite(Vp) && Vp > V_penalty)
            return (false, ap, c_floor, e_bar, V_penalty)
        end

        e = solve_e_bisect(m, n, P, ψ, γ, κ, one_m_κ, e_bar; c_floor=c_floor)
        if !isfinite(e)
            return (false, ap, c_floor, e_bar, V_penalty)
        end
        e = max(e, e_bar)

        c = m - e
        if !(isfinite(c) && c >= c_floor)
            return (false, ap, c_floor, e, V_penalty)
        end

        val = u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ, scale;
            c_floor=c_floor, V_penalty=V_penalty) + β * Vp

        ok = isfinite(val) && val > V_penalty && Vp > V_penalty && c >= c_floor && e >= e_bar
        return (ok, ap, c, e, ok ? val : V_penalty)
    end

    PEN = 1e30

    function obj(ap::Float64)
        if ap < ap_lo || ap > ap_hi
            return PEN
        end

        m = M - ap
        if !isfinite(m) || m < e_bar + c_floor
            return PEN
        end

        Vp = Vcont_safe(a_grid, EV_next_a, ap, V_penalty)
        if !isfinite(Vp) || Vp <= V_penalty
            return PEN
        end

        e = solve_e_bisect(m, n, P, ψ, γ, κ, one_m_κ, e_bar; c_floor=c_floor)
        if !isfinite(e)
            return PEN
        end
        e = max(e, e_bar)

        c = m - e
        if !isfinite(c) || c < c_floor
            return PEN
        end

        val = u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ, scale;
            c_floor=c_floor, V_penalty=V_penalty) + β * Vp

        return (isfinite(val) ? -val : PEN)
    end

    res = Optim.optimize(obj, ap_lo, ap_hi, Brent())
    ap = Optim.minimizer(res)

    m = M - ap
    Vp = Vcont_safe(a_grid, EV_next_a, ap, V_penalty)
    e = solve_e_bisect(m, n, P, ψ, γ, κ, one_m_κ, e_bar; c_floor=c_floor)
    if !isfinite(e)
        return (false, ap, c_floor, e_bar, V_penalty)
    end
    e = max(e, e_bar)

    c = m - e
    if !(isfinite(c) && c >= c_floor && isfinite(Vp) && Vp > V_penalty && isfinite(m) && m >= e_bar + c_floor)
        return (false, ap, c_floor, e, V_penalty)
    end

    val = u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ, scale;
        c_floor=c_floor, V_penalty=V_penalty) + β * Vp

    ok = isfinite(val) && val > V_penalty
    return (ok, ap, c, e, ok ? val : V_penalty)
end

function infertile_step!(
    V_endo::Vector{Float64},
    a_endo::Vector{Float64},
    ap_endo::Vector{Float64},
    EV_next_a::AbstractVector{Float64},
    c_next_a::AbstractVector{Float64},
    V_current_a::AbstractVector{Float64},
    c_current_a::AbstractVector{Float64},
    a_p_current_a::AbstractVector{Float64},
    e_current_a::AbstractVector{Float64},
    parameters::NamedTuple,
    w_bar::Float64,
    P::Float64,
    e_bar::Float64,
    n::Float64,
)

    @unpack w_grid, a_size, a_grid, a_min, aR_grid, inv_R, EGM_fac, inv_κ, γ_by_κ, P_grid, q_bar_P_grid = parameters
    @unpack β, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ = parameters
    @unpack n_grid, ψ, γ, κ = parameters
    @unpack c_floor, V_penalty = parameters

    scale = (n / P)^one_m_κ
    e_para = (ψ * scale)^inv_κ

    if !isfinite(e_para)
        @inbounds for a_i in 1:a_size
            V_current_a[a_i] = V_penalty
            c_current_a[a_i] = c_floor
            a_p_current_a[a_i] = a_min
            e_current_a[a_i] = e_bar
        end
        return nothing
    end

    nv = 0
    @inbounds for ap_i in 1:a_size
        Vp = EV_next_a[ap_i]
        if !isfinite(Vp) || Vp <= V_penalty
            continue
        end
        cnext = c_next_a[ap_i]
        if !isfinite(cnext) || cnext < c_floor
            continue
        end
        c = EGM_fac * cnext
        if !isfinite(c) || c < c_floor
            continue
        end
        e_foc = e_para * c^γ_by_κ
        e = max(e_foc, e_bar)
        if !isfinite(e)
            continue
        end
        ap = a_grid[ap_i]
        m = c + ap + e
        a_today = (m - w_bar) * inv_R
        if !isfinite(a_today)
            continue
        end
        nv += 1
        V_endo[nv] = Vp
        ap_endo[nv] = ap
        a_endo[nv] = a_today
    end

    if nv == 0
        @inbounds for a_i in 1:a_size
            V_current_a[a_i] = V_penalty
            c_current_a[a_i] = c_floor
            a_p_current_a[a_i] = a_min
            e_current_a[a_i] = e_bar
        end
        return nothing
    end

    V_endo_v = @view V_endo[1:nv]
    a_endo_v = @view a_endo[1:nv]
    ap_endo_v = @view ap_endo[1:nv]

    @inbounds for j in 2:nv
        if a_endo_v[j] <= a_endo_v[j-1]
            a_endo_v[j] = nextfloat(a_endo_v[j-1])
        end
    end

    ibind = searchsortedlast(a_grid, a_endo_v[1])

    if ibind > 0
        @inbounds for a_i in 1:ibind
            aR = aR_grid[a_i]
            M = aR + w_bar
            ap_hi = min(a_grid[end], M - (e_bar + c_floor))
            ap_lo = a_min
            ok, ap, c, e, Vnow = opt_ap_binding_infertile(M, n, P, e_bar, EV_next_a, parameters; ap_lo=ap_lo, ap_hi=ap_hi)

            if !ok
                V_current_a[a_i] = V_penalty
                c_current_a[a_i] = c_floor
                a_p_current_a[a_i] = a_min
                e_current_a[a_i] = e_bar
            else
                V_current_a[a_i] = Vnow
                c_current_a[a_i] = c
                a_p_current_a[a_i] = ap
                e_current_a[a_i] = e
            end
        end
    end

    @inbounds for a_i in (ibind+1):a_size
        a = a_grid[a_i]
        aR = aR_grid[a_i]
        ap = lininterp(a_endo_v, ap_endo_v, a)
        M = aR + w_bar
        m = M - ap
        if m < e_bar + c_floor
            V_current_a[a_i] = V_penalty
            c_current_a[a_i] = c_floor
            a_p_current_a[a_i] = ap
            e_current_a[a_i] = e_bar
        else
            Vp = Vcont_safe(a_grid, EV_next_a, ap, V_penalty)
            if !isfinite(Vp) || Vp <= V_penalty
                V_current_a[a_i] = V_penalty
                c_current_a[a_i] = c_floor
                a_p_current_a[a_i] = ap
                e_current_a[a_i] = e_bar
                continue
            end
            e = solve_e_bisect(m, n, P, ψ, γ, κ, one_m_κ, e_bar; c_floor=c_floor)
            if !isfinite(e) || e < e_bar
                V_current_a[a_i] = V_penalty
                c_current_a[a_i] = c_floor
                a_p_current_a[a_i] = ap
                e_current_a[a_i] = e_bar
                continue
            end
            c = m - e
            V_current_a[a_i] = u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ, scale; c_floor=c_floor, V_penalty=V_penalty) + β * Vp
            c_current_a[a_i] = c
            a_p_current_a[a_i] = ap
            e_current_a[a_i] = e
        end
    end
    return nothing
end

function fertile_step!(
    V_endo::Vector{Float64},
    a_endo::Vector{Float64},
    ap_endo::Vector{Float64},
    V_next_a::AbstractVector{Float64},
    c_next_a::AbstractVector{Float64},
    V_next_aK::AbstractVector{Float64},
    c_next_aK::AbstractVector{Float64},
    V_current_a::AbstractVector{Float64},
    c_current_a::AbstractVector{Float64},
    a_p_current_a::AbstractVector{Float64},
    e_current_a::AbstractVector{Float64},
    K_current_a::AbstractVector{Float64},
    V1::Vector{Float64},
    c1::Vector{Float64},
    ap1::Vector{Float64},
    e1::Vector{Float64},
    parameters::NamedTuple,
    w_bar::Float64,
    P::Float64,
    e_bar::Float64,
    n::Float64,
)
    fill!(K_current_a, 0.0)

    if n == 0.0
        retired_step!(V_endo, a_endo, ap_endo,
            V_next_a, c_next_a,
            V_current_a, c_current_a, a_p_current_a,
            parameters, w_bar)
        retired_step!(V_endo, a_endo, ap_endo,
            V_next_aK, c_next_aK,
            V1, c1, ap1,
            parameters, w_bar)
        @inbounds for i in eachindex(V_current_a)
            if V1[i] > V_current_a[i]
                V_current_a[i] = V1[i]
                c_current_a[i] = c1[i]
                a_p_current_a[i] = ap1[i]
                K_current_a[i] = 1.0
            end
        end
    else
        infertile_step!(V_endo, a_endo, ap_endo,
            V_next_a, c_next_a,
            V_current_a, c_current_a, a_p_current_a, e_current_a,
            parameters, w_bar, P, e_bar, n)
        infertile_step!(V_endo, a_endo, ap_endo,
            V_next_aK, c_next_aK,
            V1, c1, ap1, e1,
            parameters, w_bar, P, e_bar, n)
        @inbounds for i in eachindex(V_current_a)
            if V1[i] > V_current_a[i]
                V_current_a[i] = V1[i]
                c_current_a[i] = c1[i]
                a_p_current_a[i] = ap1[i]
                e_current_a[i] = e1[i]
                K_current_a[i] = 1.0
            end
        end
    end
    return nothing
end

function fill_EV_Euc!(
    EV_next::AbstractArray{Float64,4},
    c_next_CE::AbstractArray{Float64,4},
    V_next::AbstractArray{Float64,N},
    policy_c_next::AbstractArray{Float64,N},
    parameters::NamedTuple,
    age_i::Int64,
) where {N}

    @unpack a_size, n_size, ϵ_size, inf_size, Γ_ret, Γ_inf, Γ, γ, m_inv_γ = parameters
    @unpack c_floor, V_penalty = parameters

    if N == 3
        @inbounds for a_p_i in 1:a_size, n_i in 1:n_size, ϵ_i in 1:ϵ_size
            @views Γt = Γ_ret[:, :, ϵ_i]
            @views Vt = V_next[a_p_i, :, :]
            @views ct = policy_c_next[a_p_i, :, :]
            EV = 0.0
            Euc = 0.0
            bad = false
            for k in eachindex(Γt)
                p = Γt[k]
                if p == 0.0
                    continue
                end
                Vp = Vt[k]
                c = ct[k]
                if !isfinite(Vp) || Vp <= V_penalty || !isfinite(c) || c <= c_floor
                    bad = true
                    break
                end
                EV += p * Vp
                Euc += p * (c^(-γ))
            end
            if bad
                EV_next[a_p_i, n_i, ϵ_i, 2] = V_penalty
                c_next_CE[a_p_i, n_i, ϵ_i, 2] = NaN
            else
                EV_next[a_p_i, n_i, ϵ_i, 2] = EV
                c_next_CE[a_p_i, n_i, ϵ_i, 2] = Euc^m_inv_γ
            end
        end
    elseif N == 4
        @inbounds for a_p_i in 1:a_size, n_i in 1:n_size, ϵ_i in 1:ϵ_size
            @views Γt = Γ_inf[:, :, :, n_i, ϵ_i]
            @views Vt = V_next[a_p_i, :, :, :]
            @views ct = policy_c_next[a_p_i, :, :, :]
            EV = 0.0
            Euc = 0.0
            bad = false
            for k in eachindex(Γt)
                p = Γt[k]
                if p == 0.0
                    continue
                end
                Vp = Vt[k]
                c = ct[k]
                if !isfinite(Vp) || Vp <= V_penalty || !isfinite(c) || c <= c_floor
                    bad = true
                    break
                end
                EV += p * Vp
                Euc += p * (c^(-γ))
            end
            if bad
                EV_next[a_p_i, n_i, ϵ_i, 2] = V_penalty
                c_next_CE[a_p_i, n_i, ϵ_i, 2] = NaN
            else
                EV_next[a_p_i, n_i, ϵ_i, 2] = EV
                c_next_CE[a_p_i, n_i, ϵ_i, 2] = Euc^m_inv_γ
            end
        end
    else
        @inbounds for a_p_i in 1:a_size, n_i in 1:n_size, ϵ_i in 1:ϵ_size, inf_i in 1:inf_size
            @views Γt = Γ[:, :, :, :, n_i, ϵ_i, inf_i, age_i]
            @views Vt = V_next[a_p_i, :, :, :, :]
            @views ct = policy_c_next[a_p_i, :, :, :, :]
            EV = 0.0
            Euc = 0.0
            bad = false
            for k in eachindex(Γt)
                p = Γt[k]
                if p == 0.0
                    continue
                end
                Vp = Vt[k]
                c = ct[k]
                if !isfinite(Vp) || Vp <= V_penalty || !isfinite(c) || c < c_floor
                    bad = true
                    break
                end
                EV += p * Vp
                Euc += p * (c^(-γ))
            end
            if bad
                EV_next[a_p_i, n_i, ϵ_i, inf_i] = V_penalty
                c_next_CE[a_p_i, n_i, ϵ_i, inf_i] = NaN
            else
                EV_next[a_p_i, n_i, ϵ_i, inf_i] = EV
                c_next_CE[a_p_i, n_i, ϵ_i, inf_i] = Euc^m_inv_γ
            end
        end
    end
    return nothing
end

function solve_value_and_policy_function!(variables::Mutable_Variables, parameters::NamedTuple)
    """
    Compute value and policy functions
    """

    # unpack parameters
    @unpack age_size, age_grid, age_max, age_min, age_ret, age_inf = parameters
    @unpack ν_size, ν_grid, ν_Γ = parameters
    @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
    @unpack n_size, n_grid, n_Γ, n_max = parameters
    @unpack a_size, a_grid, a_min = parameters
    @unpack aR_grid, w_grid, one_m_γ, inv_one_m_γ, EGM_fac, P_grid, q_bar_P_grid = parameters
    @unpack inf_size, inf_grid = parameters
    @unpack h_grid = parameters
    @unpack b, r, R, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

    # extract working variables
    V_endo = variables.EGM_ws.V_endo
    a_endo = variables.EGM_ws.a_endo
    ap_endo = variables.EGM_ws.ap_endo
    V1 = variables.EGM_ws.V1
    c1 = variables.EGM_ws.c1
    ap1 = variables.EGM_ws.ap1
    e1 = variables.EGM_ws.e1
    EV_next = variables.EGM_ws.EV_next
    c_next_CE = variables.EGM_ws.c_next_CE

    println("Solving the HH problem in the last period (at age $age_max)")
    @views V_end = variables.V[:, 1, :, :, 2, end]
    @views c_end = variables.policy_c[:, 1, :, :, 2, end]
    for ϵ_i in 1:ϵ_size, ν_i in 1:ν_size
        @inbounds w_bar = w_grid[ϵ_i, ν_i, end]
        @views V_end_a = V_end[:, ϵ_i, ν_i]
        @views c_end_a = c_end[:, ϵ_i, ν_i]
        terminal_step!(V_end_a, c_end_a, parameters, w_bar)
    end

    println("Solving the HH problem after retirement until the second last period (from age $age_ret to $(age_max-1))")
    for age_i in (age_size-1):(-1):(age_ret-age_min+1)
        @views V_next = variables.V[:, 1, :, :, 2, age_i+1]
        @views c_next = variables.policy_c[:, 1, :, :, 2, age_i+1]
        @views V_current = variables.V[:, 1, :, :, 2, age_i]
        @views c_current = variables.policy_c[:, 1, :, :, 2, age_i]
        @views a_p_current = variables.policy_a_p[:, 1, :, :, 2, age_i]
        for ν_i in 1:ν_size, ϵ_i in 1:ϵ_size
            @inbounds w_bar = w_grid[ϵ_i, ν_i, age_i]
            @views V_next_a = V_next[:, ϵ_i, ν_i]
            @views c_next_a = c_next[:, ϵ_i, ν_i]
            @views V_current_a = V_current[:, ϵ_i, ν_i]
            @views c_current_a = c_current[:, ϵ_i, ν_i]
            @views a_p_current_a = a_p_current[:, ϵ_i, ν_i]
            retired_step!(V_endo, a_endo, ap_endo, V_next_a, c_next_a, V_current_a, c_current_a, a_p_current_a, parameters, w_bar)
        end
    end

    println("Solving the HH problem just before retirement (at age $(age_ret-1))")
    age_i = age_ret - age_min
    V_next = @view variables.V[:, 1, :, :, 2, age_i+1]
    c_next = @view variables.policy_c[:, 1, :, :, 2, age_i+1]
    fill_EV_Euc!(EV_next, c_next_CE, V_next, c_next, parameters, age_i)
    V_current = @view variables.V[:, :, :, :, 2, age_i]
    c_current = @view variables.policy_c[:, :, :, :, 2, age_i]
    a_p_current = @view variables.policy_a_p[:, :, :, :, 2, age_i]
    e_current = @view variables.policy_e[:, :, :, :, 2, age_i]
    for ϵ_i in 1:ϵ_size, ν_i in 1:ν_size
        @inbounds w_bar = w_grid[ϵ_i, ν_i, age_i]
        for n_i in 1:n_size
            @inbounds P = P_grid[n_i, ϵ_i, ν_i, age_i]
            @inbounds e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
            @inbounds n = n_grid[n_i]
            EV_next_a = @view EV_next[:, n_i, ϵ_i, 2]
            c_next_CE_a = @view c_next_CE[:, n_i, ϵ_i, 2]
            V_current_a = @view V_current[:, n_i, ϵ_i, ν_i]
            c_current_a = @view c_current[:, n_i, ϵ_i, ν_i]
            a_p_current_a = @view a_p_current[:, n_i, ϵ_i, ν_i]
            e_current_a = @view e_current[:, n_i, ϵ_i, ν_i]
            if n_i == 1
                retired_step!(V_endo, a_endo, ap_endo, EV_next_a, c_next_CE_a, V_current_a, c_current_a, a_p_current_a, parameters, w_bar)
            else
                infertile_step!(V_endo, a_endo, ap_endo, EV_next_a, c_next_CE_a, V_current_a, c_current_a, a_p_current_a, e_current_a, parameters, w_bar, P, e_bar, n)
            end
        end
    end

    println("Solving the HH problem from the infertile age to before retirement (from age $age_inf to $(age_ret-2))")
    for age_i in (age_ret-age_min-1):(-1):(age_inf-age_min+1)
        V_next = @view variables.V[:, :, :, :, 2, age_i+1]
        c_next = @view variables.policy_c[:, :, :, :, 2, age_i+1]
        fill_EV_Euc!(EV_next, c_next_CE, V_next, c_next, parameters, age_i)
        V_current = @view variables.V[:, :, :, :, 2, age_i]
        c_current = @view variables.policy_c[:, :, :, :, 2, age_i]
        a_p_current = @view variables.policy_a_p[:, :, :, :, 2, age_i]
        e_current = @view variables.policy_e[:, :, :, :, 2, age_i]
        for ν_i in 1:ν_size, ϵ_i in 1:ϵ_size
            @inbounds w_bar = w_grid[ϵ_i, ν_i, age_i]
            for n_i in 1:n_size
                @inbounds P = P_grid[n_i, ϵ_i, ν_i, age_i]
                @inbounds e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
                @inbounds n = n_grid[n_i]
                EV_next_a = @view EV_next[:, n_i, ϵ_i, 2]
                c_next_CE_a = @view c_next_CE[:, n_i, ϵ_i, 2]
                V_current_a = @view V_current[:, n_i, ϵ_i, ν_i]
                c_current_a = @view c_current[:, n_i, ϵ_i, ν_i]
                a_p_current_a = @view a_p_current[:, n_i, ϵ_i, ν_i]
                e_current_a = @view e_current[:, n_i, ϵ_i, ν_i]
                if n_i == 1
                    retired_step!(V_endo, a_endo, ap_endo, EV_next_a, c_next_CE_a, V_current_a, c_current_a, a_p_current_a, parameters, w_bar)
                else
                    infertile_step!(V_endo, a_endo, ap_endo, EV_next_a, c_next_CE_a, V_current_a, c_current_a, a_p_current_a, e_current_a, parameters, w_bar, P, e_bar, n)
                end
            end
        end
    end

    println("Solving the HH problem just before the infertile age (at age $(age_inf-1))")
    age_i = age_inf - age_min
    V_next = @view variables.V[:, :, :, :, 2, age_i+1]
    c_next = @view variables.policy_c[:, :, :, :, 2, age_i+1]
    fill_EV_Euc!(EV_next, c_next_CE, V_next, c_next, parameters, age_i)
    V_current = @view variables.V[:, :, :, :, :, age_i]
    c_current = @view variables.policy_c[:, :, :, :, :, age_i]
    a_p_current = @view variables.policy_a_p[:, :, :, :, :, age_i]
    e_current = @view variables.policy_e[:, :, :, :, :, age_i]
    K_current = @view variables.policy_K[:, :, :, :, :, age_i]
    for ϵ_i in 1:ϵ_size, ν_i in 1:ν_size
        @inbounds w_bar = w_grid[ϵ_i, ν_i, age_i]
        for n_i in 1:n_size
            @inbounds P = P_grid[n_i, ϵ_i, ν_i, age_i]
            @inbounds e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
            @inbounds n = n_grid[n_i]
            for inf_i in 1:inf_size
                EV_next_a = @view EV_next[:, n_i, ϵ_i, 2]
                c_next_CE_a = @view c_next_CE[:, n_i, ϵ_i, 2]
                V_current_a = @view V_current[:, n_i, ϵ_i, ν_i, inf_i]
                c_current_a = @view c_current[:, n_i, ϵ_i, ν_i, inf_i]
                a_p_current_a = @view a_p_current[:, n_i, ϵ_i, ν_i, inf_i]
                e_current_a = @view e_current[:, n_i, ϵ_i, ν_i, inf_i]
                K_current_a = @view K_current[:, n_i, ϵ_i, ν_i, inf_i]
                if inf_i == 1
                    if n == n_max
                        infertile_step!(V_endo, a_endo, ap_endo, EV_next_a, c_next_CE_a, V_current_a, c_current_a, a_p_current_a, e_current_a, parameters, w_bar, P, e_bar, n)
                    else
                        EV_next_aK = @view EV_next[:, n_i+1, ϵ_i, 2]
                        c_next_CE_aK = @view c_next_CE[:, n_i+1, ϵ_i, 2]
                        fertile_step!(V_endo, a_endo, ap_endo,
                            EV_next_a, c_next_CE_a, EV_next_aK, c_next_CE_aK,
                            V_current_a, c_current_a, a_p_current_a, e_current_a, K_current_a,
                            V1, c1, ap1, e1,
                            parameters, w_bar, P, e_bar, n)
                    end
                else
                    if n == 0.0
                        retired_step!(V_endo, a_endo, ap_endo, EV_next_a, c_next_CE_a, V_current_a, c_current_a, a_p_current_a, parameters, w_bar)
                    else
                        infertile_step!(V_endo, a_endo, ap_endo, EV_next_a, c_next_CE_a, V_current_a, c_current_a, a_p_current_a, e_current_a, parameters, w_bar, P, e_bar, n)
                    end
                end
            end
        end
    end

    println("Solving the HH problem until two periods before the infertile age (from age $age_min to $(age_inf-2))")
    for age_i in (age_inf-age_min-1):(-1):1
        V_next = @view variables.V[:, :, :, :, :, age_i+1]
        c_next = @view variables.policy_c[:, :, :, :, :, age_i+1]
        fill_EV_Euc!(EV_next, c_next_CE, V_next, c_next, parameters, age_i)
        V_current = @view variables.V[:, :, :, :, :, age_i]
        c_current = @view variables.policy_c[:, :, :, :, :, age_i]
        a_p_current = @view variables.policy_a_p[:, :, :, :, :, age_i]
        e_current = @view variables.policy_e[:, :, :, :, :, age_i]
        K_current = @view variables.policy_K[:, :, :, :, :, age_i]
        for ϵ_i in 1:ϵ_size, ν_i in 1:ν_size
            @inbounds w_bar = w_grid[ϵ_i, ν_i, age_i]
            for n_i in 1:n_size
                @inbounds P = P_grid[n_i, ϵ_i, ν_i, age_i]
                @inbounds e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
                @inbounds n = n_grid[n_i]
                for inf_i in 1:inf_size
                    EV_next_a = @view EV_next[:, n_i, ϵ_i, inf_i]
                    c_next_CE_a = @view c_next_CE[:, n_i, ϵ_i, inf_i]
                    V_current_a = @view V_current[:, n_i, ϵ_i, ν_i, inf_i]
                    c_current_a = @view c_current[:, n_i, ϵ_i, ν_i, inf_i]
                    a_p_current_a = @view a_p_current[:, n_i, ϵ_i, ν_i, inf_i]
                    e_current_a = @view e_current[:, n_i, ϵ_i, ν_i, inf_i]
                    K_current_a = @view K_current[:, n_i, ϵ_i, ν_i, inf_i]
                    if inf_i == 1
                        if n == n_max
                            infertile_step!(V_endo, a_endo, ap_endo, EV_next_a, c_next_CE_a, V_current_a, c_current_a, a_p_current_a, e_current_a, parameters, w_bar, P, e_bar, n)
                        else
                            EV_next_aK = @view EV_next[:, n_i+1, ϵ_i, inf_i]
                            c_next_CE_aK = @view c_next_CE[:, n_i+1, ϵ_i, inf_i]
                            fertile_step!(V_endo, a_endo, ap_endo,
                                EV_next_a, c_next_CE_a, EV_next_aK, c_next_CE_aK,
                                V_current_a, c_current_a, a_p_current_a, e_current_a, K_current_a,
                                V1, c1, ap1, e1,
                                parameters, w_bar, P, e_bar, n)
                        end
                    else
                        if n == 0.0
                            retired_step!(V_endo, a_endo, ap_endo, EV_next_a, c_next_CE_a, V_current_a, c_current_a, a_p_current_a, parameters, w_bar)
                        else
                            infertile_step!(V_endo, a_endo, ap_endo, EV_next_a, c_next_CE_a, V_current_a, c_current_a, a_p_current_a, e_current_a, parameters, w_bar, P, e_bar, n)
                        end
                    end
                end
            end
        end
    end

    return nothing
end

function save_JLD_function!(variables::Mutable_Variables, parameters::NamedTuple; filename::String)
    V = variables.V
    policy_a_p = variables.policy_a_p
    # policy_x = variables.policy_x
    # policy_l = variables.policy_l
    policy_c = variables.policy_c
    policy_e = variables.policy_e
    policy_K = variables.policy_K
    @save filename parameters V policy_a_p policy_c policy_e policy_K
    return nothing
end

# solve stationary equilibrium #

params_NC = parameters_function(edu_ind=:NC)
vars_NC = variables_function(params_NC)
solve_value_and_policy_function!(vars_NC, params_NC)
# save_JLD_function!(vars_NC, params_NC, filename="workspace_benchmark_C.jld2")

params_C = parameters_function(edu_ind=:C)
vars_C = variables_function(params_C)
solve_value_and_policy_function!(vars_C, params_C)
# save_JLD_function!(vars_C, params_C, filename="workspace_benchmark_NC.jld2")

# simuation #

using Random
using Random123
using QuantEcon: Categorical
using Polyester

@inline function make_thread_rngs(; seed::UInt64=UInt64(1124))
    nt = Threads.nthreads()
    key = (seed, UInt64(0))
    return [Philox4x(UInt64, key) for _ in 1:nt]
end

@inline base_counter(h_id::UInt64, t_id::UInt64) = (h_id << 44) | (t_id << 24)

@inline function ap_policy_to_index(
    rng,
    a_grid::AbstractVector{Float64},
    ap::Float64,
    a_size::Int,
)::Int
    if ap ≤ a_grid[1]
        return 1
    elseif ap ≥ a_grid[end]
        return a_size
    else
        j = searchsortedlast(a_grid, ap)
        a0 = a_grid[j]
        a1 = a_grid[j+1]
        w = (ap - a0) / (a1 - a0)
        return (rand(rng) < w) ? (j + 1) : j
    end
end

struct SimulatedPanel
    a_state::Matrix{Int}
    n_state::Matrix{Int}
    ϵ_state::Matrix{Int}
    ν_state::Matrix{Int}
    f_state::Matrix{Int}
    ap_choice::Matrix{Int}
    c_choice::Matrix{Float64}
    e_choice::Matrix{Float64}
    K_choice::Matrix{Int8}
    ΔK_choice::Matrix{Int8}
end

function initialize_panel(; num_households::Int, num_periods::Int)
    T, N = num_periods, num_households
    return SimulatedPanel(
        zeros(Int, T, N),
        zeros(Int, T, N),
        zeros(Int, T, N),
        zeros(Int, T, N),
        zeros(Int, T, N),
        zeros(Int, T, N),
        zeros(Float64, T, N),
        zeros(Float64, T, N),
        zeros(Int8, T, N),
        zeros(Int8, T, N),
    )
end

function simulate_household_panel!(
    panel::SimulatedPanel,
    variables::Mutable_Variables,
    parameters::NamedTuple;
    seed::UInt64=UInt64(1124),
)
    @unpack a_grid, a_ind_zero, a_size, ϵ_size, n_size, n_Γ, ϵ_Γ, ϵ_G, ν_Γ, inf_grid, age_ret, age_inf, age_min, p_fertile_at_entry = parameters

    T, N = size(panel.a_state)

    ϵ_init = Categorical(ϵ_G)
    ν_iid = Categorical(ν_Γ)
    ϵ_cat = [Categorical(ϵ_Γ[i, :]) for i in 1:ϵ_size]
    n_cat = [Categorical(n_Γ[i, :]) for i in 1:n_size]

    rngs = make_thread_rngs(seed=seed)

    @batch for h in 1:N

        rng = rngs[Threads.threadid()]
        h_id = UInt64(h)

        set_counter!(rng, base_counter(h_id, UInt64(0)))

        panel.a_state[1, h] = a_ind_zero
        panel.n_state[1, h] = 1
        panel.f_state[1, h] = (rand(rng) < p_fertile_at_entry) ? 1 : 2
        panel.ϵ_state[1, h] = rand(rng, ϵ_init)
        panel.ν_state[1, h] = rand(rng, ν_iid)

        @inbounds for t in 1:(T-1)

            set_counter!(rng, base_counter(h_id, UInt64(t)))

            a_i = panel.a_state[t, h]
            n_i = panel.n_state[t, h]
            ϵ_i = panel.ϵ_state[t, h]
            ν_i = panel.ν_state[t, h]
            f_i = panel.f_state[t, h]

            c_raw = variables.policy_c[a_i, n_i, ϵ_i, ν_i, f_i, t]
            e_raw = variables.policy_e[a_i, n_i, ϵ_i, ν_i, f_i, t]
            panel.c_choice[t, h] = c_raw
            panel.e_choice[t, h] = e_raw

            K_raw = variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, t]
            K01 = (K_raw >= 0.5) ? Int8(1) : Int8(0)
            panel.K_choice[t, h] = K01

            Δa_i = a_i != a_size ? a_i + 1 : a_size
            ΔK_raw = variables.policy_K[Δa_i, n_i, ϵ_i, ν_i, f_i, t]
            ΔK01 = (ΔK_raw >= 0.5) ? Int8(1) : Int8(0)
            panel.ΔK_choice[t, h] = ΔK01

            ap_raw = variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, t]
            ap_idx = ap_policy_to_index(rng, a_grid, ap_raw, a_size)
            panel.ap_choice[t, h] = ap_idx
            panel.a_state[t+1, h] = ap_idx

            if t < age_ret - age_min
                n_eff = min(n_size, n_i + Int(K01))
                panel.n_state[t+1, h] = rand(rng, n_cat[n_eff])
            else
                panel.n_state[t+1, h] = 1
            end

            if t <= age_ret - age_min
                panel.ϵ_state[t+1, h] = rand(rng, ϵ_cat[ϵ_i])
                panel.ν_state[t+1, h] = rand(rng, ν_iid)
            else
                panel.ϵ_state[t+1, h] = ϵ_i
                panel.ν_state[t+1, h] = ν_i
            end

            if (f_i == 1) && (t <= (age_inf - age_min))
                panel.f_state[t+1, h] = (rand(rng) < inf_grid[t]) ? 2 : 1
            else
                panel.f_state[t+1, h] = 2
            end
        end
    end
    return nothing
end

panel_NC = initialize_panel(num_households=300_000, num_periods=params_NC.age_size)
simulate_household_panel!(panel_NC, vars_NC, params_NC; seed=UInt64(1))

panel_C = initialize_panel(num_households=300_000, num_periods=params_C.age_size)
simulate_household_panel!(panel_C, vars_C, params_C; seed=UInt64(2))

mean_K_NC = vec(mean(panel_NC.K_choice, dims=2))
mean_K_C  = vec(mean(panel_C.K_choice, dims=2))
plot(params_NC.age_grid, mean_K_NC, label="NC", lw=2)
plot!(params_C.age_grid, mean_K_C, label="C", lw=2,
    xlabel="Age", ylabel="Cumulative fertility")

