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

function adda_cooper(N::Integer, ρ::Real, σ::Real; μ::Real=0.0)
    """
    Approximation of an autoregression process with a Markov chain proposed by Adda and Cooper (2003)
    """
    σ_ϵ = σ / sqrt(1.0 - ρ^2.0)
    q = quantile(Normal(), range(0, 1; length=N + 1))
    ϵ = σ_ϵ .* q .+ μ
    z = [N * σ_ϵ * (pdf(Normal(), (ϵ[i] - μ) / σ_ϵ) - pdf(Normal(), (ϵ[i+1] - μ) / σ_ϵ)) + μ for i = 1:N]

    Π = zeros(N, N)
    dist_σ = Normal(μ, σ_ϵ)
    if ρ == 0.0
        Π .= 1.0 / N
    else
        for i = 1:N, j = 1:N
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
    for i = 1:n_max
        d = Binomial(i, 1.0 - p)
        for j = 0:i
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
    model_inf[findall(model_age .== (age_inf + 1))[]:end] .= 1.0
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
    age_inf = model_age[findlast(model_inf .< 1.0)[]]
    model_inf[findall(model_age .== (age_inf + 1))[]:end] .= 1.0
    return model_age, model_inf, age_inf
end

function h_function(data_h::Array{Float64,1}, age_min::Integer, age_ret::Integer)
    """
    Curve fit of life-cycle component, intrapolated on ages up to age_ret 
    """
    model(t, ω) = ω[1] .+ ω[2] * t .+ ω[3] * t .^ 2 .+ ω[4] * t .^ 3
    ω_int = [0.5, 0.5, 0.5, 0.5]
    model_age = collect(age_min:age_ret)
    fit = curve_fit(model, model_age, data_h, ω_int)
    model_h = fit.param[1] .+ fit.param[2] * model_age .+ fit.param[3] * model_age .^ 2 .+ fit.param[4] * model_age .^ 3
    return model_h
end

function u_CRRA_e(c::Float64, e::Float64, one_m_γ::Float64, inv_one_m_γ::Float64, one_m_κ::Float64, ψ_inv_one_m_κ::Float64)
    c <= 0.0 && return -1.0e16
    return c^one_m_γ * inv_one_m_γ + e^one_m_κ * ψ_inv_one_m_κ
end

@inline function u_CRRA(c::Float64, one_m_γ::Float64, inv_one_m_γ::Float64)
    c <= 0.0 && return -1.0e16
    return c^one_m_γ * inv_one_m_γ
end

@inline function u_log(c::Float64)
    c <= 0.0 && return -1.0e16
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
    κ::Real=3.50,                   # preference curvature
    ψ::Real=0.14,                   # preference scale
    μ::Real=0.35,                   # production share
    θ::Real=0.70,                   # elasticity of substitution in production
    q_bar::Real=0.34,               # lower bound on children's consumption
    ψ_1::Real=0.91,                 # HH economies to money input to production
    ψ_2::Real=0.54,                 # HH economies to time input to production
    p::Real=0.02,                   # prob that a child becomes independent
    #====================#
    # numerical solution #
    #====================#
    age_min::Integer=18,            # min age
    age_max::Integer=80,            # max age
    age_edu::Integer=22,            # education age
    age_inf::Integer=45,            # infertile age
    age_ret::Integer=65,            # retirement age
    n_max::Integer=4,               # max number of kids
    ϵ_size::Integer=7,              # number of persistent shock
    ν_size::Integer=5,              # number of transitory shock
    a_max::Real=120,                # max of asset holding
    a_size_neg::Integer=11,         # number of negative asset
    a_size::Integer=50,             # number of asset
    a_degree::Integer=2,            # curvature of asset gridpoints
    q_x::Real=1.0,                  # price of monetary input $x$
    #=================#
    # case indicators #
    #=================#
    inf_scale::Real=1.0,            # scale of infertility risk
    edu_h_ind::Real=0.0,            # education indicator
    edu_h_scale::Real=1.0,          # scale of edu-dependent life-cycle income
    edu_σ_ϵ_scale::Real=1.0,        # scale of persistent income uncertainty    
    edu_σ_ν_scale::Real=1.0,        # scale of transitory income uncertainty
    lr_exp::Integer=0               # long-run response experiment
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
    EGM_fac = (β * R)^m_inv_γ
    one_m_κ = 1.0 - κ
    inv_κ = 1.0 / κ
    inv_one_m_κ = 1.0 / (1.0 - κ)
    ψ_inv_one_m_κ = ψ / (1.0 - κ)
    γ_by_κ = γ / κ

    # infertility parameters: taken from Trussell and Wilson (1985, Population Studies)
    data_inf = [0.07, 0.131, 0.231, 0.345, 0.576, 0.952]
    data_age = [20, 25, 30, 35, 40, 45]
    if inf_scale == 1.0
        age_grid, inf_grid = infertility_risk_function(data_age, data_inf, age_min, age_max, age_inf)
    else
        data_inf_new = data_inf * inf_scale
        age_grid, inf_grid, age_inf = infertility_risk_low_function(data_age, data_inf_new, age_min, age_max)
    end
    age_size = length(age_grid)
    inf_size = 2

    # education
    a_min = edu_h_ind == 0.0 ? 0.0 : -20.0
    d_κ = 10.0 # need to be updated 
    d_ι = edu_h_ind == 0.0 ? 0.0 : log(edu_h_scale)

    # transition of child dependence
    n_grid = collect(0:n_max)
    n_size = length(n_grid)
    n_Γ = binomial_matrix_function(n_max, p)

    # life-cycle income
    data_h = [
        1.269726,
        1.843193,
        2.291949,
        2.072828,
        2.182088,
        2.182351,
        2.360361,
        2.472632,
        2.46309,
        2.528297,
        2.535929,
        2.583853,
        2.706488,
        2.678228,
        2.690194,
        2.783265,
        2.731313,
        2.736841,
        2.771061,
        2.606311,
        2.728029,
        2.761517,
        2.699816,
        2.712689,
        2.704541,
        2.795719,
        2.734808,
        2.76406,
        2.772031,
        2.775883,
        2.851987,
        2.831737,
        2.832914,
        2.841952,
        2.809018,
        2.858408,
        2.812527,
        2.823082,
        2.762016,
        2.866724,
        2.894183,
        2.728887,
        2.808006,
        2.709969,
        2.63592,
        2.805997,
        2.752519,
        2.464245,
        2.595775,
        2.558831,
        2.34828,
        2.861841,
        2.58545,
        2.382832,
        1.55195,
        2.453288,
        2.336329,
        2.183347,
        2.987182,
        2.529096,
        3.138722,
        3.772982,
        2.509402
    ]
    h_grid = h_function(data_h, age_min, age_max)
    h_grid[(age_ret-age_min+2):end] .= h_grid[age_ret-age_min+1]
    h_size = length(h_grid)
    h_grid[5:end] .= (h_grid[5:end] .+ d_ι)

    h_scale_lr_exp = 1.05
    if lr_exp == 1 # 20-24
        h_i_l = (lr_exp - 1) * 5 + 3
        h_i_r = (lr_exp - 1) * 5 + 7
        h_grid[h_i_l:h_i_r] .= h_scale_lr_exp .* h_grid[h_i_l:h_i_r]
    elseif lr_exp == 2 # 25-29
        h_i_l = (lr_exp - 1) * 5 + 3
        h_i_r = (lr_exp - 1) * 5 + 7
        h_grid[h_i_l:h_i_r] .= h_scale_lr_exp .* h_grid[h_i_l:h_i_r]
    elseif lr_exp == 3 # 30-34
        h_i_l = (lr_exp - 1) * 5 + 3
        h_i_r = (lr_exp - 1) * 5 + 7
        h_grid[h_i_l:h_i_r] .= h_scale_lr_exp .* h_grid[h_i_l:h_i_r]
    elseif lr_exp == 4 # 35-39
        h_i_l = (lr_exp - 1) * 5 + 3
        h_i_r = (lr_exp - 1) * 5 + 7
        h_grid[h_i_l:h_i_r] .= h_scale_lr_exp .* h_grid[h_i_l:h_i_r]
    elseif lr_exp == 5 # 40-44
        h_i_l = (lr_exp - 1) * 5 + 3
        h_i_r = (lr_exp - 1) * 5 + 7
        h_grid[h_i_l:h_i_r] .= h_scale_lr_exp .* h_grid[h_i_l:h_i_r]
    elseif lr_exp == 6 # 45-49
        h_i_l = (lr_exp - 1) * 5 + 3
        h_i_r = (lr_exp - 1) * 5 + 7
        h_grid[h_i_l:h_i_r] .= h_scale_lr_exp .* h_grid[h_i_l:h_i_r]
    end

    # persistent income shock
    ϵ_grid = tauchen_grid(ϵ_size, ρ, σ_ϵ; m=2.0)
    ϵ_Γ = tauchen_transition_matrix(ϵ_grid, ρ, σ_ϵ * edu_σ_ϵ_scale)
    # ϵ_grid, ϵ_Γ = adda_cooper(ϵ_size, ρ, σ_ϵ * edu_σ_ϵ_scale)
    ϵ_G = stationary_distributions(MarkovChain(ϵ_Γ, ϵ_grid))[1]

    # transitory income shock
    ν_grid = tauchen_grid(ν_size, 0.0, σ_ν; m=2.0)
    ν_Γ = tauchen_transition_matrix(ν_grid, 0.0, σ_ν * edu_σ_ν_scale)
    # ν_grid, ν_Γ = adda_cooper(ν_size, 0.0, σ_ν * edu_σ_ν_scale)
    ν_Γ = ν_Γ[1, :]
    ν_G = ν_Γ

    Γ = zeros(n_size, ν_size, ϵ_size, n_size, ϵ_size)
    for ϵ_i in 1:ϵ_size, n_i in 1:n_size, ϵ_p_i in 1:ϵ_size, ν_p_i in 1:ν_size, n_p_i in 1:n_size
        Γ[n_p_i, ν_p_i, ϵ_p_i, n_i, ϵ_i] = n_Γ[n_i, n_p_i] * ν_Γ[ν_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i]
    end

    w_grid = Array{Float64}(undef, ϵ_size, ν_size, h_size)
    for h_i = 1:h_size, ν_i in 1:ν_size, ϵ_i in 1:ϵ_size
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
    if edu_h_ind == 0.0
        a_grid = ((range(0.0, stop=a_size - 1, length=a_size) / (a_size - 1)) .^ a_degree) * a_max
        a_ind_zero = 1
    else
        a_grid_neg = collect(range(a_min, 0.0, length=a_size_neg))
        a_grid_pos = ((range(0.0, stop=a_size - 1, length=a_size) / (a_size - 1)) .^ a_degree) * a_max
        a_grid = vcat(a_grid_neg, a_grid_pos[2:end])
        a_ind_zero = a_size_neg
        a_size = length(a_grid)
    end

    aR_grid = a_grid * R

    # # child quality inputs
    l_grid = collect(0.0:0.5:1.0)
    l_size = length(l_grid)
    x_grid = edu_h_ind == 0.0 ? a_grid : a_grid_pos
    x_size = length(x_grid)

    # indices
    ind_max_ret = collect(Iterators.product(1:ν_size, 1:ϵ_size, 1:a_size))
    ind_ret_inf = collect(Iterators.product(1:ν_size, 1:ϵ_size, 1:n_size, 1:a_size))
    ind_ret_inf_EV = collect(Iterators.product(1:ϵ_size, 1:n_size, 1:a_size))
    ind_inf_min = collect(Iterators.product(1:inf_size, 1:ν_size, 1:ϵ_size, 1:n_size, 1:a_size))
    ind_inf_min_EV = collect(Iterators.product(1:inf_size, 1:ϵ_size, 1:n_size, 1:a_size))

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
        d_κ=d_κ,
        d_ι=d_ι,
        n_max=n_max,
        n_size=n_size,
        n_grid=n_grid,
        n_Γ=n_Γ,
        h_size=h_size,
        h_grid=h_grid,
        data_h=data_h,
        ϵ_size=ϵ_size,
        ϵ_grid=ϵ_grid,
        ϵ_Γ=ϵ_Γ,
        ϵ_G=ϵ_G,
        ν_size=ν_size,
        ν_grid=ν_grid,
        ν_Γ=ν_Γ,
        ν_G=ν_G,
        Γ=Γ,
        w_grid=w_grid,
        P_grid=P_grid,
        q_bar_P_grid=q_bar_P_grid,
        a_min=a_min,
        a_max=a_max,
        a_ind_zero=a_ind_zero,
        a_size=a_size,
        a_size_neg=a_size_neg,
        a_grid=a_grid,
        aR_grid=aR_grid,
        a_degree=a_degree,
        l_size=l_size,
        l_grid=l_grid,
        x_size=x_size,
        x_grid=x_grid,
        q_x=q_x,
        inf_scale=inf_scale,
        edu_h_ind=edu_h_ind,
        edu_h_scale=edu_h_scale,
        edu_σ_ϵ_scale=edu_σ_ϵ_scale,
        edu_σ_ν_scale=edu_σ_ν_scale,
    )
end

mutable struct Mutable_Variables
    """
    Construct a type for mutable variables
    """
    V::Array{Float64,6}
    policy_c::Array{Float64,6}
    policy_a_p::Array{Float64,6}
    policy_e::Array{Float64,6}
    policy_K::Array{Int64,6}
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
    policy_K = ones(Int, a_size, n_size, ϵ_size, ν_size, inf_size, age_size)

    # return outputs
    variables = Mutable_Variables(V, policy_c, policy_a_p, policy_e, policy_K)
    return variables
end

@inline function lininterp(x::AbstractVector{Float64}, y::AbstractVector{Float64}, xq::Float64)
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
        return y0 + w * (y1 - y0)
    end
end

@inline function terminal_step!(
    V_end::AbstractArray{Float64,3},
    policy_c_end::AbstractArray{Float64,3},
    parameters::NamedTuple,
    ν_i::Integer,
    ϵ_i::Integer
)
    V_end_a = @views V_end[:, ϵ_i, ν_i]
    policy_c_end_a = @views policy_c_end[:, ϵ_i, ν_i]

    @unpack w_grid, aR_grid, one_m_γ, inv_one_m_γ = parameters
    @inbounds w_bar = w_grid[ϵ_i, ν_i, end]

    @inbounds for i in eachindex(aR_grid)
        c = aR_grid[i] + w_bar
        policy_c_end_a[i] = c
        V_end_a[i] = u_CRRA(c, one_m_γ, inv_one_m_γ)
    end

    return nothing
end

function retired_step!(
    a_endo::Vector{Float64},
    V_next::AbstractArray{Float64,3},
    policy_c_next::AbstractArray{Float64,3},
    V_current::AbstractArray{Float64,3},
    policy_c_current::AbstractArray{Float64,3},
    policy_a_p_current::AbstractArray{Float64,3},
    parameters::NamedTuple,
    ν_i::Integer,
    ϵ_i::Integer,
    age_i::Integer
)
    V_next_a = @view V_next[:, ϵ_i, ν_i]
    policy_c_next_a = @views policy_c_next[:, ϵ_i, ν_i]

    V_current_a = @view V_current[:, ϵ_i, ν_i]
    policy_c_current_a = @views policy_c_current[:, ϵ_i, ν_i]
    policy_a_p_current_a = @views policy_a_p_current[:, ϵ_i, ν_i]

    @unpack w_grid, a_size, a_grid, a_min, aR_grid, inv_R, EGM_fac, β, one_m_γ, inv_one_m_γ = parameters
    @inbounds w_bar = w_grid[ϵ_i, ν_i, age_i]

    @inbounds for ap_i in 1:a_size
        ap = a_grid[ap_i]
        c = EGM_fac * policy_c_next_a[ap_i]
        m = c + ap
        a_endo[ap_i] = (m - w_bar) * inv_R
    end

    @inbounds for ap_i in 2:a_size
        if a_endo[ap_i] <= a_endo[ap_i-1]
            a_endo[ap_i] = nextfloat(a_endo[ap_i-1])
        end
    end

    ibind = searchsortedlast(a_grid, a_endo[1])

    if ibind > 0
        @inbounds for a_i in 1:ibind
            aR = aR_grid[a_i]
            ap = a_min
            c = aR + w_bar - ap
            policy_a_p_current_a[a_i] = ap
            policy_c_current_a[a_i] = c
            V_current_a[a_i] = u_CRRA(c, one_m_γ, inv_one_m_γ) + β * V_next_a[1]
        end
    end

    @inbounds for a_i in (ibind+1):a_size
        a = a_grid[a_i]
        aR = aR_grid[a_i]
        ap = lininterp(a_endo, a_grid, a)
        c = aR + w_bar - ap
        policy_a_p_current_a[a_i] = ap
        policy_c_current_a[a_i] = c
        Vp = lininterp(a_grid, V_next_a, ap)
        V_current_a[a_i] = u_CRRA(c, one_m_γ, inv_one_m_γ) + β * Vp
    end

    return nothing
end

@inline function solve_e_bisect(
    m::Float64,
    n::Int64,
    P::Float64,
    ψ::Float64,
    γ::Float64,
    κ::Float64,
    one_m_κ::Float64,
    e_bar::Float64;
    c_floor::Float64=1e-12,
    maxit::Int64=80,
    tol::Float64=1e-12
)::Float64

    e_hi = m - c_floor
    if e_hi < e_bar
        return e_hi
    end

    # Define g(e) = marginal benefit of e - marginal cost in terms of forgone c
    # g(e) = ζ*(n/P)^(1-κ) * e^(-κ) - (m - e)^(-γ)
    A = ψ * (n / P)^one_m_κ
    @inline function g(e::Float64)::Float64
        c = m - e
        return A * e^(-κ) - c^(-γ)
    end

    # Check lower bound (constraint) corner: if g(e_bar) <= 0, optimum at e = e_bar
    gl = g(e_bar)
    if gl <= 0.0
        return e_bar
    end

    # Check upper bound corner: if g(e_hi) >= 0, utility still wants more e -> choose e_hi
    gh = g(e_hi)
    if gh >= 0.0
        return e_hi
    end

    # Now we have g(e_min) > 0 and g(e_hi) < 0 -> unique root inside (e_min, e_hi)
    lo = e_bar
    hi = e_hi
    mid = 0.5 * (lo + hi)
    @inbounds for _ in 1:maxit
        mid = 0.5 * (lo + hi)
        gm = g(mid)
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

function infertile_step!(
    a_endo::Vector{Float64},
    V_next::AbstractArray{Float64,3},
    policy_c_next::AbstractArray{Float64,3},
    V_current_n::AbstractArray{Float64,3},
    policy_c_current_n::AbstractArray{Float64,3},
    policy_a_p_current_n::AbstractArray{Float64,3},
    policy_e_current_n::AbstractArray{Float64,3},
    parameters::NamedTuple,
    n_i::Integer,
    ν_i::Integer,
    ϵ_i::Integer,
    age_i::Integer
)
    V_next_a = @view V_next[:, ϵ_i, ν_i]
    policy_c_next_a = @views policy_c_next[:, ϵ_i, ν_i]

    V_current_n_a = @view V_current_n[:, ϵ_i, ν_i]
    policy_c_current_n_a = @views policy_c_current_n[:, ϵ_i, ν_i]
    policy_a_p_current_n_a = @views policy_a_p_current_n[:, ϵ_i, ν_i]
    policy_e_current_n_a = @views policy_e_current_n[:, ϵ_i, ν_i]

    @unpack w_grid, a_size, a_grid, a_min, aR_grid, inv_R, EGM_fac, inv_κ, γ_by_κ, P_grid, q_bar_P_grid = parameters 
    @unpack β, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ = parameters
    @unpack n_grid, ψ, γ, κ, one_m_κ = parameters

    @inbounds w_bar = w_grid[ϵ_i, ν_i, age_i]
    @inbounds P = P_grid[n_i, ϵ_i, ν_i, age_i]
    @inbounds e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
    @inbounds n = n_grid[n_i]

    @inbounds for ap_i in 1:a_size
        ap = a_grid[ap_i]
        c = EGM_fac * policy_c_next_a[ap_i]
        e_foc = (ψ * (n / P)^one_m_κ)^inv_κ * c^γ_by_κ 
        e = max(e_foc, e_bar)
        m = c + ap + e
        a_endo[ap_i] = (m - w_bar) * inv_R
    end

    @inbounds for ap_i in 2:a_size
        if a_endo[ap_i] <= a_endo[ap_i-1]
            a_endo[ap_i] = nextfloat(a_endo[ap_i-1])
        end
    end

    ibind = searchsortedlast(a_grid, a_endo[1])

    if ibind > 0
        @inbounds for a_i in 1:ibind
            aR = aR_grid[a_i]
            M = aR + w_bar
            ap = a_min
            m = M - ap
            if m <= e_bar
                e = e_bar
                c = m - e
                policy_c_current_n_a[a_i] = c
                policy_a_p_current_n_a[a_i] = ap
                policy_e_current_n_a[a_i] = e
                V_current_n_a[a_i] = u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ)                
            else
                e = solve_e_bisect(m, n, P, ψ, γ, κ, one_m_κ, e_bar)
                c = m - e
                policy_a_p_current_n_a[a_i] = ap
                policy_c_current_n_a[a_i] = c
                policy_e_current_n_a[a_i] = e
                Vp = lininterp(a_grid, V_next_a, ap)
                V_current_n_a[a_i] = u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ) + β * Vp
            end
        end
    end

    @inbounds for a_i in (ibind+1):a_size
        a = a_grid[a_i]
        aR = aR_grid[a_i]
        ap = lininterp(a_endo, a_grid, a)
        M = aR + w_bar
        m = M - ap
        if m <= e_bar
            e = e_bar
            c = m - e
            policy_c_current_n_a[a_i] = c
            policy_a_p_current_n_a[a_i] = ap
            policy_e_current_n_a[a_i] = e
            V_current_n_a[a_i] = u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ)
        else
            e = solve_e_bisect(m, n, P, ψ, γ, κ, one_m_κ, e_bar)
            c = m - e
            policy_a_p_current_n_a[a_i] = ap
            policy_c_current_n_a[a_i] = c
            policy_e_current_n_a[a_i] = e
            Vp = lininterp(a_grid, V_next_a, ap)
            V_current_n_a[a_i] = u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ) + β * Vp
        end
    end

    return nothing
end

function fill_EV_Euc!(
    EV_next::AbstractArray{Float64,3},
    Euc_next::AbstractArray{Float64,3},
    V_next::AbstractArray{Float64,4},
    policy_c_next::AbstractArray{Float64,4},
    parameters::NamedTuple
)
    @unpack a_size, n_size, ϵ_size, Γ, γ = parameters
    @inbounds for a_p_i in 1:a_size, n_i in 1:n_size, ϵ_i in 1:ϵ_size
        @views Γt = Γ[:, :, :, n_i, ϵ_i]
        @views Vt = V_next[a_p_i, :, :, :]
        @views ct = policy_c_next[a_p_i, :, :, :]
        EV  = 0.0
        Euc = 0.0
        @simd for k in eachindex(Γt)
            p = Γt[k]
            EV += p * Vt[k]
            c = ct[k]
            Euc += p * (c^(-γ))
        end
        EV_next[a_p_i, n_i, ϵ_i]  = EV
        Euc_next[a_p_i, n_i, ϵ_i] = Euc
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
    @unpack aR_grid, w_grid, one_m_γ, inv_one_m_γ, EGM_fac = parameters
    @unpack inf_size, inf_grid = parameters
    @unpack h_grid = parameters
    @unpack b, r, R, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

    # container
    # c_a = Array{Float64}(undef, (a_size, a_size))
    # for a_i = 1:a_size, a_p_i = 1:a_size
    #     c_a[a_i, a_p_i] = (1.0 + r) * a_grid[a_i] - a_grid[a_p_i]
    # end
    # EV = Array{Float64}(undef, (a_size, n_size, ϵ_size))
    # EV_inf = Array{Float64}(undef, (a_size, n_size, ϵ_size, inf_size))

    println("Solving the HH problem in the last period (at age $age_max)")
    @views V_end = variables.V[:, 1, :, :, 2, end]
    @views policy_c_end = variables.policy_c[:, 1, :, :, 2, end]
    for ν_i in 1:ν_size, ϵ_i in 1:ϵ_size
        terminal_step!(V_end, policy_c_end, parameters, ν_i, ϵ_i)
    end

    println("Solving the HH problem after retirement until the second last period (from age $age_ret to $(age_max-1))")
    a_endo = Vector{Float64}(undef, a_size)
    for age_i = (age_size-1):(-1):(age_ret-age_min+1)
        @views V_next = variables.V[:, 1, :, :, 2, age_i+1]
        @views policy_c_next = variables.policy_c[:, 1, :, :, 2, age_i+1]
        @views V_current = variables.V[:, 1, :, :, 2, age_i]
        @views policy_c_current = variables.policy_c[:, 1, :, :, 2, age_i]
        @views policy_a_p_current = variables.policy_a_p[:, 1, :, :, 2, age_i]
        for ν_i in 1:ν_size, ϵ_i in 1:ϵ_size
            retired_step!(a_endo, V_next, policy_c_next, V_current, policy_c_current, policy_a_p_current, parameters, ν_i, ϵ_i, age_i)
        end
    end

    println("Solving the HH problem just before retirement (at age $(age_ret-1))")
    age_i = age_ret - age_min
    @views V_next = variables.V[:, 1, :, :, 2, age_i+1]
    @views policy_c_next = variables.policy_c[:, 1, :, :, 2, age_i+1]
    @views V_current = variables.V[:, :, :, :, 2, age_i]
    @views policy_c_current = variables.policy_c[:, :, :, :, 2, age_i]
    @views policy_a_p_current = variables.policy_a_p[:, :, :, :, 2, age_i]
    @views policy_e_current = variables.policy_e[:, :, :, :, 2, age_i]
    for ν_i in 1:ν_size, ϵ_i in 1:ϵ_size
        for n_i in 1:n_size
            V_current_n = @view V_current[:, n_i, :, :]
            policy_c_current_n = @views policy_c_current[:, n_i, :, :]
            policy_a_p_current_n = @views policy_a_p_current[:, n_i, :, :]
            policy_e_current_n = @views policy_e_current[:, n_i, :, :]
            if n_i == 1
                retired_step!(a_endo, V_next, policy_c_next, V_current_n, policy_c_current_n, policy_a_p_current_n, parameters, ν_i, ϵ_i, age_i)
            else
                infertile_step!(a_endo, V_next, policy_c_next, V_current_n, policy_c_current_n, policy_a_p_current_n, policy_e_current_n, parameters, n_i, ν_i, ϵ_i, age_i)
            end
        end
    end

    println("Solving the HH problem from the infertile age to before retirement (from age $age_inf to $(age_ret-1))")
    EV_next  = Array{Float64}(undef, a_size, n_size, ϵ_size)
    Euc_next = Array{Float64}(undef, a_size, n_size, ϵ_size)
    for age_i = (age_ret-age_min-1):(-1):(age_inf-age_min+1)
        @views V_next = variables.V[:, :, :, :, 2, age_i+1]
        @views policy_c_next = variables.policy_c[:, :, :, :, 2, age_i+1]
        @views V_current = variables.V[:, :, :, :, 2, age_i]
        @views policy_c_current = variables.policy_c[:, :, :, :, 2, age_i]
        @views policy_a_p_current = variables.policy_a_p[:, :, :, :, 2, age_i]
        @views policy_e_current = variables.policy_e[:, :, :, :, 2, age_i]
        fill_EV_Euc!(EV_next, Euc_next, V_next, policy_c_next, parameters)
    end

    # loop over all states
    for age_i = age_size:(-1):1 # (age_inf-age_min)

        age = age_grid[age_i]
        h = h_grid[age_i]
        println("Solving the problem of HH at age $age...")

        if age == age_max # terminal condition
            Threads.@threads for (ν_i, ϵ_i, a_i) in ind_max_ret
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                w_bar = exp(h + ϵ + ν) * b
                # a = a_grid[a_i]
                # @inbounds variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = utility_function((1.0 + r) * a + w_bar, 0.0, 0.0, γ, ψ, κ, q_bar)
                @inbounds variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = utility_function(c_a[a_i, 1] + w_bar, 0.0, 0.0, γ, ψ, κ, q_bar)
            end

        elseif age_ret < age < age_max # after retirement
            Threads.@threads for (ν_i, ϵ_i, a_i) in ind_max_ret
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                w_bar = exp(h + ϵ + ν) * b
                # a = a_grid[a_i]
                V_best = -10^16
                best_a_p_i = 1
                for a_p_i in 1:a_size
                    # a_p = a_grid[a_p_i]
                    # c = (1.0 + r) * a + w_bar - a_p
                    @inbounds c = c_a[a_i, a_p_i] + w_bar
                    if c > 0.0
                        temp = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar) + β * variables.V[a_p_i, 1, ϵ_i, ν_i, 2, age_i+1]
                        if temp > V_best
                            V_best = temp
                            best_a_p_i = a_p_i
                        end
                    end
                end
                @inbounds variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = V_best
                @inbounds variables.policy_a_p[a_i, 1, ϵ_i, ν_i, 2, age_i] = best_a_p_i
            end
        elseif age == age_ret # at retirement age
            Threads.@threads for (ν_i, ϵ_i, n_i, a_i) in ind_ret_inf
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                w = exp(h + ϵ + ν)
                n = n_grid[n_i]
                # a = a_grid[a_i]
                if n == 0
                    V_best = -10^16
                    best_a_p_i = 1
                    for a_p_i in 1:a_size
                        # a_p = a_grid[a_p_i]
                        # c = (1.0 + r) * a + w - a_p
                        @inbounds c = c_a[a_i, a_p_i] + w
                        if c > 0.0
                            @inbounds temp = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar) + β * variables.V[a_p_i, 1, ϵ_i, ν_i, 2, age_i+1]
                            if temp > V_best
                                V_best = temp
                                best_a_p_i = a_p_i
                            end
                        end
                    end
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_best
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_a_p_i
                else
                    V_best = -10^16
                    best_a_p_i, best_x_i, best_l_i = 1, 1, 1
                    for a_p_i in 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                        # a_p = a_grid[a_p_i]
                        x = x_grid[x_i]
                        l = l_grid[l_i]
                        # c = (1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x
                        @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                        if c > 0.0
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            if q >= q_bar
                                @inbounds temp = utility_function(c, n, q, γ, ψ, κ, q_bar) + β * variables.V[a_p_i, 1, ϵ_i, ν_i, 2, age_i+1]
                                if temp > V_best
                                    V_best = temp
                                    best_a_p_i = a_p_i
                                    best_x_i = x_i
                                    best_l_i = l_i
                                end
                            end
                        end
                    end
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_best
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_a_p_i
                    @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_x_i
                    @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_l_i
                end
            end
        elseif age_inf < age < age_ret # berween infertile age and retirement age
            @inbounds EV .= 0.0
            Threads.@threads for (ϵ_i, n_i, a_p_i) in ind_ret_inf_EV
                for ν_p_i in 1:ν_size, ϵ_p_i = 1:ϵ_size, n_p_i = 1:n_size
                    @inbounds EV[a_p_i, n_i, ϵ_i] += n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[a_p_i, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                end
            end
            Threads.@threads for (ν_i, ϵ_i, n_i, a_i) in ind_ret_inf
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                n = n_grid[n_i]
                # a = a_grid[a_i]
                w = exp(h + ϵ + ν)
                if n == 0
                    V_best = -10^16
                    best_a_p_i = 1
                    for a_p_i = 1:a_size
                        # a_p = a_grid[a_p_i]
                        # c = (1.0 + r) * a + w - a_p
                        @inbounds c = c_a[a_i, a_p_i] + w
                        if c > 0.0
                            @inbounds temp = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar) + β * EV[a_p_i, 1, ϵ_i]
                            if temp > V_best
                                V_best = temp
                                best_a_p_i = a_p_i
                            end
                        end
                    end
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_best
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_a_p_i
                else
                    V_best = -10^16
                    best_a_p_i, best_x_i, best_l_i = 1, 1, 1
                    for a_p_i in 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                        # a_p = a_grid[a_p_i]
                        x = x_grid[x_i]
                        l = l_grid[l_i]
                        # c = (1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x
                        @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                        if c > 0.0
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            if q >= q_bar
                                @inbounds temp = utility_function(c, n, q, γ, ψ, κ, q_bar) + β * EV[a_p_i, n_i, ϵ_i]
                                if temp > V_best
                                    V_best = temp
                                    best_a_p_i = a_p_i
                                    best_x_i = x_i
                                    best_l_i = l_i
                                end
                            end
                        end
                    end
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_best
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_a_p_i
                    @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_x_i
                    @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_l_i
                end
            end
        elseif age == age_inf # about to be infertile
            @inbounds EV .= 0.0
            Threads.@threads for (ϵ_i, n_i, a_p_i) in ind_ret_inf_EV
                for ν_p_i in 1:ν_size, ϵ_p_i = 1:ϵ_size, n_p_i = 1:n_size
                    @inbounds EV[a_p_i, n_i, ϵ_i] += n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[a_p_i, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                end
            end
            Threads.@threads for (f_i, ν_i, ϵ_i, n_i, a_i) in ind_inf_min
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                n = n_grid[n_i]
                # a = a_grid[a_i]
                w = exp(h + ϵ + ν)
                if n == 0

                    if f_i == 1

                        V_best_0, V_best_1 = -10^16, -10^16
                        best_0_a_p_i, best_1_a_p_i = 1, 1
                        for a_p_i = 1:a_size
                            # a_p = a_grid[a_p_i]
                            # c = (1.0 + r) * a + w - a_p
                            @inbounds c = c_a[a_i, a_p_i] + w
                            if c > 0.0
                                u_c = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar)
                                @inbounds temp_0 = u_c + β * EV[a_p_i, 1, ϵ_i]
                                @inbounds temp_1 = u_c + β * EV[a_p_i, 2, ϵ_i]
                                if temp_0 > V_best_0
                                    V_best_0 = temp_0
                                    best_0_a_p_i = a_p_i
                                end
                                if temp_1 > V_best_1
                                    V_best_1 = temp_1
                                    best_1_a_p_i = a_p_i
                                end
                            end
                        end

                        if V_best_0 >= V_best_1
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_0
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_a_p_i
                        else
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_1
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_a_p_i
                            @inbounds variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        end

                    else

                        V_best = -10^16
                        best_a_p_i = 1
                        for a_p_i = 1:a_size
                            # a_p = a_grid[a_p_i]
                            # c = (1.0 + r) * a + w - a_p
                            @inbounds c = c_a[a_i, a_p_i] + w
                            if c > 0.0
                                @inbounds temp = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar) + β * EV[a_p_i, 1, ϵ_i]
                                if temp > V_best
                                    V_best = temp
                                    best_a_p_i = a_p_i
                                end
                            end
                        end
                        @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best
                        @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_a_p_i

                    end

                elseif n == n_max

                    V_best = -10^16
                    best_a_p_i, best_x_i, best_l_i = 1, 1, 1
                    for a_p_i in 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                        # a_p = a_grid[a_p_i]
                        x = x_grid[x_i]
                        l = l_grid[l_i]
                        # c = (1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x
                        @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                        if c > 0.0
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            if q >= q_bar
                                @inbounds temp = utility_function(c, n, q, γ, ψ, κ, q_bar) + β * EV[a_p_i, n_i, ϵ_i]
                                if temp > V_best
                                    V_best = temp
                                    best_a_p_i = a_p_i
                                    best_x_i = x_i
                                    best_l_i = l_i
                                end
                            end
                        end
                    end
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_a_p_i
                    @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_x_i
                    @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_l_i

                else

                    if f_i == 1

                        V_best_0, V_best_1 = -10^16, -10^16
                        best_0_a_p_i, best_1_a_p_i = 1, 1
                        best_0_x_i, best_1_x_i = 1, 1
                        best_0_l_i, best_1_l_i = 1, 1
                        for a_p_i = 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                            # a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            # c = (1.0 + r) * a + w - a_p
                            @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                            if c > 0.0
                                q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                                if q >= q_bar
                                    u_c = utility_function(c, n, q, γ, ψ, κ, q_bar)
                                    @inbounds temp_0 = u_c + β * EV[a_p_i, n_i, ϵ_i]
                                    @inbounds temp_1 = u_c + β * EV[a_p_i, n_i+1, ϵ_i]
                                    if temp_0 > V_best_0
                                        V_best_0 = temp_0
                                        best_0_a_p_i = a_p_i
                                        best_0_x_i = x_i
                                        best_0_l_i = l_i
                                    end
                                    if temp_1 > V_best_1
                                        V_best_1 = temp_1
                                        best_1_a_p_i = a_p_i
                                        best_1_x_i = x_i
                                        best_1_l_i = l_i
                                    end
                                end
                            end
                        end

                        if V_best_0 >= V_best_1
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_0
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_a_p_i
                            @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_x_i
                            @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_l_i

                        else
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_1
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_a_p_i
                            @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_x_i
                            @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_l_i
                            @inbounds variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        end

                    else

                        V_best = -10^16
                        best_a_p_i, best_x_i, best_l_i = 1, 1, 1
                        for a_p_i in 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                            # a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            # c = (1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x
                            @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                            if c > 0.0
                                q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                                if q >= q_bar
                                    @inbounds temp = utility_function(c, n, q, γ, ψ, κ, q_bar) + β * EV[a_p_i, n_i, ϵ_i]
                                    if temp > V_best
                                        V_best = temp
                                        best_a_p_i = a_p_i
                                        best_x_i = x_i
                                        best_l_i = l_i
                                    end
                                end
                            end
                        end
                        @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best
                        @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_a_p_i
                        @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_x_i
                        @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_l_i

                    end
                end
            end
        else # fertile age
            @inbounds EV_inf .= 0.0
            Threads.@threads for (f_i, ϵ_i, n_i, a_p_i) in ind_inf_min_EV
                for ν_p_i in 1:ν_size, ϵ_p_i = 1:ϵ_size, n_p_i = 1:n_size
                    if f_i == 1
                        @inbounds EV_inf[a_p_i, n_i, ϵ_i, f_i] += (1.0 - inf_grid[age_i+1]) * n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[a_p_i, n_p_i, ϵ_p_i, ν_p_i, 1, age_i+1]
                        @inbounds EV_inf[a_p_i, n_i, ϵ_i, f_i] += inf_grid[age_i+1] * n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[a_p_i, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                    else
                        @inbounds EV_inf[a_p_i, n_i, ϵ_i, f_i] += n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[a_p_i, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                    end
                end
            end
            Threads.@threads for (f_i, ν_i, ϵ_i, n_i, a_i) in ind_inf_min
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                n = n_grid[n_i]
                # a = a_grid[a_i]
                w = exp(h + ϵ + ν)
                if n == 0

                    if f_i == 1

                        V_best_0, V_best_1 = -10^16, -10^16
                        best_0_a_p_i, best_1_a_p_i = 1, 1
                        for a_p_i = 1:a_size
                            # a_p = a_grid[a_p_i]
                            # c = (1.0 + r) * a + w - a_p
                            @inbounds c = c_a[a_i, a_p_i] + w
                            if c > 0.0
                                u_c = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar)
                                @inbounds temp_0 = u_c + β * EV_inf[a_p_i, 1, ϵ_i, 1]
                                @inbounds temp_1 = u_c + β * EV_inf[a_p_i, 2, ϵ_i, 1]
                                if temp_0 > V_best_0
                                    V_best_0 = temp_0
                                    best_0_a_p_i = a_p_i
                                end
                                if temp_1 > V_best_1
                                    V_best_1 = temp_1
                                    best_1_a_p_i = a_p_i
                                end
                            end
                        end

                        if V_best_0 >= V_best_1
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_0
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_a_p_i
                        else
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_1
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_a_p_i
                            @inbounds variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        end

                    else

                        V_best = -10^16
                        best_a_p_i = 1
                        for a_p_i = 1:a_size
                            # a_p = a_grid[a_p_i]
                            # c = (1.0 + r) * a + w - a_p
                            @inbounds c = c_a[a_i, a_p_i] + w
                            if c > 0.0
                                @inbounds temp = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar) + β * EV_inf[a_p_i, 1, ϵ_i, 2]
                                if temp > V_best
                                    V_best = temp
                                    best_a_p_i = a_p_i
                                end
                            end
                        end
                        @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best
                        @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_a_p_i

                    end

                elseif n == n_max

                    V_best = -10^16
                    best_a_p_i, best_x_i, best_l_i = 1, 1, 1
                    for a_p_i in 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                        # a_p = a_grid[a_p_i]
                        x = x_grid[x_i]
                        l = l_grid[l_i]
                        # c = (1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x
                        @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                        if c > 0.0
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            if q >= q_bar
                                @inbounds temp = utility_function(c, n, q, γ, ψ, κ, q_bar) + β * EV_inf[a_p_i, n_i, ϵ_i, f_i]
                                if temp > V_best
                                    V_best = temp
                                    best_a_p_i = a_p_i
                                    best_x_i = x_i
                                    best_l_i = l_i
                                end
                            end
                        end
                    end
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_a_p_i
                    @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_x_i
                    @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_l_i

                else
                    if f_i == 1

                        V_best_0, V_best_1 = -10^16, -10^16
                        best_0_a_p_i, best_1_a_p_i = 1, 1
                        best_0_x_i, best_1_x_i = 1, 1
                        best_0_l_i, best_1_l_i = 1, 1
                        for a_p_i = 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                            # a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            # c = (1.0 + r) * a + w - a_p
                            @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                            if c > 0.0
                                q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                                if q >= q_bar
                                    u_c = utility_function(c, n, q, γ, ψ, κ, q_bar)
                                    @inbounds temp_0 = u_c + β * EV_inf[a_p_i, n_i, ϵ_i, f_i]
                                    @inbounds temp_1 = u_c + β * EV_inf[a_p_i, n_i+1, ϵ_i, f_i]
                                    if temp_0 > V_best_0
                                        V_best_0 = temp_0
                                        best_0_a_p_i = a_p_i
                                        best_0_x_i = x_i
                                        best_0_l_i = l_i
                                    end
                                    if temp_1 > V_best_1
                                        V_best_1 = temp_1
                                        best_1_a_p_i = a_p_i
                                        best_1_x_i = x_i
                                        best_1_l_i = l_i
                                    end
                                end
                            end
                        end

                        if V_best_0 >= V_best_1
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_0
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_a_p_i
                            @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_x_i
                            @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_l_i

                        else
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_1
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_a_p_i
                            @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_x_i
                            @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_l_i
                            @inbounds variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        end

                    else

                        V_best = -10^16
                        best_a_p_i, best_x_i, best_l_i = 1, 1, 1
                        for a_p_i in 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                            # a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            # c = (1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x
                            @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                            if c > 0.0
                                q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                                if q >= q_bar
                                    @inbounds temp = utility_function(c, n, q, γ, ψ, κ, q_bar) + β * EV_inf[a_p_i, n_i, ϵ_i, 2]
                                    if temp > V_best
                                        V_best = temp
                                        best_a_p_i = a_p_i
                                        best_x_i = x_i
                                        best_l_i = l_i
                                    end
                                end
                            end
                        end
                        @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best
                        @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_a_p_i
                        @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_x_i
                        @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_l_i
                    end
                end
            end
        end
    end
end

function solve_value_and_policy_edu_function!(variables::Mutable_Variables, parameters::NamedTuple)
    """
    Compute value and policy functions with education choice
    """

    # unpack parameters
    @unpack age_size, age_grid, age_max, age_min, age_ret, age_inf, age_edu = parameters
    @unpack ν_size, ν_grid, ν_Γ = parameters
    @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
    @unpack n_size, n_grid, n_Γ, n_max = parameters
    @unpack a_size, a_grid, a_ind_zero = parameters
    @unpack l_size, l_grid, x_size, x_grid = parameters
    @unpack inf_size, inf_grid = parameters
    @unpack h_grid = parameters
    @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x, d_κ = parameters

    # index 
    ind_max_ret = collect(Iterators.product(1:ν_size, 1:ϵ_size, 1:a_size))
    ind_ret_inf = collect(Iterators.product(1:ν_size, 1:ϵ_size, 1:n_size, 1:a_size))
    ind_ret_inf_EV = collect(Iterators.product(1:ϵ_size, 1:n_size, 1:a_size))
    ind_inf_min = collect(Iterators.product(1:inf_size, 1:ν_size, 1:ϵ_size, 1:n_size, 1:a_size))
    ind_inf_min_EV = collect(Iterators.product(1:inf_size, 1:ϵ_size, 1:n_size, 1:a_size))
    ind_edu = collect(Iterators.product(1:inf_size, 1:ν_size, 1:ϵ_size, 1:a_size))
    ind_edu_EV = collect(Iterators.product(1:inf_size, 1:ϵ_size, 1:a_size))

    # container
    c_a = Array{Float64}(undef, (a_size, a_size))
    for a_i = 1:a_size, a_p_i = 1:a_size
        c_a[a_i, a_p_i] = (1.0 + r) * a_grid[a_i] - a_grid[a_p_i]
    end
    EV = Array{Float64}(undef, (a_size, n_size, ϵ_size))
    EV_inf = Array{Float64}(undef, (a_size, n_size, ϵ_size, inf_size))
    EV_edu = Array{Float64}(undef, (a_size, ϵ_size, inf_size))

    # loop over all states
    for age_i = age_size:(-1):1 # (age_inf-age_min)
        # for age_i = age_size:(-1):(age_ret-age_min+2)
        age = age_grid[age_i]
        h = h_grid[age_i]
        println("Solving the problem of HH at age $age...")
        if age == age_max # terminal condition
            Threads.@threads for (ν_i, ϵ_i, a_i) in ind_max_ret
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                w_bar = exp(h + ϵ + ν) * b
                # a = a_grid[a_i]
                # @inbounds variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = utility_function((1.0 + r) * a + w_bar, 0.0, 0.0, γ, ψ, κ, q_bar)
                @inbounds variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = utility_function(c_a[a_i, 1] + w_bar, 0.0, 0.0, γ, ψ, κ, q_bar)
            end
        elseif age_ret < age < age_max # after retirement
            Threads.@threads for (ν_i, ϵ_i, a_i) in ind_max_ret
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                w_bar = exp(h + ϵ + ν) * b
                # a = a_grid[a_i]
                V_best = -10^16
                best_a_p_i = 1
                for a_p_i in 1:a_size
                    # a_p = a_grid[a_p_i]
                    # c = (1.0 + r) * a + w_bar - a_p
                    @inbounds c = c_a[a_i, a_p_i] + w_bar
                    if c > 0.0
                        temp = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar) + β * variables.V[a_p_i, 1, ϵ_i, ν_i, 2, age_i+1]
                        if temp > V_best
                            V_best = temp
                            best_a_p_i = a_p_i
                        end
                    end
                end
                @inbounds variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = V_best
                @inbounds variables.policy_a_p[a_i, 1, ϵ_i, ν_i, 2, age_i] = best_a_p_i
            end
        elseif age == age_ret # at retirement age
            Threads.@threads for (ν_i, ϵ_i, n_i, a_i) in ind_ret_inf
                # println("($ν_i, $ϵ_i, $n_i, $a_i)")
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                w = exp(h + ϵ + ν)
                n = n_grid[n_i]
                # a = a_grid[a_i]
                if n == 0
                    V_best = -10^16
                    best_a_p_i = 1
                    for a_p_i in 1:a_size
                        # a_p = a_grid[a_p_i]
                        # c = (1.0 + r) * a + w - a_p
                        @inbounds c = c_a[a_i, a_p_i] + w
                        if c > 0.0
                            @inbounds temp = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar) + β * variables.V[a_p_i, 1, ϵ_i, ν_i, 2, age_i+1]
                            if temp > V_best
                                V_best = temp
                                best_a_p_i = a_p_i
                            end
                        end
                    end
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_best
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_a_p_i
                else
                    V_best = -10^16
                    best_a_p_i, best_x_i, best_l_i = 1, 1, 1
                    for a_p_i in 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                        # println("($a_p_i, $x_i, $l_i)")
                        # a_p = a_grid[a_p_i]
                        x = x_grid[x_i]
                        l = l_grid[l_i]
                        # c = (1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x
                        @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                        if c > 0.0
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            if q >= q_bar
                                @inbounds temp = utility_function(c, n, q, γ, ψ, κ, q_bar) + β * variables.V[a_p_i, 1, ϵ_i, ν_i, 2, age_i+1]
                                if temp > V_best
                                    V_best = temp
                                    best_a_p_i = a_p_i
                                    best_x_i = x_i
                                    best_l_i = l_i
                                end
                            end
                        end
                    end
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_best
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_a_p_i
                    @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_x_i
                    @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_l_i
                end
            end
        elseif age_inf < age < age_ret # berween infertile age and retirement age
            @inbounds EV .= 0.0
            Threads.@threads for (ϵ_i, n_i, a_p_i) in ind_ret_inf_EV
                for ν_p_i in 1:ν_size, ϵ_p_i = 1:ϵ_size, n_p_i = 1:n_size
                    @inbounds EV[a_p_i, n_i, ϵ_i] += n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[a_p_i, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                end
            end
            Threads.@threads for (ν_i, ϵ_i, n_i, a_i) in ind_ret_inf
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                n = n_grid[n_i]
                # a = a_grid[a_i]
                w = exp(h + ϵ + ν)
                if n == 0
                    V_best = -10^16
                    best_a_p_i = 1
                    for a_p_i = 1:a_size
                        # a_p = a_grid[a_p_i]
                        # c = (1.0 + r) * a + w - a_p
                        @inbounds c = c_a[a_i, a_p_i] + w
                        if c > 0.0
                            @inbounds temp = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar) + β * EV[a_p_i, 1, ϵ_i]
                            if temp > V_best
                                V_best = temp
                                best_a_p_i = a_p_i
                            end
                        end
                    end
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_best
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_a_p_i
                else
                    V_best = -10^16
                    best_a_p_i, best_x_i, best_l_i = 1, 1, 1
                    for a_p_i in 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                        # a_p = a_grid[a_p_i]
                        x = x_grid[x_i]
                        l = l_grid[l_i]
                        # c = (1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x
                        @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                        if c > 0.0
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            if q >= q_bar
                                @inbounds temp = utility_function(c, n, q, γ, ψ, κ, q_bar) + β * EV[a_p_i, n_i, ϵ_i]
                                if temp > V_best
                                    V_best = temp
                                    best_a_p_i = a_p_i
                                    best_x_i = x_i
                                    best_l_i = l_i
                                end
                            end
                        end
                    end
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_best
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_a_p_i
                    @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_x_i
                    @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_l_i
                end
            end
        elseif age == age_inf # about to be infertile
            @inbounds EV .= 0.0
            Threads.@threads for (ϵ_i, n_i, a_p_i) in ind_ret_inf_EV
                for ν_p_i in 1:ν_size, ϵ_p_i = 1:ϵ_size, n_p_i = 1:n_size
                    @inbounds EV[a_p_i, n_i, ϵ_i] += n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[a_p_i, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                end
            end
            Threads.@threads for (f_i, ν_i, ϵ_i, n_i, a_i) in ind_inf_min
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                n = n_grid[n_i]
                # a = a_grid[a_i]
                w = exp(h + ϵ + ν)
                if n == 0

                    if f_i == 1

                        V_best_0, V_best_1 = -10^16, -10^16
                        best_0_a_p_i, best_1_a_p_i = 1, 1
                        for a_p_i = 1:a_size
                            # a_p = a_grid[a_p_i]
                            # c = (1.0 + r) * a + w - a_p
                            @inbounds c = c_a[a_i, a_p_i] + w
                            if c > 0.0
                                u_c = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar)
                                @inbounds temp_0 = u_c + β * EV[a_p_i, 1, ϵ_i]
                                @inbounds temp_1 = u_c + β * EV[a_p_i, 2, ϵ_i]
                                if temp_0 > V_best_0
                                    V_best_0 = temp_0
                                    best_0_a_p_i = a_p_i
                                end
                                if temp_1 > V_best_1
                                    V_best_1 = temp_1
                                    best_1_a_p_i = a_p_i
                                end
                            end
                        end

                        if V_best_0 >= V_best_1
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_0
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_a_p_i
                        else
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_1
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_a_p_i
                            @inbounds variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        end

                    else

                        V_best = -10^16
                        best_a_p_i = 1
                        for a_p_i = 1:a_size
                            # a_p = a_grid[a_p_i]
                            # c = (1.0 + r) * a + w - a_p
                            @inbounds c = c_a[a_i, a_p_i] + w
                            if c > 0.0
                                @inbounds temp = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar) + β * EV[a_p_i, 1, ϵ_i]
                                if temp > V_best
                                    V_best = temp
                                    best_a_p_i = a_p_i
                                end
                            end
                        end
                        @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best
                        @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_a_p_i

                    end

                elseif n == n_max

                    V_best = -10^16
                    best_a_p_i, best_x_i, best_l_i = 1, 1, 1
                    for a_p_i in 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                        # a_p = a_grid[a_p_i]
                        x = x_grid[x_i]
                        l = l_grid[l_i]
                        # c = (1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x
                        @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                        if c > 0.0
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            if q >= q_bar
                                @inbounds temp = utility_function(c, n, q, γ, ψ, κ, q_bar) + β * EV[a_p_i, n_i, ϵ_i]
                                if temp > V_best
                                    V_best = temp
                                    best_a_p_i = a_p_i
                                    best_x_i = x_i
                                    best_l_i = l_i
                                end
                            end
                        end
                    end
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_a_p_i
                    @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_x_i
                    @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_l_i

                else

                    if f_i == 1

                        V_best_0, V_best_1 = -10^16, -10^16
                        best_0_a_p_i, best_1_a_p_i = 1, 1
                        best_0_x_i, best_1_x_i = 1, 1
                        best_0_l_i, best_1_l_i = 1, 1
                        for a_p_i = 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                            # a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            # c = (1.0 + r) * a + w - a_p
                            @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                            if c > 0.0
                                q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                                if q >= q_bar
                                    u_c = utility_function(c, n, q, γ, ψ, κ, q_bar)
                                    @inbounds temp_0 = u_c + β * EV[a_p_i, n_i, ϵ_i]
                                    @inbounds temp_1 = u_c + β * EV[a_p_i, n_i+1, ϵ_i]
                                    if temp_0 > V_best_0
                                        V_best_0 = temp_0
                                        best_0_a_p_i = a_p_i
                                        best_0_x_i = x_i
                                        best_0_l_i = l_i
                                    end
                                    if temp_1 > V_best_1
                                        V_best_1 = temp_1
                                        best_1_a_p_i = a_p_i
                                        best_1_x_i = x_i
                                        best_1_l_i = l_i
                                    end
                                end
                            end
                        end

                        if V_best_0 >= V_best_1
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_0
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_a_p_i
                            @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_x_i
                            @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_l_i

                        else
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_1
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_a_p_i
                            @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_x_i
                            @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_l_i
                            @inbounds variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        end

                    else

                        V_best = -10^16
                        best_a_p_i, best_x_i, best_l_i = 1, 1, 1
                        for a_p_i in 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                            # a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            # c = (1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x
                            @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                            if c > 0.0
                                q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                                if q >= q_bar
                                    @inbounds temp = utility_function(c, n, q, γ, ψ, κ, q_bar) + β * EV[a_p_i, n_i, ϵ_i]
                                    if temp > V_best
                                        V_best = temp
                                        best_a_p_i = a_p_i
                                        best_x_i = x_i
                                        best_l_i = l_i
                                    end
                                end
                            end
                        end
                        @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best
                        @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_a_p_i
                        @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_x_i
                        @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_l_i

                    end
                end
            end
        elseif age_edu <= age < age_inf # b/w education and fertile age

            @inbounds EV_inf .= 0.0
            Threads.@threads for (f_i, ϵ_i, n_i, a_p_i) in ind_inf_min_EV
                for ν_p_i in 1:ν_size, ϵ_p_i = 1:ϵ_size, n_p_i = 1:n_size
                    if f_i == 1
                        @inbounds EV_inf[a_p_i, n_i, ϵ_i, f_i] += (1.0 - inf_grid[age_i+1]) * n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[a_p_i, n_p_i, ϵ_p_i, ν_p_i, 1, age_i+1]
                        @inbounds EV_inf[a_p_i, n_i, ϵ_i, f_i] += inf_grid[age_i+1] * n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[a_p_i, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                    else
                        @inbounds EV_inf[a_p_i, n_i, ϵ_i, f_i] += n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[a_p_i, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                    end
                end
            end
            Threads.@threads for (f_i, ν_i, ϵ_i, n_i, a_i) in ind_inf_min
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                n = n_grid[n_i]
                # a = a_grid[a_i]
                w = exp(h + ϵ + ν)
                if n == 0

                    if f_i == 1

                        V_best_0, V_best_1 = -10^16, -10^16
                        best_0_a_p_i, best_1_a_p_i = 1, 1
                        for a_p_i = 1:a_size
                            # a_p = a_grid[a_p_i]
                            # c = (1.0 + r) * a + w - a_p
                            @inbounds c = c_a[a_i, a_p_i] + w
                            if c > 0.0
                                u_c = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar)
                                @inbounds temp_0 = u_c + β * EV_inf[a_p_i, 1, ϵ_i, 1]
                                @inbounds temp_1 = u_c + β * EV_inf[a_p_i, 2, ϵ_i, 1]
                                if temp_0 > V_best_0
                                    V_best_0 = temp_0
                                    best_0_a_p_i = a_p_i
                                end
                                if temp_1 > V_best_1
                                    V_best_1 = temp_1
                                    best_1_a_p_i = a_p_i
                                end
                            end
                        end

                        if V_best_0 >= V_best_1
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_0
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_a_p_i
                        else
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_1
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_a_p_i
                            @inbounds variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        end

                    else

                        V_best = -10^16
                        best_a_p_i = 1
                        for a_p_i = 1:a_size
                            # a_p = a_grid[a_p_i]
                            # c = (1.0 + r) * a + w - a_p
                            @inbounds c = c_a[a_i, a_p_i] + w
                            if c > 0.0
                                @inbounds temp = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar) + β * EV_inf[a_p_i, 1, ϵ_i, 2]
                                if temp > V_best
                                    V_best = temp
                                    best_a_p_i = a_p_i
                                end
                            end
                        end
                        @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best
                        @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_a_p_i

                    end

                elseif n == n_max

                    V_best = -10^16
                    best_a_p_i, best_x_i, best_l_i = 1, 1, 1
                    for a_p_i in 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                        # a_p = a_grid[a_p_i]
                        x = x_grid[x_i]
                        l = l_grid[l_i]
                        # c = (1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x
                        @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                        if c > 0.0
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            if q >= q_bar
                                @inbounds temp = utility_function(c, n, q, γ, ψ, κ, q_bar) + β * EV_inf[a_p_i, n_i, ϵ_i, f_i]
                                if temp > V_best
                                    V_best = temp
                                    best_a_p_i = a_p_i
                                    best_x_i = x_i
                                    best_l_i = l_i
                                end
                            end
                        end
                    end
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_a_p_i
                    @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_x_i
                    @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_l_i

                else
                    if f_i == 1

                        V_best_0, V_best_1 = -10^16, -10^16
                        best_0_a_p_i, best_1_a_p_i = 1, 1
                        best_0_x_i, best_1_x_i = 1, 1
                        best_0_l_i, best_1_l_i = 1, 1
                        for a_p_i = 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                            # a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            # c = (1.0 + r) * a + w - a_p
                            @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                            if c > 0.0
                                q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                                if q >= q_bar
                                    u_c = utility_function(c, n, q, γ, ψ, κ, q_bar)
                                    @inbounds temp_0 = u_c + β * EV_inf[a_p_i, n_i, ϵ_i, f_i]
                                    @inbounds temp_1 = u_c + β * EV_inf[a_p_i, n_i+1, ϵ_i, f_i]
                                    if temp_0 > V_best_0
                                        V_best_0 = temp_0
                                        best_0_a_p_i = a_p_i
                                        best_0_x_i = x_i
                                        best_0_l_i = l_i
                                    end
                                    if temp_1 > V_best_1
                                        V_best_1 = temp_1
                                        best_1_a_p_i = a_p_i
                                        best_1_x_i = x_i
                                        best_1_l_i = l_i
                                    end
                                end
                            end
                        end

                        if V_best_0 >= V_best_1
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_0
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_a_p_i
                            @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_x_i
                            @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_0_l_i

                        else
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best_1
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_a_p_i
                            @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_x_i
                            @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_1_l_i
                            @inbounds variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        end

                    else

                        V_best = -10^16
                        best_a_p_i, best_x_i, best_l_i = 1, 1, 1
                        for a_p_i in 1:a_size, x_i in 1:x_size, l_i in 1:l_size
                            # a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            # c = (1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x
                            @inbounds c = c_a[a_i, a_p_i] + (1.0 - l) * w - q_x * x
                            if c > 0.0
                                q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                                if q >= q_bar
                                    @inbounds temp = utility_function(c, n, q, γ, ψ, κ, q_bar) + β * EV_inf[a_p_i, n_i, ϵ_i, 2]
                                    if temp > V_best
                                        V_best = temp
                                        best_a_p_i = a_p_i
                                        best_x_i = x_i
                                        best_l_i = l_i
                                    end
                                end
                            end
                        end
                        @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_best
                        @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_a_p_i
                        @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_x_i
                        @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = best_l_i
                    end
                end
            end
        elseif age == age_min

            @inbounds EV_edu .= 0.0
            ϵ_Γ_4 = ϵ_Γ^4
            no_inf_risk = (1.0 - inf_grid[age_i+1]) * (1.0 - inf_grid[age_i+2]) * (1.0 - inf_grid[age_i+3]) * (1.0 - inf_grid[age_i+4])
            inf_risk = 1.0 - no_inf_risk
            Threads.@threads for (f_i, ϵ_i, a_p_i) in ind_edu_EV
                for ν_p_i in 1:ν_size, ϵ_p_i = 1:ϵ_size
                    if f_i == 1
                        @inbounds EV_edu[a_p_i, ϵ_i, f_i] += no_inf_risk * ϵ_Γ_4[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[a_p_i, 1, ϵ_p_i, ν_p_i, 1, age_i+4]
                        @inbounds EV_edu[a_p_i, ϵ_i, f_i] += inf_risk * ϵ_Γ_4[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[a_p_i, 1, ϵ_p_i, ν_p_i, 2, age_i+4]
                    else
                        @inbounds EV_edu[a_p_i, ϵ_i, f_i] += ϵ_Γ_4[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[a_p_i, 1, ϵ_p_i, ν_p_i, 2, age_i+4]
                    end
                end
            end

            Threads.@threads for (f_i, ν_i, ϵ_i, a_i) in ind_edu
                V_best = -10^16
                best_a_p_i = 1
                for a_p_i = 1:a_size
                    a_p = a_grid[a_p_i]
                    r_kernel = (1.0 + r)^4.0 + (1.0 + r)^3.0 + (1.0 + r)^2.0 + (1.0 + r)
                    c = (-a_p - d_κ) / r_kernel
                    if c > 0.0
                        u_c = utility_function(c, 0.0, 0.0, γ, ψ, κ, q_bar)
                        @inbounds temp = u_c + β * u_c + (β^2.0) * u_c + (β^3.0) * u_c + (β^4.0) * EV_edu[a_p_i, ϵ_i, f_i]
                        if temp > V_best
                            V_best = temp
                            best_a_p_i = a_p_i
                        end
                    end
                end
                @inbounds variables.V[a_i, 1, ϵ_i, ν_i, f_i, age_i:(age_i+3)] .= V_best
                @inbounds variables.policy_a_p[a_i, 1, ϵ_i, ν_i, f_i, age_i:(age_i+3)] .= best_a_p_i
            end
        end
    end
end

function save_JLD_function!(variables::Mutable_Variables, parameters::NamedTuple; filename::String)
    V = variables.V
    policy_a_p = variables.policy_a_p
    policy_x = variables.policy_x
    policy_l = variables.policy_l
    policy_K = variables.policy_K
    @save filename parameters V policy_a_p policy_x policy_l policy_K
    return nothing
end

#==============================#
# solve stationary equilibrium #
#==============================#
parameters = parameters_function()
variables = variables_function(parameters)
solve_value_and_policy_function!(variables, parameters)
save_JLD_function!(variables, parameters, filename="workspace_benchmark.jld2")

parameters_lr_exp_1 = parameters_function(lr_exp=1)
variables_lr_exp_1 = variables_function(parameters_lr_exp_1)
solve_value_and_policy_function!(variables_lr_exp_1, parameters_lr_exp_1)
save_JLD_function!(variables_lr_exp_1, parameters_lr_exp_1, filename="workspace_lr_exp_1.jld2")

parameters_lr_exp_2 = parameters_function(lr_exp=2)
variables_lr_exp_2 = variables_function(parameters_lr_exp_2)
solve_value_and_policy_function!(variables_lr_exp_2, parameters_lr_exp_2)
save_JLD_function!(variables_lr_exp_2, parameters_lr_exp_2, filename="workspace_lr_exp_2.jld2")

parameters_lr_exp_3 = parameters_function(lr_exp=3)
variables_lr_exp_3 = variables_function(parameters_lr_exp_3)
solve_value_and_policy_function!(variables_lr_exp_3, parameters_lr_exp_3)
save_JLD_function!(variables_lr_exp_3, parameters_lr_exp_3, filename="workspace_lr_exp_3.jld2")

parameters_lr_exp_4 = parameters_function(lr_exp=4)
variables_lr_exp_4 = variables_function(parameters_lr_exp_4)
solve_value_and_policy_function!(variables_lr_exp_4, parameters_lr_exp_4)
save_JLD_function!(variables_lr_exp_4, parameters_lr_exp_4, filename="workspace_lr_exp_4.jld2")

parameters_lr_exp_5 = parameters_function(lr_exp=5)
variables_lr_exp_5 = variables_function(parameters_lr_exp_5)
solve_value_and_policy_function!(variables_lr_exp_5, parameters_lr_exp_5)
save_JLD_function!(variables_lr_exp_5, parameters_lr_exp_5, filename="workspace_lr_exp_5.jld2")

parameters_lr_exp_6 = parameters_function(lr_exp=6)
variables_lr_exp_6 = variables_function(parameters_lr_exp_6)
solve_value_and_policy_function!(variables_lr_exp_6, parameters_lr_exp_6)
save_JLD_function!(variables_lr_exp_6, parameters_lr_exp_6, filename="workspace_lr_exp_6.jld2")

# parameters_low_inf = parameters_function(inf_scale=0.5)
# variables_low_inf = variables_function(parameters_low_inf)
# solve_value_and_policy_function!(variables_low_inf, parameters_low_inf)
# save_JLD_function!(variables_low_inf, parameters_low_inf, filename = "workspace_low_inf.jld2")

# parameters_no_inf = parameters_function(inf_scale=0.2)
# variables_no_inf = variables_function(parameters_no_inf)
# solve_value_and_policy_function!(variables_no_inf, parameters_no_inf)
# save_JLD_function!(variables_no_inf, parameters_no_inf, filename = "workspace_no_inf.jld2")

# parameters_edu_h = parameters_function(edu_h_ind=1.0, edu_h_scale=1.25, edu_σ_ϵ_scale=3.00, edu_σ_ν_scale=3.00)
# variables_edu_h = variables_function(parameters_edu_h)
# solve_value_and_policy_edu_function!(variables_edu_h, parameters_edu_h)
# save_JLD_function!(variables_edu_h, parameters_edu_h, filename = "workspace_edu_h.jld2")

# parameters_edu_h_low_σ = parameters_function(edu_h_ind=1.0, edu_h_scale=1.25, edu_σ_ϵ_scale=1.00, edu_σ_ν_scale=1.00)
# variables_edu_h_low_σ = variables_function(parameters_edu_h_low_σ)
# solve_value_and_policy_edu_function!(variables_edu_h_low_σ, parameters_edu_h_low_σ)
# save_JLD_function!(variables_edu_h_low_σ, parameters_edu_h_low_σ, filename = "workspace_edu_h_low_σ.jld2")

#===========#
# simuation #
#===========#
function simulation_function(; num_hh::Int=50000, filename::String)
    """
    simulate variable panels for a given set of policy functions
    """
    # load workspace
    @load filename parameters V policy_a_p policy_x policy_l policy_K

    # set seed
    Random.seed!(1124)

    # simulation periods
    num_periods = parameters.age_size

    # variable panels
    panel_a = ones(Int, num_hh, num_periods)
    panel_a_p = ones(Int, num_hh, num_periods)
    panel_x = ones(num_hh, num_periods)
    panel_l = ones(Int, num_hh, num_periods)
    panel_n = ones(Int, num_hh, num_periods)
    panel_K = ones(Int, num_hh, num_periods)
    shock_ϵ = zeros(Int, num_hh, num_periods)
    shock_ν = zeros(Int, num_hh, num_periods)
    shock_f = zeros(Int, num_hh, num_periods)

    # short-run fertility response
    panel_ΔK = ones(Int, num_hh, num_periods)

    # loop over HHs and Time periods
    for period_i in 1:num_periods
        println("Simulating period = $period_i")
        Threads.@threads for hh_i in 1:num_hh
            if period_i == 1

                @inbounds begin
                    # initiate states
                    panel_a[hh_i, period_i] = 1
                    panel_n[hh_i, period_i] = 1
                    shock_ϵ[hh_i, period_i] = rand(Categorical(vec(parameters.ϵ_G)))
                    shock_ν[hh_i, period_i] = rand(Categorical(vec(parameters.ν_G)))
                    shock_f[hh_i, period_i] = rand(Categorical(vec([1.0 - parameters.inf_grid[period_i], parameters.inf_grid[period_i]])))

                    # actions
                    panel_a_p[hh_i, period_i] = policy_a_p[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                    panel_K[hh_i, period_i] = policy_K[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                    panel_ΔK[hh_i, period_i] = policy_K[panel_a[hh_i, period_i] > parameters.a_size ? panel_a[hh_i, period_i] : panel_a[hh_i, period_i] + 1, panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                end

            elseif 1 < period_i < findall(parameters.age_grid .== (parameters.age_inf + 1))[]

                @inbounds begin
                    # initiate states
                    panel_a[hh_i, period_i] = panel_a_p[hh_i, period_i-1]
                    panel_n[hh_i, period_i] = rand(Categorical(vec(parameters.n_Γ[(panel_n[hh_i, period_i-1]+panel_K[hh_i, period_i-1]-1), :])))
                    shock_ϵ[hh_i, period_i] = rand(Categorical(vec(parameters.ϵ_Γ[shock_ϵ[hh_i, period_i-1], :])))
                    shock_ν[hh_i, period_i] = rand(Categorical(vec(parameters.ν_Γ)))
                    shock_f[hh_i, period_i] = shock_f[hh_i, period_i-1] == 2 ? 2 : rand(Categorical(vec([1.0 - parameters.inf_grid[period_i], parameters.inf_grid[period_i]])))

                    # actions
                    panel_a_p[hh_i, period_i] = policy_a_p[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                    panel_x[hh_i, period_i] = policy_x[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                    panel_l[hh_i, period_i] = policy_l[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                    panel_K[hh_i, period_i] = policy_K[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                    panel_ΔK[hh_i, period_i] = policy_K[panel_a[hh_i, period_i] > parameters.a_size ? panel_a[hh_i, period_i] : panel_a[hh_i, period_i] + 1, panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                end

            elseif findall(parameters.age_grid .== (parameters.age_inf + 1))[] <= period_i <= findall(parameters.age_grid .== parameters.age_ret)[]

                @inbounds begin
                    # initiate states
                    panel_a[hh_i, period_i] = panel_a_p[hh_i, period_i-1]
                    panel_n[hh_i, period_i] = rand(Categorical(vec(parameters.n_Γ[(panel_n[hh_i, period_i-1]+panel_K[hh_i, period_i-1]-1), :])))
                    shock_ϵ[hh_i, period_i] = rand(Categorical(vec(parameters.ϵ_Γ[shock_ϵ[hh_i, period_i-1], :])))
                    shock_ν[hh_i, period_i] = rand(Categorical(vec(parameters.ν_Γ)))
                    shock_f[hh_i, period_i] = shock_f[hh_i, period_i-1] == 2 ? 2 : rand(Categorical(vec([1.0 - parameters.inf_grid[period_i], parameters.inf_grid[period_i]])))

                    # actions
                    panel_a_p[hh_i, period_i] = policy_a_p[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                    panel_x[hh_i, period_i] = policy_x[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                    panel_l[hh_i, period_i] = policy_l[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                    panel_K[hh_i, period_i] = policy_K[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                    panel_ΔK[hh_i, period_i] = policy_K[panel_a[hh_i, period_i] > parameters.a_size ? panel_a[hh_i, period_i] : panel_a[hh_i, period_i] + 1, panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                end

            elseif findall(parameters.age_grid .== (parameters.age_ret))[] < period_i < parameters.age_size

                @inbounds begin
                    # initiate states
                    panel_a[hh_i, period_i] = panel_a_p[hh_i, period_i-1]
                    panel_n[hh_i, period_i] = 1
                    shock_ϵ[hh_i, period_i] = shock_ϵ[hh_i, period_i-1]
                    shock_ν[hh_i, period_i] = shock_ν[hh_i, period_i-1]
                    shock_f[hh_i, period_i] = 2

                    # actions
                    panel_a_p[hh_i, period_i] = policy_a_p[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                end

            else

                @inbounds begin
                    # initiate states
                    panel_a[hh_i, period_i] = panel_a_p[hh_i, period_i-1]
                    panel_n[hh_i, period_i] = 1
                    shock_ϵ[hh_i, period_i] = shock_ϵ[hh_i, period_i-1]
                    shock_ν[hh_i, period_i] = shock_ν[hh_i, period_i-1]
                    shock_f[hh_i, period_i] = 2

                    # actions
                    panel_a_p[hh_i, period_i] = parameters.a_ind_zero
                end
            end
        end
    end
    return panel_a, panel_a_p, panel_x, panel_l, panel_n, panel_K, shock_ϵ, shock_ν, shock_f, panel_ΔK
end

num_hh = 200_000

panel_a, panel_a_p, panel_x, panel_l, panel_n, panel_K, shock_ϵ, shock_ν, shock_f, panel_ΔK = simulation_function(num_hh=num_hh, filename="workspace_benchmark.jld2")

panel_a_lr_exp_1, panel_a_p_lr_exp_1, panel_x_lr_exp_1, panel_l_lr_exp_1, panel_n_lr_exp_1, panel_K_lr_exp_1, shock_ϵ_lr_exp_1, shock_ν_lr_exp_1, shock_f_lr_exp_1, panel_ΔK_lr_exp_1 = simulation_function(num_hh=num_hh, filename="workspace_lr_exp_1.jld2")

panel_a_lr_exp_2, panel_a_p_lr_exp_2, panel_x_lr_exp_2, panel_l_lr_exp_2, panel_n_lr_exp_2, panel_K_lr_exp_2, shock_ϵ_lr_exp_2, shock_ν_lr_exp_2, shock_f_lr_exp_2, panel_ΔK_lr_exp_2 = simulation_function(num_hh=num_hh, filename="workspace_lr_exp_2.jld2")

panel_a_lr_exp_3, panel_a_p_lr_exp_3, panel_x_lr_exp_3, panel_l_lr_exp_3, panel_n_lr_exp_3, panel_K_lr_exp_3, shock_ϵ_lr_exp_3, shock_ν_lr_exp_3, shock_f_lr_exp_3, panel_ΔK_lr_exp_3 = simulation_function(num_hh=num_hh, filename="workspace_lr_exp_3.jld2")

panel_a_lr_exp_4, panel_a_p_lr_exp_4, panel_x_lr_exp_4, panel_l_lr_exp_4, panel_n_lr_exp_4, panel_K_lr_exp_4, shock_ϵ_lr_exp_4, shock_ν_lr_exp_4, shock_f_lr_exp_4, panel_ΔK_lr_exp_4 = simulation_function(num_hh=num_hh, filename="workspace_lr_exp_4.jld2")

panel_a_lr_exp_5, panel_a_p_lr_exp_5, panel_x_lr_exp_5, panel_l_lr_exp_5, panel_n_lr_exp_5, panel_K_lr_exp_5, shock_ϵ_lr_exp_5, shock_ν_lr_exp_5, shock_f_lr_exp_5, panel_ΔK_lr_exp_5 = simulation_function(num_hh=num_hh, filename="workspace_lr_exp_5.jld2")

panel_a_lr_exp_6, panel_a_p_lr_exp_6, panel_x_lr_exp_6, panel_l_lr_exp_6, panel_n_lr_exp_6, panel_K_lr_exp_6, shock_ϵ_lr_exp_6, shock_ν_lr_exp_6, shock_f_lr_exp_6, panel_ΔK_lr_exp_6 = simulation_function(num_hh=num_hh, filename="workspace_lr_exp_6.jld2")

# panel_a_low_inf, panel_a_p_low_inf, panel_x_low_inf, panel_l_low_inf, panel_n_low_inf, panel_K_low_inf, shock_ϵ_low_inf, shock_ν_low_inf, shock_f_low_inf = simulation_function(num_hh = num_hh, filename = "workspace_low_inf.jld2")

# panel_a_no_inf, panel_a_p_no_inf, panel_x_no_inf, panel_l_no_inf, panel_n_no_inf, panel_K_no_inf, shock_ϵ_no_inf, shock_ν_no_inf, shock_f_no_inf = simulation_function(num_hh = num_hh, filename = "workspace_no_inf.jld2")

# panel_a_edu_h, panel_a_p_edu_h, panel_x_edu_h, panel_l_edu_h, panel_n_edu_h, panel_K_edu_h, shock_ϵ_edu_h, shock_ν_edu_h, shock_f_edu_h = simulation_function(num_hh = num_hh, filename = "workspace_edu_h.jld2")

# panel_a_edu_h_low_σ, panel_a_p_edu_h_low_σ, panel_x_edu_h_low_σ, panel_l_edu_h_low_σ, panel_n_edu_h_low_σ, panel_K_edu_h_low_σ, shock_ϵ_edu_h_low_σ, shock_ν_edu_h_low_σ, shock_f_edu_h_low_σ = simulation_function(num_hh = num_hh, filename = "workspace_edu_h_low_σ.jld2")

#====================#
# simulation results #
#====================#
# plot_h_edu_mixed = plot(
#     box=:on,
#     size=[800, 600],
#     xlim=[18, 62],
#     xticks=18:4:62,
#     ylim=[-0.2, 4.2],
#     xtickfont=font(16, "Computer Modern", :black),
#     ytickfont=font(16, "Computer Modern", :black),
#     legendfont=font(16, "Computer Modern", :black),
#     guidefont=font(18, "Computer Modern", :black),
#     titlefont=font(18, "Computer Modern", :black),
#     margin=4mm,
#     xlabel="Age",
#     ylabel="Unit Wage"
# )
# plot_h_edu_mixed = plot!(
#     parameters.age_min:parameters.age_ret,
#     parameters.h_grid[1:(parameters.age_ret-parameters.age_min+1)],
#     label="Benchmark",
#     lw=3,
#     lc=:blue
# )
# plot_h_edu_mixed = plot!(
#     parameters.age_min:parameters.age_ret,
#     vcat(zeros(4), parameters_edu_h.h_grid[5:(parameters.age_ret-parameters.age_min+1)]),
#     label="Education",
#     lw=3,
#     lc=:red,
#     ls=:dashdot
# )
# savefig(plot_h_edu_mixed, string("plot_h_edu_mixed.pdf"))

# # should be average complete fertility rate, not conception right?
# avg_conception_rate = zeros(parameters.age_size)
# avg_conception_rate_low_inf = zeros(parameters.age_size)
# avg_conception_rate_no_inf = zeros(parameters.age_size)
# avg_conception_rate_edu_h = zeros(parameters.age_size)
# avg_conception_rate_edu_h_low_σ = zeros(parameters.age_size)
# Threads.@threads for t = 1:parameters.age_size
#     avg_conception_rate[t] = sum(panel_K[:, t] .- 1) / (num_hh / 1000)
#     avg_conception_rate_low_inf[t] = sum(panel_K_low_inf[:, t] .- 1) / (num_hh / 1000)
#     avg_conception_rate_no_inf[t] = sum(panel_K_no_inf[:, t] .- 1) / (num_hh / 1000)
#     avg_conception_rate_edu_h[t] = sum(panel_K_edu_h[:, t] .- 1) / (num_hh / 1000)
#     avg_conception_rate_edu_h_low_σ[t] = sum(panel_K_edu_h_low_σ[:, t] .- 1) / (num_hh / 1000)
# end
# plot_conception_dist_by_age_mixed = plot(
#     box=:on,
#     size=[800, 600],
#     xlim=[18, 45],
#     xticks=18:3:45,
#     # ylim=[-5, 60],
#     xtickfont=font(16, "Computer Modern", :black),
#     ytickfont=font(16, "Computer Modern", :black),
#     legendfont=font(16, "Computer Modern", :black),
#     guidefont=font(18, "Computer Modern", :black),
#     titlefont=font(18, "Computer Modern", :black),
#     margin=4mm,
#     xlabel="Age",
#     ylabel="Avg Conception Rate"
# )
# plot_conception_dist_by_age_mixed = plot!(
#     parameters.age_min:parameters.age_inf,
#     avg_conception_rate[1:(parameters.age_inf-parameters.age_min+1)],
#     label="Benchmark",
#     lw=3,
#     lc=:blue
# )
# plot_conception_dist_by_age_mixed = plot!(
#     parameters.age_min:parameters.age_inf,
#     avg_conception_rate_low_inf[1:(parameters.age_inf-parameters.age_min+1)],
#     label="Low Infertility Risk",
#     lw=3,
#     lc=:black,
#     ls=:dash
# )
# plot_conception_dist_by_age_mixed = plot!(
#     parameters.age_min:parameters.age_inf,
#     avg_conception_rate_edu_h[1:(parameters.age_inf-parameters.age_min+1)],
#     label="Education (Low Risk)",
#     lw=3,
#     lc=:green,
#     ls=:dot
# )
# plot_conception_dist_by_age_mixed = plot!(
#     parameters.age_min:parameters.age_inf,
#     avg_conception_rate_edu_h_low_σ[1:(parameters.age_inf-parameters.age_min+1)],
#     label="Education",
#     lw=3,
#     lc=:red,
#     ls=:dashdot
# )
# savefig(plot_conception_dist_by_age_mixed, string("plot_conception_dist_by_age_mixed.pdf"))

# # infertility risk
# plot_inf_risk_mixed = plot(
#     box=:on,
#     size=[800, 600],
#     xlim=[18, 54],
#     xticks=18:3:54,
#     ylim=[-0.05, 1.05],
#     xtickfont=font(16, "Computer Modern", :black),
#     ytickfont=font(16, "Computer Modern", :black),
#     legendfont=font(16, "Computer Modern", :black),
#     guidefont=font(18, "Computer Modern", :black),
#     titlefont=font(18, "Computer Modern", :black),
#     margin=4mm,
#     xlabel="Age",
#     ylabel="Probability"
# )
# plot_inf_risk_mixed = scatter!(
#     parameters.data_age,
#     parameters.data_inf,
#     label="Trussell and Wilson (1985)",
#     markersize=7,
#     markercolor=:red,
#     markerstrokewidth=0
# )
# plot_inf_risk_mixed = plot!(
#     parameters.age_min:parameters.age_ret,
#     parameters.inf_grid[1:(parameters.age_ret-parameters.age_min+1)],
#     label="Benchmark",
#     lw=3,
#     lc=:blue
# )
# savefig(plot_inf_risk_mixed, string("plot_inf_risk_data.pdf"))

# plot_inf_risk_mixed = plot(
#     box=:on,
#     size=[800, 600],
#     xlim=[18, 54],
#     xticks=18:3:54,
#     ylim=[-0.05, 1.05],
#     xtickfont=font(16, "Computer Modern", :black),
#     ytickfont=font(16, "Computer Modern", :black),
#     legendfont=font(16, "Computer Modern", :black),
#     guidefont=font(18, "Computer Modern", :black),
#     titlefont=font(18, "Computer Modern", :black),
#     margin=4mm,
#     xlabel="Age",
#     ylabel="Probability")
# plot_inf_risk_mixed = scatter!(
#     parameters.data_age,
#     parameters.data_inf,
#     label="Trussell and Wilson (1985)",
#     markersize=7,
#     markercolor=:red,
#     markerstrokewidth=0
# )
# plot_inf_risk_mixed = plot!(
#     parameters.age_min:parameters.age_ret,
#     parameters.inf_grid[1:(parameters.age_ret-parameters.age_min+1)],
#     label="Benchmark",
#     lw=3,
#     lc=:blue
# )
# plot_inf_risk_mixed = plot!(
#     parameters.age_min:parameters.age_ret,
#     parameters_low_inf.inf_grid[1:(parameters.age_ret-parameters.age_min+1)],
#     label="Low Infertility Risk",
#     lw=3,
#     lc=:black,
#     ls=:dash
# )
# savefig(plot_inf_risk_mixed, string("plot_inf_risk_mixed.pdf"))

# inf_acc = zeros(parameters.age_size)
# inf_acc_low_risk = zeros(parameters_low_inf.age_size)
# inf_acc_no_risk = zeros(parameters_no_inf.age_size)
# inf_grid_temp = 1.0 .- parameters.inf_grid
# inf_grid_low_inf_temp = 1.0 .- parameters_low_inf.inf_grid
# inf_grid_no_inf_temp = 1.0 .- parameters_no_inf.inf_grid
# for age_i = 1:parameters.age_size
#     inf_acc[age_i] = 1.0 - reduce(*, inf_grid_temp[1:age_i])
#     inf_acc_low_risk[age_i] = 1.0 - reduce(*, inf_grid_low_inf_temp[1:age_i])
#     inf_acc_no_risk[age_i] = 1.0 - reduce(*, inf_grid_no_inf_temp[1:age_i])
# end
# plot_inf_risk_mixed = plot(
#     box=:on,
#     size=[800, 600],
#     xlim=[18, 54],
#     xticks=18:3:54,
#     ylim=[-0.05, 1.05],
#     xtickfont=font(16, "Computer Modern", :black),
#     ytickfont=font(16, "Computer Modern", :black),
#     legendfont=font(16, "Computer Modern", :black),
#     guidefont=font(18, "Computer Modern", :black),
#     titlefont=font(18, "Computer Modern", :black),
#     margin=4mm,
#     legend=:bottomright,
#     xlabel="Age",
#     ylabel="Fraction"
# )
# plot_inf_risk_mixed = plot!(
#     parameters.age_min:parameters.age_ret,
#     inf_acc[1:(parameters.age_ret-parameters.age_min+1)],
#     label="Benchmark",
#     lw=3,
#     lc=:blue
# )
# plot_inf_risk_mixed = plot!(
#     parameters.age_min:parameters.age_ret,
#     inf_acc_low_risk[1:(parameters.age_ret-parameters.age_min+1)],
#     label="Low Infertility Risk",
#     lw=3,
#     lc=:black,
#     ls=:dash
# )
# plot_inf_risk_mixed = plot!(
#     parameters.age_min:parameters.age_ret,
#     inf_acc_no_risk[1:(parameters.age_ret-parameters.age_min+1)],
#     label="Low Infertility Risk (Tenth Less)",
#     lw=3,
#     lc=:red,
#     ls=:dashdot
# )
# savefig(plot_inf_risk_mixed, string("plot_inf_risk_mixed_acc.pdf"))

# plot_inf_risk_mixed = plot(
#     box=:on,
#     size=[800, 600],
#     xlim=[18, 54],
#     xticks=18:3:54,
#     ylim=[-0.05, 1.05],
#     xtickfont=font(16, "Computer Modern", :black),
#     ytickfont=font(16, "Computer Modern", :black),
#     legendfont=font(16, "Computer Modern", :black),
#     guidefont=font(18, "Computer Modern", :black),
#     titlefont=font(18, "Computer Modern", :black),
#     margin=4mm,
#     legend=:bottomright,
#     xlabel="Age",
#     ylabel="Fraction"
# )
# plot_inf_risk_mixed = plot!(
#     parameters.age_min:parameters.age_ret,
#     inf_acc[1:(parameters.age_ret-parameters.age_min+1)],
#     label="Benchmark",
#     lw=3,
#     lc=:blue
# )
# savefig(plot_inf_risk_mixed, string("plot_inf_risk_mixed_acc_ben.pdf"))

#====================#
# Short-Run Response #
#====================#
using DataFrames
using CategoricalArrays
using FixedEffectModels

I, J = size(panel_ΔK)
panel_dK = panel_ΔK .- panel_K
panel_a_adj = copy(panel_a)
panel_a_adj[panel_a.==parameters.a_size] .= parameters.a_size - 1
panel_da = parameters.a_grid[panel_a_adj.+1] .- parameters.a_grid[panel_a_adj]
df = DataFrame(
    i=repeat(1:I, J),
    age=repeat(parameters.age_grid, inner=I),
    dK=vec(panel_dK),
    ai=vec(panel_a),
    da=vec(panel_da),
)
df = df[df.ai.!=parameters.a_size, :]
df.age_group = map(df.age) do j
    if 15 <= j <= 19
        "15-19"
    elseif 20 <= j <= 24
        "20-24"
    elseif 25 <= j <= 29
        "25-29"
    elseif 30 <= j <= 34
        "30-34"
    elseif 35 <= j <= 39
        "35-39"
    elseif 40 <= j <= 44
        "40-44"
    elseif 45 <= j <= 50
        "45-50"
    else
        missing
    end
end
dropmissing!(df, [:age_group])
df.age_group = CategoricalArray(df.age_group, ordered=true)

# run the age-group fixed effects model
m = reg(df, @formula(dK ~ 0 + age_group + da & age_group), Vcov.cluster(:i))
# m = reg(df, @formula(dK ~ 0 + da & age_group ), Vcov.cluster(:i))
println(m)
out = DataFrame(term=coefnames(m), b=coef(m), se=stderror(m))
alpha = out[contains.(out.term, "da").&&contains.(out.term, "age_group"), :]

# plot 
x = 1:nrow(alpha)
y = alpha.b
yerr = 1.96 .* alpha.se
age_group = ["15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-50"]
plot_sr = plot(
    x, y,
    # seriestype = :scatter,
    yerror=yerr,
    xlabel="Age group",
    ylabel="α_J",
    xticks=(x, age_group),
    legend=false,
    # xrotation = 45,
)
savefig(plot_sr, string("plot_short_run.pdf"))

#===================#
# Long-Run Response #
#===================#
plot_lr_h = plot(parameters.age_grid, parameters.h_grid, labels="Benchmark", legend=:bottomright)
plot_lr_h = plot!(parameters_lr_exp_1.age_grid, parameters_lr_exp_1.h_grid, labels="20-24")
plot_lr_h = plot!(parameters_lr_exp_2.age_grid, parameters_lr_exp_2.h_grid, labels="25-29")
plot_lr_h = plot!(parameters_lr_exp_3.age_grid, parameters_lr_exp_3.h_grid, labels="30-34")
plot_lr_h = plot!(parameters_lr_exp_4.age_grid, parameters_lr_exp_4.h_grid, labels="35-39")
plot_lr_h = plot!(parameters_lr_exp_5.age_grid, parameters_lr_exp_5.h_grid, labels="40-44")
plot_lr_h = plot!(parameters_lr_exp_6.age_grid, parameters_lr_exp_6.h_grid, labels="45-49")
savefig(plot_lr_h, string("plot_long_run_h.pdf"))

using DataFrames
using CategoricalArrays
using FixedEffectModels

panel_n_max = maximum(parameters.n_grid[panel_n], dims=2)
panel_n_1_max = maximum(parameters.n_grid[panel_n_lr_exp_1], dims=2)
panel_n_2_max = maximum(parameters.n_grid[panel_n_lr_exp_2], dims=2)
panel_n_3_max = maximum(parameters.n_grid[panel_n_lr_exp_3], dims=2)
panel_n_4_max = maximum(parameters.n_grid[panel_n_lr_exp_4], dims=2)
panel_n_5_max = maximum(parameters.n_grid[panel_n_lr_exp_5], dims=2)
panel_n_6_max = maximum(parameters.n_grid[panel_n_lr_exp_6], dims=2)
panel_dn_1_max = panel_n_1_max - panel_n_max
panel_dn_2_max = panel_n_2_max - panel_n_max
panel_dn_3_max = panel_n_3_max - panel_n_max
panel_dn_4_max = panel_n_4_max - panel_n_max
panel_dn_5_max = panel_n_5_max - panel_n_max
panel_dn_6_max = panel_n_6_max - panel_n_max
panel_dN = [panel_dn_1_max; panel_dn_2_max; panel_dn_3_max; panel_dn_4_max; panel_dn_5_max; panel_dn_6_max]

I, J = size(panel_n)
age_group_lr = ["20-24", "25-29", "30-34", "35-39", "40-44", "45-49"]
df_lr = DataFrame(
    i=repeat(1:I, 6),
    exp=repeat(age_group_lr, inner=I),
    dN=vec(panel_dN),
)

# run the long-run model
m_lr = reg(df_lr, @formula(dN ~ 0 + exp), Vcov.cluster(:i))
println(m_lr)
beta = DataFrame(term=coefnames(m_lr), b=coef(m_lr), se=stderror(m_lr))

# plot 
x = 1:nrow(beta)
y = beta.b .* 1000
yerr = 1.96 .* beta.se .* 1000
plot_lr = plot(
    x, y,
    # seriestype = :scatter,
    yerror=yerr,
    xlabel="Age group",
    ylabel="β_J",
    xticks=(x, age_group_lr),
    legend=false,
    # xrotation = 45,
)
savefig(plot_lr, string("plot_long_run.pdf"))