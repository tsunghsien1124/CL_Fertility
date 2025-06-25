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

function utility_function(c::Real, n::Real, q::Real, γ::Real, ψ::Real, κ::Real, q_bar::Real)
    """
    utility function with child quality
        c - consumption
        n - number of kids
        q - child quality
    """
    if c > 0.0
        if n == 0
            return 1.0 / ((1.0 - γ) * c^(γ - 1.0))
        else
            if q >= q_bar
                return 1.0 / ((1.0 - γ) * c^(γ - 1.0)) + ψ / ((1.0 - κ) * (n * q)^(κ - 1.0))
            else
                return -10^16
            end
        end
    else
        return -10^16
    end
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
    q_bar::Real=1.80,               # lower bound on children's consumption # 0.34
    ψ_1::Real=0.91,                 # HH economies to money inpur ro production
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
    ν_size::Integer=3,              # number of transitory shock
    a_max::Real=800,                # max of asset holding
    a_size_neg::Integer=5,          # number of negative asset
    a_size::Integer=50,             # number of asset
    a_degree::Integer=8,            # curvature of asset gridpoints
    q_x::Real=1.0,                  # price of monetary input $x$
    #=================#
    # case indicators #
    #=================#
    edu_ind::Integer=0,             # indicator of edu-dependent life-cycle income
    inf_scale::Real=1.0,            # scale of infertility risk
    σ_ϵ_scale::Real=1.0             # scale of persistent income uncertainty           
)
    """
    Contruct an immutable object containg all paramters
    """
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
    a_min = edu_ind == 0 ? 0.0 : -10.0
    d_κ = 0.1 # need to be updated 
    d_ι = log(1.25)

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
    if edu_ind == 1
        h_grid[5:end] .= (h_grid[5:end] .+ d_ι)
    end

    # persistent income shock
    ϵ_grid = tauchen_grid(ϵ_size, ρ, σ_ϵ)
    ϵ_Γ = tauchen_transition_matrix(ϵ_grid, ρ, σ_ϵ * σ_ϵ_scale)
    ϵ_G = stationary_distributions(MarkovChain(ϵ_Γ, ϵ_grid))[1]

    # transitory income shock
    ν_grid, ν_Γ = adda_cooper(ν_size, 0.0, σ_ν)
    ν_Γ = ν_Γ[1, :]
    ν_G = ν_Γ

    # asset holding
    if edu_ind == 0
        a_grid = ((range(0.0, stop=a_size - 1, length=a_size) / (a_size - 1)) .^ a_degree) * a_max
        a_ind_zero = 1
    else
        a_grid_neg = collect(range(a_min, 0.0, length=a_size_neg))
        a_grid_pos = ((range(0.0, stop=a_size - 1, length=a_size) / (a_size - 1)) .^ a_degree) * a_max
        a_grid = vcat(a_grid_neg, a_grid_pos[2:end])
        a_ind_zero = a_size_neg
        a_size = length(a_grid)
    end

    # child quality inputs
    l_grid = collect(0.0:0.5:1.0)
    l_size = length(l_grid)
    x_grid = edu_ind == 0 ? a_grid : a_grid_pos
    x_size = length(x_grid)

    # return values
    return (
        r=r,
        β=β,
        γ=γ,
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
        a_min=a_min,
        a_max=a_max,
        a_ind_zero=a_ind_zero,
        a_size=a_size,
        a_grid=a_grid,
        a_degree=a_degree,
        l_size=l_size,
        l_grid=l_grid,
        x_size=x_size,
        x_grid=x_grid,
        q_x=q_x,
        edu_ind=edu_ind,
        inf_scale=inf_scale,
        σ_ϵ_scale=σ_ϵ_scale,
    )
end

mutable struct Mutable_Variables
    """
    Construct a type for mutable variables
    """
    V::Array{Float64,6}
    policy_a_p::Array{Int64,6}
    policy_x::Array{Int64,6}
    policy_l::Array{Int64,6}
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
    policy_a_p = ones(Int, a_size, n_size, ϵ_size, ν_size, inf_size, age_size)
    policy_x = ones(Int, a_size, n_size, ϵ_size, ν_size, inf_size, age_size)
    policy_l = ones(Int, a_size, n_size, ϵ_size, ν_size, inf_size, age_size)
    policy_K = ones(Int, a_size, n_size, ϵ_size, ν_size, inf_size, age_size)

    # return outputs
    variables = Mutable_Variables(V, policy_a_p, policy_x, policy_l, policy_K)
    return variables
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
    @unpack a_size, a_grid = parameters
    @unpack l_size, l_grid, x_size, x_grid = parameters
    @unpack inf_size, inf_grid = parameters
    @unpack h_grid = parameters
    @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

    # index 
    ind_max_ret = collect(Iterators.product(1:ν_size, 1:ϵ_size, 1:a_size))
    ind_ret_inf = collect(Iterators.product(1:ν_size, 1:ϵ_size, 1:n_size, 1:a_size))
    ind_ret_inf_EV = collect(Iterators.product(1:ϵ_size, 1:n_size, 1:a_size))
    ind_inf_min = collect(Iterators.product(1:inf_size, 1:ν_size, 1:ϵ_size, 1:n_size, 1:a_size))
    ind_inf_min_EV = collect(Iterators.product(1:inf_size, 1:ϵ_size, 1:n_size, 1:a_size))

    # container
    c_a = Array{Float64}(undef, (a_size, a_size))
    for a_i = 1:a_size, a_p_i = 1:a_size
        c_a[a_i, a_p_i] = (1.0 + r) * a_grid[a_i] - a_grid[a_p_i]
    end
    EV = Array{Float64}(undef, (a_size, n_size, ϵ_size))
    EV_inf = Array{Float64}(undef, (a_size, n_size, ϵ_size, inf_size))

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
                @inbounds variables.V[a_i, 1, ϵ_i, ν_i, f_i, age_i] = V_best
                @inbounds variables.policy_a_p[a_i, 1, ϵ_i, ν_i, f_i, age_i] = best_a_p_i
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
save_JLD_function!(variables, parameters, filename = "workspace_benchmark.jld2")

parameters_edu_μ = parameters_function(edu_ind=1)
variables_edu_μ = variables_function(parameters_edu_μ)
solve_value_and_policy_edu_function!(variables_edu_μ, parameters_edu_μ)
save_JLD_function!(variables_edu_μ, parameters_edu_μ, filename = "workspace_edu_μ.jld2")

parameters_low_inf = parameters_function(inf_scale=0.6)
variables_low_inf = variables_function(parameters_low_inf)
solve_value_and_policy_edu_function!(variables_low_inf, parameters_low_inf)
save_JLD_function!(variables_low_inf, parameters_low_inf, filename = "workspace_low_inf.jld2")

#===========#
# simuation #
#===========#
# load workspace
@load "workspace.jld2" parameters V policy_a_p policy_x policy_l policy_K
@load "workspace_edu.jld2" parameters_edu V_edu policy_a_p_edu policy_x_edu policy_l_edu policy_K_edu
@load "workspace_no_inf_risk.jld2" parameters_no_inf_risk V_no_inf_risk policy_a_p_no_inf_risk policy_x_no_inf_risk policy_l_no_inf_risk policy_K_no_inf_risk

# size of simulation
num_hh = 50000
num_periods = parameters.age_size

# set seed
Random.seed!(1124)

# endogenous state or choice variables
panel_a = ones(Int, num_hh, num_periods)
panel_a_p = ones(Int, num_hh, num_periods)
panel_x = ones(num_hh, num_periods)
panel_l = ones(Int, num_hh, num_periods)
panel_n = ones(Int, num_hh, num_periods)
panel_K = ones(Int, num_hh, num_periods)

panel_a_edu = ones(Int, num_hh, num_periods)
panel_a_p_edu = ones(Int, num_hh, num_periods)
panel_x_edu = ones(num_hh, num_periods)
panel_l_edu = ones(Int, num_hh, num_periods)
panel_n_edu = ones(Int, num_hh, num_periods)
panel_K_edu = ones(Int, num_hh, num_periods)

panel_a_no_inf_risk = ones(Int, num_hh, num_periods)
panel_a_p_no_inf_risk = ones(Int, num_hh, num_periods)
panel_x_no_inf_risk = ones(num_hh, num_periods)
panel_l_no_inf_risk = ones(Int, num_hh, num_periods)
panel_n_no_inf_risk = ones(Int, num_hh, num_periods)
panel_K_no_inf_risk = ones(Int, num_hh, num_periods)

# exogenous variables
shock_n = zeros(Int, num_hh, num_periods)
shock_ϵ = zeros(Int, num_hh, num_periods)
shock_ν = zeros(Int, num_hh, num_periods)
shock_f = zeros(Int, num_hh, num_periods)

shock_n_edu = zeros(Int, num_hh, num_periods)
shock_ϵ_edu = zeros(Int, num_hh, num_periods)
shock_ν_edu = zeros(Int, num_hh, num_periods)
shock_f_edu = zeros(Int, num_hh, num_periods)

shock_n_no_inf_risk = zeros(Int, num_hh, num_periods)
shock_ϵ_no_inf_risk = zeros(Int, num_hh, num_periods)
shock_ν_no_inf_risk = zeros(Int, num_hh, num_periods)
shock_f_no_inf_risk = zeros(Int, num_hh, num_periods)

# Loop over HHs and Time periods
for period_i in 1:parameters.age_size
    println("Simulating period = $period_i")
    Threads.@threads for hh_i in 1:num_hh
        # for hh_i in 1:num_hh
        if period_i == 1

            @inbounds begin
                # initiate states
                panel_a[hh_i, period_i] = 1
                panel_n[hh_i, period_i] = 1
                shock_ϵ[hh_i, period_i] = rand(Categorical(vec(parameters.ϵ_G)))
                shock_ν[hh_i, period_i] = rand(Categorical(vec(parameters.ν_G)))
                shock_f[hh_i, period_i] = rand(Categorical(vec([1.0 - parameters.inf_grid[period_i], parameters.inf_grid[period_i]])))

                panel_a_edu[hh_i, period_i] = parameters_edu.a_ind_zero
                panel_n_edu[hh_i, period_i] = 1
                shock_ϵ_edu[hh_i, period_i] = shock_ϵ[hh_i, period_i]
                shock_ν_edu[hh_i, period_i] = shock_ν[hh_i, period_i]
                shock_f_edu[hh_i, period_i] = shock_f[hh_i, period_i]

                panel_a_no_inf_risk[hh_i, period_i] = 1
                panel_n_no_inf_risk[hh_i, period_i] = 1
                shock_ϵ_no_inf_risk[hh_i, period_i] = shock_ϵ[hh_i, period_i]
                shock_ν_no_inf_risk[hh_i, period_i] = shock_ν[hh_i, period_i]
                shock_f_no_inf_risk[hh_i, period_i] = rand(Categorical(vec([1.0 - parameters_no_inf_risk.inf_grid[period_i], parameters_no_inf_risk.inf_grid[period_i]])))

                # actions
                panel_a_p[hh_i, period_i] = policy_a_p[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                # panel_x[hh_i, period_i] = policy_x[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                # panel_l[hh_i, period_i] = policy_l[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                panel_K[hh_i, period_i] = policy_K[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]

                panel_a_p_edu[hh_i, period_i] = policy_a_p_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]
                # panel_x_edu[hh_i, period_i] = policy_x_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]
                # panel_l_edu[hh_i, period_i] = policy_l_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]
                # panel_K_edu[hh_i, period_i] = policy_K_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]

                panel_a_p_no_inf_risk[hh_i, period_i] = policy_a_p_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
                # panel_x_no_inf_risk[hh_i, period_i] = policy_x_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
                # panel_l_no_inf_risk[hh_i, period_i] = policy_l_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
                # panel_K_no_inf_risk[hh_i, period_i] = policy_K_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
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

                if period_i == 5
                    # initiate states
                    panel_a_edu[hh_i, period_i] = panel_a_p_edu[hh_i, 1]
                    panel_n_edu[hh_i, period_i] = 1
                    shock_ϵ_edu[hh_i, period_i] = shock_ϵ[hh_i, period_i]
                    shock_ν_edu[hh_i, period_i] = shock_ν[hh_i, period_i]
                    shock_f_edu[hh_i, period_i] = shock_f[hh_i, period_i]

                    panel_a_no_inf_risk[hh_i, period_i] = panel_a_p_no_inf_risk[hh_i, period_i-1]
                    panel_n_no_inf_risk[hh_i, period_i] = 1
                    shock_ϵ_no_inf_risk[hh_i, period_i] = shock_ϵ[hh_i, period_i]
                    shock_ν_no_inf_risk[hh_i, period_i] = shock_ν[hh_i, period_i]
                    prob_no_inf_risk = 1.0 - (1.0 - parameters_no_inf_risk.inf_grid[1]) * (1.0 - parameters_no_inf_risk.inf_grid[2]) * (1.0 - parameters_no_inf_risk.inf_grid[3]) * (1.0 - parameters_no_inf_risk.inf_grid[4]) * (1.0 - parameters_no_inf_risk.inf_grid[5])
                    shock_f_no_inf_risk[hh_i, period_i] = rand(Categorical(vec([1.0 - prob_no_inf_risk, prob_no_inf_risk])))

                    # actions
                    panel_a_p_edu[hh_i, period_i] = policy_a_p_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]
                    # panel_x_edu[hh_i, period_i] = policy_x_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]
                    # panel_l_edu[hh_i, period_i] = policy_l_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]
                    panel_K_edu[hh_i, period_i] = policy_K_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]

                    panel_a_p_no_inf_risk[hh_i, period_i] = policy_a_p_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
                    # panel_x_no_inf_risk[hh_i, period_i] = policy_x_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
                    # panel_l_no_inf_risk[hh_i, period_i] = policy_l_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
                    panel_K_no_inf_risk[hh_i, period_i] = policy_K_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
                elseif period_i > 5
                    # initiate states
                    panel_a_edu[hh_i, period_i] = panel_a_p_edu[hh_i, period_i-1]
                    panel_n_edu[hh_i, period_i] = rand(Categorical(vec(parameters_edu.n_Γ[(panel_n_edu[hh_i, period_i-1]+panel_K_edu[hh_i, period_i-1]-1), :])))
                    shock_ϵ_edu[hh_i, period_i] = shock_ϵ[hh_i, period_i]
                    shock_ν_edu[hh_i, period_i] = shock_ν[hh_i, period_i]
                    shock_f_edu[hh_i, period_i] = shock_f[hh_i, period_i]

                    panel_a_no_inf_risk[hh_i, period_i] = panel_a_p_no_inf_risk[hh_i, period_i-1]
                    panel_n_no_inf_risk[hh_i, period_i] = rand(Categorical(vec(parameters_no_inf_risk.n_Γ[(panel_n_no_inf_risk[hh_i, period_i-1]+panel_K_no_inf_risk[hh_i, period_i-1]-1), :])))
                    shock_ϵ_no_inf_risk[hh_i, period_i] = shock_ϵ[hh_i, period_i]
                    shock_ν_no_inf_risk[hh_i, period_i] = shock_ν[hh_i, period_i]
                    shock_f_no_inf_risk[hh_i, period_i] = shock_f_no_inf_risk[hh_i, period_i-1] == 2 ? 2 : rand(Categorical(vec([1.0 - parameters_no_inf_risk.inf_grid[period_i], parameters_no_inf_risk.inf_grid[period_i]])))

                    # actions
                    panel_a_p_edu[hh_i, period_i] = policy_a_p_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]
                    panel_x_edu[hh_i, period_i] = policy_x_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]
                    panel_l_edu[hh_i, period_i] = policy_l_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]
                    panel_K_edu[hh_i, period_i] = policy_K_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]

                    panel_a_p_no_inf_risk[hh_i, period_i] = policy_a_p_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
                    panel_x_no_inf_risk[hh_i, period_i] = policy_x_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
                    panel_l_no_inf_risk[hh_i, period_i] = policy_l_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
                    panel_K_no_inf_risk[hh_i, period_i] = policy_K_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
                end
            end

        elseif findall(parameters.age_grid .== (parameters.age_inf + 1))[] <= period_i <= findall(parameters.age_grid .== parameters.age_ret)[]

            @inbounds begin
                # initiate states
                panel_a[hh_i, period_i] = panel_a_p[hh_i, period_i-1]
                panel_n[hh_i, period_i] = rand(Categorical(vec(parameters.n_Γ[(panel_n[hh_i, period_i-1]+panel_K[hh_i, period_i-1]-1), :])))
                shock_ϵ[hh_i, period_i] = rand(Categorical(vec(parameters.ϵ_Γ[shock_ϵ[hh_i, period_i-1], :])))
                shock_ν[hh_i, period_i] = rand(Categorical(vec(parameters.ν_Γ)))
                shock_f[hh_i, period_i] = shock_f[hh_i, period_i-1] == 2 ? 2 : rand(Categorical(vec([1.0 - parameters.inf_grid[period_i], parameters.inf_grid[period_i]])))

                panel_a_edu[hh_i, period_i] = panel_a_p_edu[hh_i, period_i-1]
                panel_n_edu[hh_i, period_i] = rand(Categorical(vec(parameters_edu.n_Γ[(panel_n_edu[hh_i, period_i-1]+panel_K_edu[hh_i, period_i-1]-1), :])))
                shock_ϵ_edu[hh_i, period_i] = shock_ϵ[hh_i, period_i]
                shock_ν_edu[hh_i, period_i] = shock_ν[hh_i, period_i]
                shock_f_edu[hh_i, period_i] = shock_f[hh_i, period_i]

                panel_a_no_inf_risk[hh_i, period_i] = panel_a_p_no_inf_risk[hh_i, period_i-1]
                panel_n_no_inf_risk[hh_i, period_i] = rand(Categorical(vec(parameters_no_inf_risk.n_Γ[(panel_n_no_inf_risk[hh_i, period_i-1]+panel_K_no_inf_risk[hh_i, period_i-1]-1), :])))
                shock_ϵ_no_inf_risk[hh_i, period_i] = shock_ϵ[hh_i, period_i]
                shock_ν_no_inf_risk[hh_i, period_i] = shock_ν[hh_i, period_i]
                shock_f_no_inf_risk[hh_i, period_i] = shock_f_no_inf_risk[hh_i, period_i-1] == 2 ? 2 : rand(Categorical(vec([1.0 - parameters_no_inf_risk.inf_grid[period_i], parameters_no_inf_risk.inf_grid[period_i]])))

                # actions
                panel_a_p[hh_i, period_i] = policy_a_p[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                panel_x[hh_i, period_i] = policy_x[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                panel_l[hh_i, period_i] = policy_l[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                panel_K[hh_i, period_i] = policy_K[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]

                panel_a_p_edu[hh_i, period_i] = policy_a_p_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]
                panel_x_edu[hh_i, period_i] = policy_x_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]
                panel_l_edu[hh_i, period_i] = policy_l_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]
                panel_K_edu[hh_i, period_i] = policy_K_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]

                panel_a_p_no_inf_risk[hh_i, period_i] = policy_a_p_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
                panel_x_no_inf_risk[hh_i, period_i] = policy_x_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
                panel_l_no_inf_risk[hh_i, period_i] = policy_l_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
                panel_K_no_inf_risk[hh_i, period_i] = policy_K_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
            end

        elseif findall(parameters.age_grid .== (parameters.age_ret))[] < period_i < parameters.age_size

            @inbounds begin
                # initiate states
                panel_a[hh_i, period_i] = panel_a_p[hh_i, period_i-1]
                panel_n[hh_i, period_i] = 1
                shock_ϵ[hh_i, period_i] = shock_ϵ[hh_i, period_i-1]
                shock_ν[hh_i, period_i] = shock_ν[hh_i, period_i-1]
                shock_f[hh_i, period_i] = 2

                panel_a_edu[hh_i, period_i] = panel_a_p_edu[hh_i, period_i-1]
                panel_n_edu[hh_i, period_i] = 1
                shock_ϵ_edu[hh_i, period_i] = shock_ϵ[hh_i, period_i]
                shock_ν_edu[hh_i, period_i] = shock_ν[hh_i, period_i]
                shock_f_edu[hh_i, period_i] = shock_f[hh_i, period_i]

                panel_a_no_inf_risk[hh_i, period_i] = panel_a_p_no_inf_risk[hh_i, period_i-1]
                panel_n_no_inf_risk[hh_i, period_i] = 1
                shock_ϵ_no_inf_risk[hh_i, period_i] = shock_ϵ[hh_i, period_i]
                shock_ν_no_inf_risk[hh_i, period_i] = shock_ν[hh_i, period_i]
                shock_f_no_inf_risk[hh_i, period_i] = 2

                # actions
                panel_a_p[hh_i, period_i] = policy_a_p[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
                panel_a_p_edu[hh_i, period_i] = policy_a_p_edu[panel_a_edu[hh_i, period_i], panel_n_edu[hh_i, period_i], shock_ϵ_edu[hh_i, period_i], shock_ν_edu[hh_i, period_i], shock_f_edu[hh_i, period_i], period_i]
                panel_a_p_no_inf_risk[hh_i, period_i] = policy_a_p_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
            end

        else

            @inbounds begin
                # initiate states
                panel_a[hh_i, period_i] = panel_a_p[hh_i, period_i-1]
                panel_n[hh_i, period_i] = 1
                shock_ϵ[hh_i, period_i] = shock_ϵ[hh_i, period_i-1]
                shock_ν[hh_i, period_i] = shock_ν[hh_i, period_i-1]
                shock_f[hh_i, period_i] = 2

                panel_a_edu[hh_i, period_i] = panel_a_p_edu[hh_i, period_i-1]
                panel_n_edu[hh_i, period_i] = 1
                shock_ϵ_edu[hh_i, period_i] = shock_ϵ[hh_i, period_i]
                shock_ν_edu[hh_i, period_i] = shock_ν[hh_i, period_i]
                shock_f_edu[hh_i, period_i] = shock_f[hh_i, period_i]

                panel_a_no_inf_risk[hh_i, period_i] = panel_a_p_no_inf_risk[hh_i, period_i-1]
                panel_n_no_inf_risk[hh_i, period_i] = 1
                shock_ϵ_no_inf_risk[hh_i, period_i] = shock_ϵ[hh_i, period_i]
                shock_ν_no_inf_risk[hh_i, period_i] = shock_ν[hh_i, period_i]
                shock_f_no_inf_risk[hh_i, period_i] = 2

                # actions
                panel_a_p[hh_i, period_i] = 1
                panel_a_p_edu[hh_i, period_i] = parameters_edu.a_ind_zero
                panel_a_p_no_inf_risk[hh_i, period_i] = parameters_edu.a_ind_zero
            end
        end
    end
end

#====================#
# simulation results #
#====================#
plot_h_edu_mixed = plot(
    box=:on,
    size=[800, 600],
    xlim=[18, 62],
    xticks=18:4:62,
    ylim=[-0.2, 4.2],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm,
    xlabel="Age",
    ylabel="Unit Wage"
)
plot_h_edu_mixed = plot!(
    parameters.age_min:parameters.age_ret,
    parameters.h_grid[1:(parameters.age_ret-parameters.age_min+1)],
    label="Benchmark",
    lw=3,
    lc=:blue
)
plot_h_edu_mixed = plot!(
    parameters.age_min:parameters.age_ret,
    vcat(zeros(4), parameters_edu.h_grid[5:(parameters.age_ret-parameters.age_min+1)]),
    label="Education",
    lw=3,
    lc=:red,
    ls=:dashdot
)
savefig(plot_h_edu_mixed, string("plot_h_edu_mixed.pdf"))

# should be average complete fertility rate, not conception right?
avg_conception_rate = zeros(parameters.age_size)
avg_conception_rate_edu = zeros(parameters.age_size)
avg_conception_rate_no_inf_risk = zeros(parameters.age_size)
Threads.@threads for t = 1:parameters.age_size
    avg_conception_rate[t] = sum(panel_K[:, t] .- 1) / (num_hh / 1000)
    avg_conception_rate_edu[t] = sum(panel_K_edu[:, t] .- 1) / (num_hh / 1000)
    avg_conception_rate_no_inf_risk[t] = sum(panel_K_no_inf_risk[:, t] .- 1) / (num_hh / 1000)
end
plot_conception_dist_by_age_mixed = plot(
    box=:on,
    size=[800, 600],
    xlim=[18, 45],
    xticks=18:3:45,
    ylim=[-5, 60],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm,
    xlabel="Age",
    ylabel="Avg Conception Rate"
)
plot_conception_dist_by_age_mixed = plot!(
    parameters.age_min:parameters.age_inf,
    avg_conception_rate[1:(parameters.age_inf-parameters.age_min+1)],
    label="Benchmark",
    lw=3,
    lc=:blue
)
plot_conception_dist_by_age_mixed = plot!(
    parameters.age_min:parameters.age_inf,
    avg_conception_rate_edu[1:(parameters.age_inf-parameters.age_min+1)],
    label="Education",
    lw=3,
    lc=:red,
    ls=:dashdot
)
plot_conception_dist_by_age_mixed = plot!(
    parameters.age_min:parameters.age_inf,
    avg_conception_rate_no_inf_risk[1:(parameters.age_inf-parameters.age_min+1)],
    label="Education, Low Infertility Risk",
    lw=3,
    lc=:black,
    ls=:dash
)
savefig(plot_conception_dist_by_age_mixed, string("plot_conception_dist_by_age_mixed.pdf"))

# infertility risk
plot_inf_risk_mixed = plot(
    box=:on,
    size=[800, 600],
    xlim=[18, 54],
    xticks=18:3:54,
    ylim=[-0.05, 1.05],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm,
    xlabel="Age",
    ylabel="Probability"
)
plot_inf_risk_mixed = scatter!(
    parameters.data_age,
    parameters.data_inf,
    label="Trussell and Wilson (1985)",
    markersize=7,
    markercolor=:red,
    markerstrokewidth=0
)
plot_inf_risk_mixed = plot!(
    parameters.age_min:parameters.age_ret,
    parameters.inf_grid[1:(parameters.age_ret-parameters.age_min+1)],
    label="Benchmark",
    lw=3,
    lc=:blue
)
savefig(plot_inf_risk_mixed, string("plot_inf_risk_data.pdf"))

plot_inf_risk_mixed = plot(
    box=:on,
    size=[800, 600],
    xlim=[18, 54],
    xticks=18:3:54,
    ylim=[-0.05, 1.05],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm,
    xlabel="Age",
    ylabel="Probability")
plot_inf_risk_mixed = scatter!(
    parameters.data_age,
    parameters.data_inf,
    label="Trussell and Wilson (1985)",
    markersize=7,
    markercolor=:red,
    markerstrokewidth=0
)
plot_inf_risk_mixed = plot!(
    parameters.age_min:parameters.age_ret,
    parameters.inf_grid[1:(parameters.age_ret-parameters.age_min+1)],
    label="Benchmark",
    lw=3,
    lc=:blue
)
plot_inf_risk_mixed = plot!(
    parameters.age_min:parameters.age_ret,
    parameters_no_inf_risk.inf_grid[1:(parameters.age_ret-parameters.age_min+1)],
    label="Low Infertility Risk",
    lw=3,
    lc=:black,
    ls=:dash
)
savefig(plot_inf_risk_mixed, string("plot_inf_risk_mixed.pdf"))

inf_acc = zeros(parameters.age_size)
inf_acc_no_risk = zeros(parameters_no_inf_risk.age_size)
inf_grid_temp = 1.0 .- parameters.inf_grid
inf_grid_no_risk_temp = 1.0 .- parameters_no_inf_risk.inf_grid
for age_i = 1:parameters.age_size
    inf_acc[age_i] = 1.0 - reduce(*, inf_grid_temp[1:age_i])
    inf_acc_no_risk[age_i] = 1.0 - reduce(*, inf_grid_no_risk_temp[1:age_i])
end
plot_inf_risk_mixed = plot(
    box=:on,
    size=[800, 600],
    xlim=[18, 54],
    xticks=18:3:54,
    ylim=[-0.05, 1.05],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm,
    legend=:bottomright,
    xlabel="Age",
    ylabel="Fraction"
)
plot_inf_risk_mixed = plot!(
    parameters.age_min:parameters.age_ret,
    inf_acc[1:(parameters.age_ret-parameters.age_min+1)],
    label="Benchmark",
    lw=3,
    lc=:blue
)
plot_inf_risk_mixed = plot!(
    parameters.age_min:parameters.age_ret,
    inf_acc_no_risk[1:(parameters.age_ret-parameters.age_min+1)],
    label="Low Infertility Risk",
    lw=3,
    lc=:black,
    ls=:dash
)
savefig(plot_inf_risk_mixed, string("plot_inf_risk_mixed_acc.pdf"))

plot_inf_risk_mixed = plot(
    box=:on,
    size=[800, 600],
    xlim=[18, 54],
    xticks=18:3:54,
    ylim=[-0.05, 1.05],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm,
    legend=:bottomright,
    xlabel="Age",
    ylabel="Fraction"
)
plot_inf_risk_mixed = plot!(
    parameters.age_min:parameters.age_ret,
    inf_acc[1:(parameters.age_ret-parameters.age_min+1)],
    label="Benchmark",
    lw=3,
    lc=:blue
)
savefig(plot_inf_risk_mixed, string("plot_inf_risk_mixed_acc_ben.pdf"))


