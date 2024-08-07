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

function adda_cooper(N::Integer, ρ::Real, σ::Real; μ::Real=0.0)
    """
    Approximation of an autoregression process with a Markov chain proposed by Adda and Cooper (2003)
    """
    σ_ϵ = σ / sqrt(1.0 - ρ^2.0)
    ϵ = σ_ϵ .* quantile.(Normal(), [i / N for i = 0:N]) .+ μ
    z = zeros(N)
    for i = 1:N
        if i != (N + 1) / 2
            z[i] = N * σ_ϵ * (pdf(Normal(), (ϵ[i] - μ) / σ_ϵ) - pdf(Normal(), (ϵ[i+1] - μ) / σ_ϵ)) + μ
        end
    end
    Π = zeros(N, N)
    if ρ == 0.0
        Π .= 1.0 / N
    else
        for i = 1:N, j = 1:N
            f(u) = exp(-(u - μ)^2.0 / (2.0 * σ_ϵ^2.0)) * (cdf(Normal(), (ϵ[j+1] - μ * (1.0 - ρ) - ρ * u) / σ) - cdf(Normal(), (ϵ[j] - μ * (1.0 - ρ) - ρ * u) / σ))
            integral, err = quadgk(u -> f(u), ϵ[i], ϵ[i+1])
            Π[i, j] = (N / sqrt(2.0 * π * σ_ϵ^2.0)) * integral
        end
    end
    return z, Π
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
    age_inf::Integer=45,            # infertile age
    age_ret::Integer=65,            # retirement age
    n_max::Integer=10,              # max number of kids
    ϵ_size::Integer=7,              # number of persistent shock
    ν_size::Integer=2,              # number of transitory shock
    a_max::Real=800,                # max of asset holding
    a_size::Integer=50,             # number of asset
    a_degree::Integer=2,            # curvature of asset gridpoints
    q_x::Real=1.0,                  # price of monetary input $x$
    h_edu::Integer=0                # edu-dependent life-cycle income           
)
    """
    Contruct an immutable object containg all paramters
    """
    # infertility parameters: taken from Trussell and Wilson (1985, Population Studies)
    data_inf = [0.07, 0.131, 0.231, 0.345, 0.576, 0.952]
    data_age = [20, 25, 30, 35, 40, 45]
    age_grid, inf_grid = infertility_risk_function(data_age, data_inf, age_min, age_max, age_inf)
    age_size = length(age_grid)
    if h_edu == 2
        inf_grid[1:findall(age_grid .== age_inf)[]] .= 0.0
    end
    inf_size = 2

    # education
    d_size = 2
    d_min = -0.2
    d_κ = 0.1

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
    h_size = length(h_grid)
    if h_edu == 1
        h_grid[1:4] .= h_grid[1] * 0.6
        # h_grid[5:end] = h_grid[5:end]*1.05
        h_g = 0.80
        for t = 5:h_size
            h_grid[t] = h_grid[t] * h_g
            h_g += 0.015
        end
    end

    # persistent income shock
    ϵ_MC = tauchen(ϵ_size, ρ, σ_ϵ, 0.0, 3)
    ϵ_Γ = ϵ_MC.p
    ϵ_grid = collect(ϵ_MC.state_values)
    ϵ_G = stationary_distributions(MarkovChain(ϵ_Γ, ϵ_grid))[1]

    # transitory income shock
    ν_grid, ν_Γ = adda_cooper(ν_size, 0.0, σ_ν)
    ν_Γ = ν_Γ[1, :]
    ν_G = ν_Γ

    # asset holding
    a_grid = ((range(0.0, stop=a_size - 1, length=a_size) / (a_size - 1)) .^ a_degree) * a_max

    # child quality inputs
    l_grid = collect(0.0:0.5:1.0)
    l_size = length(l_grid)
    x_grid = a_grid
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
        age_inf=age_inf,
        age_ret=age_ret,
        age_size=age_size,
        age_grid=age_grid,
        inf_grid=inf_grid,
        data_age=data_age,
        data_inf=data_inf,
        inf_size=inf_size,
        d_size=d_size,
        d_min=d_min,
        d_κ=d_κ,
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
        a_max=a_max,
        a_size=a_size,
        a_grid=a_grid,
        a_degree=a_degree,
        l_size=l_size,
        l_grid=l_grid,
        x_size=x_size,
        x_grid=x_grid,
        q_x=q_x
    )
end

mutable struct Mutable_Variables
    """
    Construct a type for mutable variables
    """
    V::Array{Float64,6}
    policy_a_p::Array{Float64,6}
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

    # define value and policy functions: (a,n,ϵ,ν,f,d,t)
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
    @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
    @unpack ν_size, ν_grid, ν_Γ = parameters
    @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
    @unpack n_size, n_grid, n_Γ, n_max = parameters
    @unpack a_size, a_grid = parameters
    @unpack l_size, l_grid, x_size, x_grid = parameters
    @unpack inf_size, inf_grid = parameters
    @unpack h_grid = parameters
    @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

    # loop over all states
    for age_i = age_size:(-1):1
        age = age_grid[age_i]
        age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
        age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
        println("Solving the problem of HH at age $age...")
        if age == age_max # terminal condition
            Threads.@threads for (ν_i, ϵ_i, a_i) in collect(Iterators.product(1:ν_size, 1:ϵ_size, 1:a_size))
                # for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, a_i = 1:a_size
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                a = a_grid[a_i]
                w_bar = exp(h + ϵ + ν) * b
                @inbounds variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = utility_function((1.0 + r) * a + w_bar, 0.0, 0.0, γ, ψ, κ, q_bar)
            end
        elseif age_ret < age < age_max # after retirement
            Threads.@threads for (ν_i, ϵ_i, a_i) in collect(Iterators.product(1:ν_size, 1:ϵ_size, 1:a_size))
                # for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, a_i = 1:a_size
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                a = a_grid[a_i]
                w_bar = exp(h + ϵ + ν) * b
                EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
                V_all = utility_function.((1.0 + r) * a + w_bar .- a_grid, 0.0, 0.0, γ, ψ, κ, q_bar) .+ β * EV
                V_max_i = argmax(V_all)
                @inbounds variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
                @inbounds variables.policy_a_p[a_i, 1, ϵ_i, ν_i, 2, age_i] = V_max_i # a_grid[V_max_i]
            end
        elseif age == age_ret # at retirement age
            Threads.@threads for (ν_i, ϵ_i, n_i, a_i) in collect(Iterators.product(1:ν_size, 1:ϵ_size, 1:n_size, 1:a_size))
                # for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                n = n_grid[n_i]
                a = a_grid[a_i]
                w = exp(h + ϵ + ν)
                if n == 0
                    EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
                    V_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV
                    V_max_i = argmax(V_all)
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_max_i # a_grid[V_max_i]
                else
                    EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
                    Sol_all = zeros(a_size * x_size * l_size, 4)
                    Sol_all_i = 0
                    for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
                        Sol_all_i += 1
                        @inbounds Sol_all[Sol_all_i, 1] = a_p_i
                        @inbounds Sol_all[Sol_all_i, 2] = x_i
                        @inbounds Sol_all[Sol_all_i, 3] = l_i
                        a_p = a_grid[a_i]
                        x = x_grid[x_i]
                        l = l_grid[l_i]
                        q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                        @inbounds Sol_all[Sol_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV[a_p_i]
                    end
                    V_max_i = argmax(Sol_all[:, 4])
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Sol_all[V_max_i, 4]
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 1]) # a_grid[Int(Sol_all[V_max_i, 1])]
                    @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 2]) # x_grid[Int(Sol_all[V_max_i, 2])]
                    @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 3]) # l_grid[Int(Sol_all[V_max_i, 3])]
                end
            end
        elseif age_inf < age < age_ret # berween infertile age and retirement age
            Threads.@threads for (ν_i, ϵ_i, n_i, a_i) in collect(Iterators.product(1:ν_size, 1:ϵ_size, 1:n_size, 1:a_size))
                # for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                n = n_grid[n_i]
                a = a_grid[a_i]
                w = exp(h + ϵ + ν)
                if n == 0
                    EV = zeros(a_size)
                    for ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                        EV += ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, 1, ϵ_p_i, ν_p_i, 2, age_i+1]
                    end
                    V_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV
                    V_max_i = argmax(V_all)
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_max_i # a_grid[V_max_i]
                else
                    EV = zeros(a_size)
                    for n_p_i = 1:n_size, ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                        EV += n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                    end
                    Sol_all = zeros(a_size * x_size * l_size, 4)
                    Sol_all_i = 0
                    for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
                        Sol_all_i += 1
                        @inbounds Sol_all[Sol_all_i, 1] = a_p_i
                        @inbounds Sol_all[Sol_all_i, 2] = x_i
                        @inbounds Sol_all[Sol_all_i, 3] = l_i
                        a_p = a_grid[a_i]
                        x = x_grid[x_i]
                        l = l_grid[l_i]
                        q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                        @inbounds Sol_all[Sol_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV[a_p_i]
                    end
                    V_max_i = argmax(Sol_all[:, 4])
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Sol_all[V_max_i, 4]
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 1]) # a_grid[Int(Sol_all[V_max_i, 1])]
                    @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 2]) # x_grid[Int(Sol_all[V_max_i, 2])]
                    @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 3]) # l_grid[Int(Sol_all[V_max_i, 3])]
                end
            end
        elseif age == age_inf # about to be infertile
            Threads.@threads for (f_i, ν_i, ϵ_i, n_i, a_i) in collect(Iterators.product(1:inf_size, 1:ν_size, 1:ϵ_size, 1:n_size, 1:a_size))
                # for f_i = 1:inf_size, ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                n = n_grid[n_i]
                a = a_grid[a_i]
                w = exp(h + ϵ + ν)
                if n == 0
                    if f_i == 1
                        EV_K_0 = zeros(a_size)
                        for ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_0 += ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, 1, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        V_K_0_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV_K_0
                        V_K_0_max_i = argmax(V_K_0_all)
                        EV_K_1 = zeros(a_size)
                        for n_p_i = 1:n_size, ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_1 += n_Γ[n_i+1, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        V_K_1_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV_K_1
                        V_K_1_max_i = argmax(V_K_1_all)
                        if V_K_1_all[V_K_1_max_i] >= V_K_0_all[V_K_0_max_i]
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_1_all[V_K_1_max_i]
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_1_max_i # a_grid[V_K_1_max_i]
                            @inbounds variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        else
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_all[V_K_0_max_i]
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_max_i # a_grid[V_K_0_max_i]
                        end
                    else
                        EV_K_0 = zeros(a_size)
                        for ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_0 += ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, 1, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        V_K_0_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV_K_0
                        V_K_0_max_i = argmax(V_K_0_all)
                        @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_all[V_K_0_max_i]
                        @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_max_i # a_grid[V_K_0_max_i]
                    end
                elseif n == n_max
                    EV_K_0 = zeros(a_size)
                    for n_p_i = 1:n_size, ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                        EV_K_0 += n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                    end
                    Sol_K_0_all = zeros(a_size * x_size * l_size, 4)
                    Sol_K_0_all_i = 0
                    for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
                        Sol_K_0_all_i += 1
                        @inbounds Sol_K_0_all[Sol_K_0_all_i, 1] = a_p_i
                        @inbounds Sol_K_0_all[Sol_K_0_all_i, 2] = x_i
                        @inbounds Sol_K_0_all[Sol_K_0_all_i, 3] = l_i
                        a_p = a_grid[a_i]
                        x = x_grid[x_i]
                        l = l_grid[l_i]
                        q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                        @inbounds Sol_K_0_all[Sol_K_0_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_0[a_p_i]
                    end
                    V_K_0_max_i = argmax(Sol_K_0_all[:, 4])
                    @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_0_all[V_K_0_max_i, 4]
                    @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 1]) # a_grid[Int(Sol_K_0_all[V_K_0_max_i, 1])]
                    @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 2]) # x_grid[Int(Sol_K_0_all[V_K_0_max_i, 2])]
                    @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 3]) # l_grid[Int(Sol_K_0_all[V_K_0_max_i, 3])]
                else
                    if f_i == 1
                        EV_K_0 = zeros(a_size)
                        for n_p_i = 1:n_size, ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_0 += n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        Sol_K_0_all = zeros(a_size * x_size * l_size, 4)
                        Sol_K_0_all_i = 0
                        for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
                            Sol_K_0_all_i += 1
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 1] = a_p_i
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 2] = x_i
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 3] = l_i
                            a_p = a_grid[a_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_0[a_p_i]
                        end
                        V_K_0_max_i = argmax(Sol_K_0_all[:, 4])
                        EV_K_1 = zeros(a_size)
                        for n_p_i = 1:n_size, ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_1 += n_Γ[n_i+1, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        Sol_K_1_all = zeros(a_size * x_size * l_size, 4)
                        Sol_K_1_all_i = 0
                        for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
                            Sol_K_1_all_i += 1
                            @inbounds Sol_K_1_all[Sol_K_1_all_i, 1] = a_p_i
                            @inbounds Sol_K_1_all[Sol_K_1_all_i, 2] = x_i
                            @inbounds Sol_K_1_all[Sol_K_1_all_i, 3] = l_i
                            a_p = a_grid[a_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            @inbounds Sol_K_1_all[Sol_K_1_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_1[a_p_i]
                        end
                        V_K_1_max_i = argmax(Sol_K_1_all[:, 4])
                        if Sol_K_1_all[V_K_1_max_i, 4] >= Sol_K_0_all[V_K_0_max_i, 4]
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_1_all[V_K_1_max_i, 4]
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_1_all[V_K_1_max_i, 1]) # a_grid[Int(Sol_K_1_all[V_K_1_max_i, 1])]
                            @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_1_all[V_K_1_max_i, 2]) # x_grid[Int(Sol_K_1_all[V_K_1_max_i, 2])]
                            @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_1_all[V_K_1_max_i, 3]) # l_grid[Int(Sol_K_1_all[V_K_1_max_i, 3])]
                            @inbounds variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        else
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_0_all[V_K_0_max_i, 4]
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 1]) # a_grid[Int(Sol_K_0_all[V_K_0_max_i, 1])]
                            @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 2]) # x_grid[Int(Sol_K_0_all[V_K_0_max_i, 2])]
                            @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 3]) # l_grid[Int(Sol_K_0_all[V_K_0_max_i, 3])]
                        end
                    else
                        EV_K_0 = zeros(a_size)
                        for n_p_i = 1:n_size, ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_0 += n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        Sol_K_0_all = zeros(a_size * x_size * l_size, 4)
                        Sol_K_0_all_i = 0
                        for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
                            Sol_K_0_all_i += 1
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 1] = a_p_i
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 2] = x_i
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 3] = l_i
                            a_p = a_grid[a_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_0[a_p_i]
                        end
                        V_K_0_max_i = argmax(Sol_K_0_all[:, 4])
                        @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_0_all[V_K_0_max_i, 4]
                        @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 1]) # a_grid[Int(Sol_K_0_all[V_K_0_max_i, 1])]
                        @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 2]) # x_grid[Int(Sol_K_0_all[V_K_0_max_i, 2])]
                        @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 3]) # l_grid[Int(Sol_K_0_all[V_K_0_max_i, 3])]
                    end
                end
            end
        else # fertile age
            Threads.@threads for (f_i, ν_i, ϵ_i, n_i, a_i) in collect(Iterators.product(1:inf_size, 1:ν_size, 1:ϵ_size, 1:n_size, 1:a_size))
                # for f_i = 1:inf_size, ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                n = n_grid[n_i]
                a = a_grid[a_i]
                w = exp(h + ϵ + ν)
                inf_risk = inf_grid[age_i+1]
                if n == 0
                    if f_i == 1
                        EV_K_0 = zeros(a_size)
                        for ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_0 += (1.0 - inf_risk) * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_i, ϵ_p_i, ν_p_i, 1, age_i+1] + inf_risk * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        V_K_0_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV_K_0
                        V_K_0_max_i = argmax(V_K_0_all)
                        EV_K_1 = zeros(a_size)
                        for n_p_i = 1:n_size, ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_1 += (1.0 - inf_risk) * n_Γ[n_i+1, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 1, age_i+1] + inf_risk * n_Γ[n_i+1, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        V_K_1_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV_K_1
                        V_K_1_max_i = argmax(V_K_1_all)
                        if V_K_1_all[V_K_1_max_i] >= V_K_0_all[V_K_0_max_i]
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_1_all[V_K_1_max_i]
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_1_max_i # a_grid[V_K_1_max_i]
                            @inbounds variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        else
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_all[V_K_0_max_i]
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_max_i # a_grid[V_K_0_max_i]
                        end
                    else
                        EV_K_0 = zeros(a_size)
                        for ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_0 += ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        V_K_0_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV_K_0
                        V_K_0_max_i = argmax(V_K_0_all)
                        @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_all[V_K_0_max_i]
                        @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_max_i # a_grid[V_K_0_max_i]
                    end
                elseif n == n_max
                    if f_i == 1
                        EV_K_0 = zeros(a_size)
                        for n_p_i = 1:n_size, ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_0 += (1.0 - inf_risk) * n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 1, age_i+1] + inf_risk * n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        Sol_K_0_all = zeros(a_size * x_size * l_size, 4)
                        Sol_K_0_all_i = 0
                        for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
                            Sol_K_0_all_i += 1
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 1] = a_p_i
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 2] = x_i
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 3] = l_i
                            a_p = a_grid[a_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_0[a_p_i]
                        end
                        V_K_0_max_i = argmax(Sol_K_0_all[:, 4])
                        @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_0_all[V_K_0_max_i, 4]
                        @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 1]) # a_grid[Int(Sol_K_0_all[V_K_0_max_i, 1])]
                        @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 2]) # x_grid[Int(Sol_K_0_all[V_K_0_max_i, 2])]
                        @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 3]) # l_grid[Int(Sol_K_0_all[V_K_0_max_i, 3])]
                    else
                        EV_K_0 = zeros(a_size)
                        for n_p_i = 1:n_size, ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_0 += n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        Sol_K_0_all = zeros(a_size * x_size * l_size, 4)
                        Sol_K_0_all_i = 0
                        for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
                            Sol_K_0_all_i += 1
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 1] = a_p_i
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 2] = x_i
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 3] = l_i
                            a_p = a_grid[a_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_0[a_p_i]
                        end
                        V_K_0_max_i = argmax(Sol_K_0_all[:, 4])
                        @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_0_all[V_K_0_max_i, 4]
                        @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 1]) # a_grid[Int(Sol_K_0_all[V_K_0_max_i, 1])]
                        @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 2]) # x_grid[Int(Sol_K_0_all[V_K_0_max_i, 2])]
                        @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 3]) # l_grid[Int(Sol_K_0_all[V_K_0_max_i, 3])]
                    end
                else
                    if f_i == 1
                        EV_K_0 = zeros(a_size)
                        for n_p_i = 1:n_size, ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_0 += (1.0 - inf_risk) * n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 1, age_i+1] + inf_risk * n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        Sol_K_0_all = zeros(a_size * x_size * l_size, 4)
                        Sol_K_0_all_i = 0
                        for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
                            Sol_K_0_all_i += 1
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 1] = a_p_i
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 2] = x_i
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 3] = l_i
                            a_p = a_grid[a_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_0[a_p_i]
                        end
                        V_K_0_max_i = argmax(Sol_K_0_all[:, 4])
                        EV_K_1 = zeros(a_size)
                        for n_p_i = 1:n_size, ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_1 += (1.0 - inf_risk) * n_Γ[n_i+1, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 1, age_i+1] + inf_risk * n_Γ[n_i+1, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        Sol_K_1_all = zeros(a_size * x_size * l_size, 4)
                        Sol_K_1_all_i = 0
                        for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
                            Sol_K_1_all_i += 1
                            @inbounds Sol_K_1_all[Sol_K_1_all_i, 1] = a_p_i
                            @inbounds Sol_K_1_all[Sol_K_1_all_i, 2] = x_i
                            @inbounds Sol_K_1_all[Sol_K_1_all_i, 3] = l_i
                            a_p = a_grid[a_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            @inbounds Sol_K_1_all[Sol_K_1_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_1[a_p_i]
                        end
                        V_K_1_max_i = argmax(Sol_K_1_all[:, 4])
                        if Sol_K_1_all[V_K_1_max_i, 4] >= Sol_K_0_all[V_K_0_max_i, 4]
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_1_all[V_K_1_max_i, 4]
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_1_all[V_K_1_max_i, 1]) # a_grid[Int(Sol_K_1_all[V_K_1_max_i, 1])]
                            @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_1_all[V_K_1_max_i, 2]) # x_grid[Int(Sol_K_1_all[V_K_1_max_i, 2])]
                            @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_1_all[V_K_1_max_i, 3]) # l_grid[Int(Sol_K_1_all[V_K_1_max_i, 3])]
                            @inbounds variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        else
                            @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_0_all[V_K_0_max_i, 4]
                            @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 1]) # a_grid[Int(Sol_K_0_all[V_K_0_max_i, 1])]
                            @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 2]) # x_grid[Int(Sol_K_0_all[V_K_0_max_i, 2])]
                            @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 3])  # l_grid[Int(Sol_K_0_all[V_K_0_max_i, 3])]
                        end
                    else
                        EV_K_0 = zeros(a_size)
                        for n_p_i = 1:n_size, ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_0 += n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        Sol_K_0_all = zeros(a_size * x_size * l_size, 4)
                        Sol_K_0_all_i = 0
                        for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
                            Sol_K_0_all_i += 1
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 1] = a_p_i
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 2] = x_i
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 3] = l_i
                            a_p = a_grid[a_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            @inbounds Sol_K_0_all[Sol_K_0_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_0[a_p_i]
                        end
                        V_K_0_max_i = argmax(Sol_K_0_all[:, 4])
                        @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_0_all[V_K_0_max_i, 4]
                        @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 1]) # a_grid[Int(Sol_K_0_all[V_K_0_max_i, 1])]
                        @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 2]) # x_grid[Int(Sol_K_0_all[V_K_0_max_i, 2])]
                        @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 3])  # l_grid[Int(Sol_K_0_all[V_K_0_max_i, 3])]
                    end
                end
            end
        end
    end
end

#==============================#
# solve stationary equilibrium #
#==============================#
# function solve_all()
#     parameters = parameters_function()
#     variables = variables_function(parameters)
#     solve_value_and_policy_function!(variables, parameters)
#     return parameters, variables
# end
# @btime parameters, variables = solve_all()

# V = variables.V
# policy_a_p = variables.policy_a_p
# policy_x = variables.policy_x
# policy_l = variables.policy_l
# policy_K = variables.policy_K
# @save "workspace.jld2" parameters V policy_a_p policy_x policy_l policy_K

# parameters_h_edu = parameters_function(h_edu=1)
# variables_h_edu = variables_function(parameters_h_edu)
# solve_value_and_policy_function!(variables_h_edu, parameters_h_edu)
# V_h_edu = variables_h_edu.V
# policy_a_p_h_edu = variables_h_edu.policy_a_p
# policy_x_h_edu = variables_h_edu.policy_x
# policy_l_h_edu = variables_h_edu.policy_l
# policy_K_h_edu = variables_h_edu.policy_K
# @save "workspace_h_edu.jld2" parameters_h_edu V_h_edu policy_a_p_h_edu policy_x_h_edu policy_l_h_edu policy_K_h_edu

# parameters_no_inf_risk = parameters_function(h_edu=2)
# variables_no_inf_risk = variables_function(parameters_no_inf_risk)
# solve_value_and_policy_function!(variables_no_inf_risk, parameters_no_inf_risk)
# V_no_inf_risk = variables_no_inf_risk.V
# policy_a_p_no_inf_risk = variables_no_inf_risk.policy_a_p
# policy_x_no_inf_risk = variables_no_inf_risk.policy_x
# policy_l_no_inf_risk = variables_no_inf_risk.policy_l
# policy_K_no_inf_risk = variables_no_inf_risk.policy_K
# @save "workspace_no_inf_risk.jld2" parameters_no_inf_risk V_no_inf_risk policy_a_p_no_inf_risk policy_x_no_inf_risk policy_l_no_inf_risk policy_K_no_inf_risk
