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
    model(t, ω) = ω[1] * exp.(ω[2] * t)
    ω_int = [0.5, 0.5]
    fit = curve_fit(model, data_age, data_inf, ω_int)
    model_age = collect(age_min:age_max)
    model_inf = fit.param[1] .* exp.(fit.param[2] .* model_age)
    model_inf[findall(model_age .== (age_inf + 1))[]:end] .= 1.0
    return model_age, model_inf
end

function h_function(data_h::Array{Float64,1}, age_min::Integer, age_ret::Integer)
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
    contruct an immutable object containg all paramters
    """
    # infertility parameters: taken from Trussell and Wilson (1985, Population Studies)
    data_inf = [0.07, 0.131, 0.231, 0.345, 0.576, 0.952]
    data_age = [20, 25, 30, 35, 40, 45]
    age_grid, inf_grid = infertility_risk_function(data_age, data_inf, age_min, age_max, age_inf)
    age_size = length(age_grid)
    if h_edu == 2
        inf_grid[1:findall(age_grid .== age_inf)[]] .= 0.0
    end

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
        h_grid[1:4] .= h_grid[1]*0.6
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
    construct a type for mutable variables
    """
    V::Array{Float64,6}
    policy_a_p::Array{Float64,6}
    policy_x::Array{Float64,6}
    policy_l::Array{Float64,6}
    policy_K::Array{Float64,6}
end

function variables_function(parameters::NamedTuple)
    """
    construct a mutable object containing endogenous variables
    """
    # unpack parameters
    @unpack a_size, n_size, ϵ_size, ν_size, age_size = parameters

    # define value and policy functions: (a,n,ϵ,ν,f,t)
    V = zeros(a_size, n_size, ϵ_size, ν_size, 2, age_size)
    policy_a_p = ones(a_size, n_size, ϵ_size, ν_size, 2, age_size)
    policy_x = ones(a_size, n_size, ϵ_size, ν_size, 2, age_size)
    policy_l = ones(a_size, n_size, ϵ_size, ν_size, 2, age_size)
    policy_K = ones(a_size, n_size, ϵ_size, ν_size, 2, age_size)

    # return outputs
    variables = Mutable_Variables(V, policy_a_p, policy_x, policy_l, policy_K)
    return variables
end

function solve_value_and_policy_function!(variables::Mutable_Variables, parameters::NamedTuple)
    """
    compute value and policy functions
    """

    # unpack parameters
    @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
    @unpack ν_size, ν_grid, ν_Γ = parameters
    @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
    @unpack n_size, n_grid, n_Γ, n_max = parameters
    @unpack a_size, a_grid = parameters
    @unpack l_size, l_grid, x_size, x_grid = parameters
    @unpack inf_grid = parameters
    @unpack h_grid = parameters
    @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

    # loop over all states
    for age_i = age_size:(-1):1
        age = age_grid[age_i]
        age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
        age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
        println("Solving the problem of HH at age $age...")
        if age == age_max # terminal condition
            for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, a_i = 1:a_size
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                a = a_grid[a_i]
                w_bar = exp(h + ϵ + ν) * b
                variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = utility_function((1.0 + r) * a + w_bar, 0.0, 0.0, γ, ψ, κ, q_bar)
            end
        elseif age_ret < age < age_max # after retirement
            for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, a_i = 1:a_size
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                a = a_grid[a_i]
                w_bar = exp(h + ϵ + ν) * b
                EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
                V_all = utility_function.((1.0 + r) * a + w_bar .- a_grid, 0.0, 0.0, γ, ψ, κ, q_bar) .+ β * EV
                V_max_i = argmax(V_all)
                variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
                variables.policy_a_p[a_i, 1, ϵ_i, ν_i, 2, age_i] = V_max_i # a_grid[V_max_i]
            end
        elseif age == age_ret # at retirement age
            for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                n = n_grid[n_i]
                a = a_grid[a_i]
                w = exp(h + ϵ + ν)
                if n == 0
                    EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
                    V_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV
                    V_max_i = argmax(V_all)
                    variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
                    variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_max_i # a_grid[V_max_i]
                else
                    EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
                    Sol_all = zeros(a_size * x_size * l_size, 4)
                    Sol_all_i = 0
                    for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
                        Sol_all_i += 1
                        Sol_all[Sol_all_i, 1] = a_p_i
                        Sol_all[Sol_all_i, 2] = x_i
                        Sol_all[Sol_all_i, 3] = l_i
                        a_p = a_grid[a_p_i]
                        x = x_grid[x_i]
                        l = l_grid[l_i]
                        q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                        Sol_all[Sol_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV[a_p_i]
                    end
                    V_max_i = argmax(Sol_all[:, 4])
                    variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Sol_all[V_max_i, 4]
                    variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 1]) # a_grid[Int(Sol_all[V_max_i, 1])]
                    variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 2]) # x_grid[Int(Sol_all[V_max_i, 2])]
                    variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 3]) # l_grid[Int(Sol_all[V_max_i, 3])]
                end
            end
        elseif age_inf < age < age_ret # berween infertile age and retirement age
            for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
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
                    variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
                    variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_max_i # a_grid[V_max_i]
                else
                    EV = zeros(a_size)
                    for n_p_i = 1:n_size, ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                        EV += n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                    end
                    Sol_all = zeros(a_size * x_size * l_size, 4)
                    Sol_all_i = 0
                    for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
                        Sol_all_i += 1
                        Sol_all[Sol_all_i, 1] = a_p_i
                        Sol_all[Sol_all_i, 2] = x_i
                        Sol_all[Sol_all_i, 3] = l_i
                        a_p = a_grid[a_p_i]
                        x = x_grid[x_i]
                        l = l_grid[l_i]
                        q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                        Sol_all[Sol_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV[a_p_i]
                    end
                    V_max_i = argmax(Sol_all[:, 4])
                    variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Sol_all[V_max_i, 4]
                    variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 1]) # a_grid[Int(Sol_all[V_max_i, 1])]
                    variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 2]) # x_grid[Int(Sol_all[V_max_i, 2])]
                    variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 3]) # l_grid[Int(Sol_all[V_max_i, 3])]
                end
            end
        elseif age == age_inf # about to be infertile
            for f_i = 1:2, ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
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
                            variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_1_all[V_K_1_max_i]
                            variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_1_max_i # a_grid[V_K_1_max_i]
                            variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        else
                            variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_all[V_K_0_max_i]
                            variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_max_i # a_grid[V_K_0_max_i]
                        end
                    else
                        EV_K_0 = zeros(a_size)
                        for ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_0 += ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, 1, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        V_K_0_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV_K_0
                        V_K_0_max_i = argmax(V_K_0_all)
                        variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_all[V_K_0_max_i]
                        variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_max_i # a_grid[V_K_0_max_i]
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
                        Sol_K_0_all[Sol_K_0_all_i, 1] = a_p_i
                        Sol_K_0_all[Sol_K_0_all_i, 2] = x_i
                        Sol_K_0_all[Sol_K_0_all_i, 3] = l_i
                        a_p = a_grid[a_p_i]
                        x = x_grid[x_i]
                        l = l_grid[l_i]
                        q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                        Sol_K_0_all[Sol_K_0_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_0[a_p_i]
                    end
                    V_K_0_max_i = argmax(Sol_K_0_all[:, 4])
                    variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_0_all[V_K_0_max_i, 4]
                    variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 1]) # a_grid[Int(Sol_K_0_all[V_K_0_max_i, 1])]
                    variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 2]) # x_grid[Int(Sol_K_0_all[V_K_0_max_i, 2])]
                    variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 3]) # l_grid[Int(Sol_K_0_all[V_K_0_max_i, 3])]
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
                            Sol_K_0_all[Sol_K_0_all_i, 1] = a_p_i
                            Sol_K_0_all[Sol_K_0_all_i, 2] = x_i
                            Sol_K_0_all[Sol_K_0_all_i, 3] = l_i
                            a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            Sol_K_0_all[Sol_K_0_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_0[a_p_i]
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
                            Sol_K_1_all[Sol_K_1_all_i, 1] = a_p_i
                            Sol_K_1_all[Sol_K_1_all_i, 2] = x_i
                            Sol_K_1_all[Sol_K_1_all_i, 3] = l_i
                            a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            Sol_K_1_all[Sol_K_1_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_1[a_p_i]
                        end
                        V_K_1_max_i = argmax(Sol_K_1_all[:, 4])
                        if Sol_K_1_all[V_K_1_max_i, 4] >= Sol_K_0_all[V_K_0_max_i, 4]
                            variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_1_all[V_K_1_max_i, 4]
                            variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_1_all[V_K_1_max_i, 1]) # a_grid[Int(Sol_K_1_all[V_K_1_max_i, 1])]
                            variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_1_all[V_K_1_max_i, 2]) # x_grid[Int(Sol_K_1_all[V_K_1_max_i, 2])]
                            variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_1_all[V_K_1_max_i, 3]) # l_grid[Int(Sol_K_1_all[V_K_1_max_i, 3])]
                            variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        else
                            variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_0_all[V_K_0_max_i, 4]
                            variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 1]) # a_grid[Int(Sol_K_0_all[V_K_0_max_i, 1])]
                            variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 2]) # x_grid[Int(Sol_K_0_all[V_K_0_max_i, 2])]
                            variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 3]) # l_grid[Int(Sol_K_0_all[V_K_0_max_i, 3])]
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
                            Sol_K_0_all[Sol_K_0_all_i, 1] = a_p_i
                            Sol_K_0_all[Sol_K_0_all_i, 2] = x_i
                            Sol_K_0_all[Sol_K_0_all_i, 3] = l_i
                            a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            Sol_K_0_all[Sol_K_0_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_0[a_p_i]
                        end
                        V_K_0_max_i = argmax(Sol_K_0_all[:, 4])
                        variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_0_all[V_K_0_max_i, 4]
                        variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 1]) # a_grid[Int(Sol_K_0_all[V_K_0_max_i, 1])]
                        variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 2]) # x_grid[Int(Sol_K_0_all[V_K_0_max_i, 2])]
                        variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 3]) # l_grid[Int(Sol_K_0_all[V_K_0_max_i, 3])]
                    end
                end
            end
        else # fertile age
            for f_i = 1:2, ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
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
                            variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_1_all[V_K_1_max_i]
                            variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_1_max_i # a_grid[V_K_1_max_i]
                            variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        else
                            variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_all[V_K_0_max_i]
                            variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_max_i # a_grid[V_K_0_max_i]
                        end
                    else
                        EV_K_0 = zeros(a_size)
                        for ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_0 += ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        V_K_0_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV_K_0
                        V_K_0_max_i = argmax(V_K_0_all)
                        variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_all[V_K_0_max_i]
                        variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = V_K_0_max_i # a_grid[V_K_0_max_i]
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
                            Sol_K_0_all[Sol_K_0_all_i, 1] = a_p_i
                            Sol_K_0_all[Sol_K_0_all_i, 2] = x_i
                            Sol_K_0_all[Sol_K_0_all_i, 3] = l_i
                            a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            Sol_K_0_all[Sol_K_0_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_0[a_p_i]
                        end
                        V_K_0_max_i = argmax(Sol_K_0_all[:, 4])
                        variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_0_all[V_K_0_max_i, 4]
                        variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 1]) # a_grid[Int(Sol_K_0_all[V_K_0_max_i, 1])]
                        variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 2]) # x_grid[Int(Sol_K_0_all[V_K_0_max_i, 2])]
                        variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 3]) # l_grid[Int(Sol_K_0_all[V_K_0_max_i, 3])]
                    else
                        EV_K_0 = zeros(a_size)
                        for n_p_i = 1:n_size, ϵ_p_i = 1:ϵ_size, ν_p_i = 1:ν_size
                            EV_K_0 += n_Γ[n_i, n_p_i] * ϵ_Γ[ϵ_i, ϵ_p_i] * ν_Γ[ν_p_i] * variables.V[:, n_p_i, ϵ_p_i, ν_p_i, 2, age_i+1]
                        end
                        Sol_K_0_all = zeros(a_size * x_size * l_size, 4)
                        Sol_K_0_all_i = 0
                        for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
                            Sol_K_0_all_i += 1
                            Sol_K_0_all[Sol_K_0_all_i, 1] = a_p_i
                            Sol_K_0_all[Sol_K_0_all_i, 2] = x_i
                            Sol_K_0_all[Sol_K_0_all_i, 3] = l_i
                            a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            Sol_K_0_all[Sol_K_0_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_0[a_p_i]
                        end
                        V_K_0_max_i = argmax(Sol_K_0_all[:, 4])
                        variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_0_all[V_K_0_max_i, 4]
                        variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 1]) # a_grid[Int(Sol_K_0_all[V_K_0_max_i, 1])]
                        variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 2]) # x_grid[Int(Sol_K_0_all[V_K_0_max_i, 2])]
                        variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 3]) # l_grid[Int(Sol_K_0_all[V_K_0_max_i, 3])]
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
                            Sol_K_0_all[Sol_K_0_all_i, 1] = a_p_i
                            Sol_K_0_all[Sol_K_0_all_i, 2] = x_i
                            Sol_K_0_all[Sol_K_0_all_i, 3] = l_i
                            a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            Sol_K_0_all[Sol_K_0_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_0[a_p_i]
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
                            Sol_K_1_all[Sol_K_1_all_i, 1] = a_p_i
                            Sol_K_1_all[Sol_K_1_all_i, 2] = x_i
                            Sol_K_1_all[Sol_K_1_all_i, 3] = l_i
                            a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            Sol_K_1_all[Sol_K_1_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_1[a_p_i]
                        end
                        V_K_1_max_i = argmax(Sol_K_1_all[:, 4])
                        if Sol_K_1_all[V_K_1_max_i, 4] >= Sol_K_0_all[V_K_0_max_i, 4]
                            variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_1_all[V_K_1_max_i, 4]
                            variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_1_all[V_K_1_max_i, 1]) # a_grid[Int(Sol_K_1_all[V_K_1_max_i, 1])]
                            variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_1_all[V_K_1_max_i, 2]) # x_grid[Int(Sol_K_1_all[V_K_1_max_i, 2])]
                            variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_1_all[V_K_1_max_i, 3]) # l_grid[Int(Sol_K_1_all[V_K_1_max_i, 3])]
                            variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = 2
                        else
                            variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_0_all[V_K_0_max_i, 4]
                            variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 1]) # a_grid[Int(Sol_K_0_all[V_K_0_max_i, 1])]
                            variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 2]) # x_grid[Int(Sol_K_0_all[V_K_0_max_i, 2])]
                            variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 3])  # l_grid[Int(Sol_K_0_all[V_K_0_max_i, 3])]
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
                            Sol_K_0_all[Sol_K_0_all_i, 1] = a_p_i
                            Sol_K_0_all[Sol_K_0_all_i, 2] = x_i
                            Sol_K_0_all[Sol_K_0_all_i, 3] = l_i
                            a_p = a_grid[a_p_i]
                            x = x_grid[x_i]
                            l = l_grid[l_i]
                            q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
                            Sol_K_0_all[Sol_K_0_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV_K_0[a_p_i]
                        end
                        V_K_0_max_i = argmax(Sol_K_0_all[:, 4])
                        variables.V[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Sol_K_0_all[V_K_0_max_i, 4]
                        variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 1]) # a_grid[Int(Sol_K_0_all[V_K_0_max_i, 1])]
                        variables.policy_x[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 2]) # x_grid[Int(Sol_K_0_all[V_K_0_max_i, 2])]
                        variables.policy_l[a_i, n_i, ϵ_i, ν_i, f_i, age_i] = Int(Sol_K_0_all[V_K_0_max_i, 3])  # l_grid[Int(Sol_K_0_all[V_K_0_max_i, 3])]
                    end
                end
            end
        end
    end
end

#==============================#
# solve stationary equilibrium #
#==============================#
parameters = parameters_function()
variables = variables_function(parameters)
solve_value_and_policy_function!(variables, parameters)
V = variables.V
policy_a_p = variables.policy_a_p
policy_x = variables.policy_x
policy_l = variables.policy_l
policy_K = variables.policy_K
@save "workspace.jld2" parameters V policy_a_p policy_x policy_l policy_K

parameters_h_edu = parameters_function(h_edu = 1)
variables_h_edu = variables_function(parameters_h_edu)
solve_value_and_policy_function!(variables_h_edu, parameters_h_edu)
V_h_edu = variables_h_edu.V
policy_a_p_h_edu = variables_h_edu.policy_a_p
policy_x_h_edu = variables_h_edu.policy_x
policy_l_h_edu = variables_h_edu.policy_l
policy_K_h_edu = variables_h_edu.policy_K
@save "workspace_h_edu.jld2" parameters_h_edu V_h_edu policy_a_p_h_edu policy_x_h_edu policy_l_h_edu policy_K_h_edu

parameters_no_inf_risk = parameters_function(h_edu = 2)
variables_no_inf_risk = variables_function(parameters_no_inf_risk)
solve_value_and_policy_function!(variables_no_inf_risk, parameters_no_inf_risk)
V_no_inf_risk = variables_no_inf_risk.V
policy_a_p_no_inf_risk = variables_no_inf_risk.policy_a_p
policy_x_no_inf_risk = variables_no_inf_risk.policy_x
policy_l_no_inf_risk = variables_no_inf_risk.policy_l
policy_K_no_inf_risk = variables_no_inf_risk.policy_K
@save "workspace_no_inf_risk.jld2" parameters_no_inf_risk V_no_inf_risk policy_a_p_no_inf_risk policy_x_no_inf_risk policy_l_no_inf_risk policy_K_no_inf_risk

#===========#
# simuation #
#===========#
# load workspace
@load "workspace.jld2" parameters V policy_a_p policy_x policy_l policy_K
@load "workspace_h_edu.jld2" parameters_h_edu V_h_edu policy_a_p_h_edu policy_x_h_edu policy_l_h_edu policy_K_h_edu
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

panel_a_h_edu = ones(Int, num_hh, num_periods)
panel_a_p_h_edu = ones(Int, num_hh, num_periods)
panel_x_h_edu = ones(num_hh, num_periods)
panel_l_h_edu = ones(Int, num_hh, num_periods)
panel_n_h_edu = ones(Int, num_hh, num_periods)
panel_K_h_edu = ones(Int, num_hh, num_periods)

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

shock_n_h_edu = zeros(Int, num_hh, num_periods)
shock_ϵ_h_edu = zeros(Int, num_hh, num_periods)
shock_ν_h_edu = zeros(Int, num_hh, num_periods)
shock_f_h_edu = zeros(Int, num_hh, num_periods)

shock_n_no_inf_risk = zeros(Int, num_hh, num_periods)
shock_ϵ_no_inf_risk = zeros(Int, num_hh, num_periods)
shock_ν_no_inf_risk = zeros(Int, num_hh, num_periods)
shock_f_no_inf_risk = zeros(Int, num_hh, num_periods)

# Loop over HHs and Time periods
for period_i in 1:num_periods
    for hh_i in 1:num_hh
        if period_i == 1
            if hh_i == 1
                println("period = $period_i (1)")
            end
            # initiate states
            panel_a[hh_i, period_i] = 1
            panel_n[hh_i, period_i] = 1
            shock_ϵ[hh_i, period_i] = rand(Categorical(vec(parameters.ϵ_G)))
            shock_ν[hh_i, period_i] = rand(Categorical(vec(parameters.ν_G)))
            shock_f[hh_i, period_i] = rand(Categorical(vec([1.0 - parameters.inf_grid[period_i], parameters.inf_grid[period_i]])))

            panel_a_h_edu[hh_i, period_i] = 1
            panel_n_h_edu[hh_i, period_i] = 1
            shock_ϵ_h_edu[hh_i, period_i] = shock_ϵ[hh_i, period_i]
            shock_ν_h_edu[hh_i, period_i] = shock_ν[hh_i, period_i]
            shock_f_h_edu[hh_i, period_i] = shock_f[hh_i, period_i]

            panel_a_no_inf_risk[hh_i, period_i] = 1
            panel_n_no_inf_risk[hh_i, period_i] = 1
            shock_ϵ_no_inf_risk[hh_i, period_i] = shock_ϵ[hh_i, period_i]
            shock_ν_no_inf_risk[hh_i, period_i] = shock_ν[hh_i, period_i]
            shock_f_no_inf_risk[hh_i, period_i] = 1

            # actions
            panel_a_p[hh_i, period_i] = policy_a_p[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
            panel_x[hh_i, period_i] = policy_x[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
            panel_l[hh_i, period_i] = policy_l[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
            panel_K[hh_i, period_i] = policy_K[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]

            panel_a_p_h_edu[hh_i, period_i] = policy_a_p_h_edu[panel_a_h_edu[hh_i, period_i], panel_n_h_edu[hh_i, period_i], shock_ϵ_h_edu[hh_i, period_i], shock_ν_h_edu[hh_i, period_i], shock_f_h_edu[hh_i, period_i], period_i]
            panel_x_h_edu[hh_i, period_i] = policy_x_h_edu[panel_a_h_edu[hh_i, period_i], panel_n_h_edu[hh_i, period_i], shock_ϵ_h_edu[hh_i, period_i], shock_ν_h_edu[hh_i, period_i], shock_f_h_edu[hh_i, period_i], period_i]
            panel_l_h_edu[hh_i, period_i] = policy_l_h_edu[panel_a_h_edu[hh_i, period_i], panel_n_h_edu[hh_i, period_i], shock_ϵ_h_edu[hh_i, period_i], shock_ν_h_edu[hh_i, period_i], shock_f_h_edu[hh_i, period_i], period_i]
            panel_K_h_edu[hh_i, period_i] = policy_K_h_edu[panel_a_h_edu[hh_i, period_i], panel_n_h_edu[hh_i, period_i], shock_ϵ_h_edu[hh_i, period_i], shock_ν_h_edu[hh_i, period_i], shock_f_h_edu[hh_i, period_i], period_i]

            panel_a_p_no_inf_risk[hh_i, period_i] = policy_a_p_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
            panel_x_no_inf_risk[hh_i, period_i] = policy_x_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
            panel_l_no_inf_risk[hh_i, period_i] = policy_l_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
            panel_K_no_inf_risk[hh_i, period_i] = policy_K_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
        elseif 1 < period_i < findall(parameters.age_grid .== (parameters.age_inf + 1))[]
            if hh_i == 1
                println("period = $period_i (2)")
            end
            # initiate states
            panel_a[hh_i, period_i] = panel_a_p[hh_i, period_i-1]
            panel_n[hh_i, period_i] = rand(Categorical(vec(parameters.n_Γ[(panel_n[hh_i, period_i-1]+panel_K[hh_i, period_i-1]-1),:])))
            shock_ϵ[hh_i, period_i] = rand(Categorical(vec(parameters.ϵ_Γ[shock_ϵ[hh_i,period_i-1],:])))
            shock_ν[hh_i, period_i] = rand(Categorical(vec(parameters.ν_Γ)))
            shock_f[hh_i, period_i] = shock_f[hh_i, period_i-1] == 2 ? 2 : rand(Categorical(vec([1.0 - parameters.inf_grid[period_i], parameters.inf_grid[period_i]])))

            panel_a_h_edu[hh_i, period_i] = panel_a_p_h_edu[hh_i, period_i-1]
            panel_n_h_edu[hh_i, period_i] = rand(Categorical(vec(parameters_h_edu.n_Γ[(panel_n_h_edu[hh_i, period_i-1]+panel_K_h_edu[hh_i, period_i-1]-1),:])))
            shock_ϵ_h_edu[hh_i, period_i] = shock_ϵ[hh_i, period_i]
            shock_ν_h_edu[hh_i, period_i] = shock_ν[hh_i, period_i]
            shock_f_h_edu[hh_i, period_i] = shock_f[hh_i, period_i]

            panel_a_no_inf_risk[hh_i, period_i] = panel_a_p_no_inf_risk[hh_i, period_i-1]
            panel_n_no_inf_risk[hh_i, period_i] = rand(Categorical(vec(parameters_no_inf_risk.n_Γ[(panel_n_no_inf_risk[hh_i, period_i-1]+panel_K_no_inf_risk[hh_i, period_i-1]-1),:])))
            shock_ϵ_no_inf_risk[hh_i, period_i] = shock_ϵ[hh_i, period_i]
            shock_ν_no_inf_risk[hh_i, period_i] = shock_ν[hh_i, period_i]
            shock_f_no_inf_risk[hh_i, period_i] = 1

            # actions
            panel_a_p[hh_i, period_i] = policy_a_p[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
            panel_x[hh_i, period_i] = policy_x[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
            panel_l[hh_i, period_i] = policy_l[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
            panel_K[hh_i, period_i] = policy_K[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]

            panel_a_p_h_edu[hh_i, period_i] = policy_a_p_h_edu[panel_a_h_edu[hh_i, period_i], panel_n_h_edu[hh_i, period_i], shock_ϵ_h_edu[hh_i, period_i], shock_ν_h_edu[hh_i, period_i], shock_f_h_edu[hh_i, period_i], period_i]
            panel_x_h_edu[hh_i, period_i] = policy_x_h_edu[panel_a_h_edu[hh_i, period_i], panel_n_h_edu[hh_i, period_i], shock_ϵ_h_edu[hh_i, period_i], shock_ν_h_edu[hh_i, period_i], shock_f_h_edu[hh_i, period_i], period_i]
            panel_l_h_edu[hh_i, period_i] = policy_l_h_edu[panel_a_h_edu[hh_i, period_i], panel_n_h_edu[hh_i, period_i], shock_ϵ_h_edu[hh_i, period_i], shock_ν_h_edu[hh_i, period_i], shock_f_h_edu[hh_i, period_i], period_i]
            panel_K_h_edu[hh_i, period_i] = policy_K_h_edu[panel_a_h_edu[hh_i, period_i], panel_n_h_edu[hh_i, period_i], shock_ϵ_h_edu[hh_i, period_i], shock_ν_h_edu[hh_i, period_i], shock_f_h_edu[hh_i, period_i], period_i]

            panel_a_p_no_inf_risk[hh_i, period_i] = policy_a_p_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
            panel_x_no_inf_risk[hh_i, period_i] = policy_x_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
            panel_l_no_inf_risk[hh_i, period_i] = policy_l_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
            panel_K_no_inf_risk[hh_i, period_i] = policy_K_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
        elseif findall(parameters.age_grid .== (parameters.age_inf + 1))[] <= period_i <= findall(parameters.age_grid .== parameters.age_ret)[] 
            if hh_i == 1
                println("period = $period_i (3)")
            end
            # initiate states
            panel_a[hh_i, period_i] = panel_a_p[hh_i, period_i-1]
            panel_n[hh_i, period_i] = rand(Categorical(vec(parameters.n_Γ[(panel_n[hh_i, period_i-1]+panel_K[hh_i, period_i-1]-1),:])))
            shock_ϵ[hh_i, period_i] = rand(Categorical(vec(parameters.ϵ_Γ[shock_ϵ[hh_i,period_i-1],:])))
            shock_ν[hh_i, period_i] = rand(Categorical(vec(parameters.ν_Γ)))
            shock_f[hh_i, period_i] = 2

            panel_a_h_edu[hh_i, period_i] = panel_a_p_h_edu[hh_i, period_i-1]
            panel_n_h_edu[hh_i, period_i] = rand(Categorical(vec(parameters_h_edu.n_Γ[(panel_n_h_edu[hh_i, period_i-1]+panel_K_h_edu[hh_i, period_i-1]-1),:])))
            shock_ϵ_h_edu[hh_i, period_i] = shock_ϵ[hh_i, period_i]
            shock_ν_h_edu[hh_i, period_i] = shock_ν[hh_i, period_i]
            shock_f_h_edu[hh_i, period_i] = shock_f[hh_i, period_i]

            panel_a_no_inf_risk[hh_i, period_i] = panel_a_p_no_inf_risk[hh_i, period_i-1]
            panel_n_no_inf_risk[hh_i, period_i] = rand(Categorical(vec(parameters_no_inf_risk.n_Γ[(panel_n_no_inf_risk[hh_i, period_i-1]+panel_K_no_inf_risk[hh_i, period_i-1]-1),:])))
            shock_ϵ_no_inf_risk[hh_i, period_i] = shock_ϵ[hh_i, period_i]
            shock_ν_no_inf_risk[hh_i, period_i] = shock_ν[hh_i, period_i]
            shock_f_no_inf_risk[hh_i, period_i] = 2

            # actions
            panel_a_p[hh_i, period_i] = policy_a_p[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
            panel_x[hh_i, period_i] = policy_x[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
            panel_l[hh_i, period_i] = policy_l[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
            panel_K[hh_i, period_i] = policy_K[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]

            panel_a_p_h_edu[hh_i, period_i] = policy_a_p_h_edu[panel_a_h_edu[hh_i, period_i], panel_n_h_edu[hh_i, period_i], shock_ϵ_h_edu[hh_i, period_i], shock_ν_h_edu[hh_i, period_i], shock_f_h_edu[hh_i, period_i], period_i]
            panel_x_h_edu[hh_i, period_i] = policy_x_h_edu[panel_a_h_edu[hh_i, period_i], panel_n_h_edu[hh_i, period_i], shock_ϵ_h_edu[hh_i, period_i], shock_ν_h_edu[hh_i, period_i], shock_f_h_edu[hh_i, period_i], period_i]
            panel_l_h_edu[hh_i, period_i] = policy_l_h_edu[panel_a_h_edu[hh_i, period_i], panel_n_h_edu[hh_i, period_i], shock_ϵ_h_edu[hh_i, period_i], shock_ν_h_edu[hh_i, period_i], shock_f_h_edu[hh_i, period_i], period_i]
            panel_K_h_edu[hh_i, period_i] = policy_K_h_edu[panel_a_h_edu[hh_i, period_i], panel_n_h_edu[hh_i, period_i], shock_ϵ_h_edu[hh_i, period_i], shock_ν_h_edu[hh_i, period_i], shock_f_h_edu[hh_i, period_i], period_i]

            panel_a_p_no_inf_risk[hh_i, period_i] = policy_a_p_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
            panel_x_no_inf_risk[hh_i, period_i] = policy_x_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
            panel_l_no_inf_risk[hh_i, period_i] = policy_l_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
            panel_K_no_inf_risk[hh_i, period_i] = policy_K_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
        elseif findall(parameters.age_grid .== (parameters.age_ret))[] < period_i < parameters.age_size 
            if hh_i == 1
                println("period = $period_i (4)")
            end
            # initiate states
            panel_a[hh_i, period_i] = panel_a_p[hh_i, period_i-1]
            panel_n[hh_i, period_i] = 1
            shock_ϵ[hh_i, period_i] = shock_ϵ[hh_i, period_i-1]
            shock_ν[hh_i, period_i] = shock_ν[hh_i, period_i-1]
            shock_f[hh_i, period_i] = 2

            panel_a_h_edu[hh_i, period_i] = panel_a_p_h_edu[hh_i, period_i-1]
            panel_n_h_edu[hh_i, period_i] = 1
            shock_ϵ_h_edu[hh_i, period_i] = shock_ϵ[hh_i, period_i]
            shock_ν_h_edu[hh_i, period_i] = shock_ν[hh_i, period_i]
            shock_f_h_edu[hh_i, period_i] = shock_f[hh_i, period_i]

            panel_a_no_inf_risk[hh_i, period_i] = panel_a_p_no_inf_risk[hh_i, period_i-1]
            panel_n_no_inf_risk[hh_i, period_i] = 1
            shock_ϵ_no_inf_risk[hh_i, period_i] = shock_ϵ[hh_i, period_i]
            shock_ν_no_inf_risk[hh_i, period_i] = shock_ν[hh_i, period_i]
            shock_f_no_inf_risk[hh_i, period_i] = 2

            # actions
            panel_a_p[hh_i, period_i] = policy_a_p[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]

            panel_a_p_h_edu[hh_i, period_i] = policy_a_p_h_edu[panel_a_h_edu[hh_i, period_i], panel_n_h_edu[hh_i, period_i], shock_ϵ_h_edu[hh_i, period_i], shock_ν_h_edu[hh_i, period_i], shock_f_h_edu[hh_i, period_i], period_i]

            panel_a_p_no_inf_risk[hh_i, period_i] = policy_a_p_no_inf_risk[panel_a_no_inf_risk[hh_i, period_i], panel_n_no_inf_risk[hh_i, period_i], shock_ϵ_no_inf_risk[hh_i, period_i], shock_ν_no_inf_risk[hh_i, period_i], shock_f_no_inf_risk[hh_i, period_i], period_i]
        else
            if hh_i == 1
                println("period = $period_i (5)")
            end
            # initiate states
            panel_a[hh_i, period_i] = panel_a_p[hh_i, period_i-1]
            panel_n[hh_i, period_i] = 1
            shock_ϵ[hh_i, period_i] = shock_ϵ[hh_i, period_i-1]
            shock_ν[hh_i, period_i] = shock_ν[hh_i, period_i-1]
            shock_f[hh_i, period_i] = 2

            panel_a_h_edu[hh_i, period_i] = panel_a_p_h_edu[hh_i, period_i-1]
            panel_n_h_edu[hh_i, period_i] = 1
            shock_ϵ_h_edu[hh_i, period_i] = shock_ϵ[hh_i, period_i]
            shock_ν_h_edu[hh_i, period_i] = shock_ν[hh_i, period_i]
            shock_f_h_edu[hh_i, period_i] = shock_f[hh_i, period_i]

            panel_a_no_inf_risk[hh_i, period_i] = panel_a_p_no_inf_risk[hh_i, period_i-1]
            panel_n_no_inf_risk[hh_i, period_i] = 1
            shock_ϵ_no_inf_risk[hh_i, period_i] = shock_ϵ[hh_i, period_i]
            shock_ν_no_inf_risk[hh_i, period_i] = shock_ν[hh_i, period_i]
            shock_f_no_inf_risk[hh_i, period_i] = 2
            # actions
            panel_a_p[hh_i, period_i] = 1
            panel_a_p_h_edu[hh_i, period_i] = 1
            panel_a_p_no_inf_risk[hh_i, period_i] = 1
        end
    end
end

#=================================#
# simulation results (adjusted h) #
#=================================#
# life-cyle earning income
# plot(parameters.age_min:parameters.age_ret, parameters.h_grid[1:(parameters.age_ret-parameters.age_min+1)], legend=:topleft)
# plot!(parameters.age_min:parameters.age_ret, parameters_h_edu.h_grid[1:(parameters.age_ret-parameters.age_min+1)], legend=:topleft)
plot_h = Plots.plot(
    box=:on,
    size=[800, 500],
    xlim=[18, 65],
    xticks=18:5:65,
    ylim=[0.8, 4.2],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
plot_h = Plots.plot!(
    parameters.age_min:parameters.age_ret,
    parameters.h_grid[1:(parameters.age_ret-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:blue
)
savefig(plot_h, string("plot_h.pdf"))
plot_h_edu = plot(
    box=:on,
    size=[800, 500],
    xlim=[18, 65],
    xticks=18:5:65,
    ylim=[0.8, 4.2],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
plot_h_edu = plot!(
    parameters.age_min:parameters.age_ret,
    parameters_h_edu.h_grid[1:(parameters.age_ret-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:red,
    ls=:dashdot
)
savefig(plot_h_edu, string("plot_h_edu.pdf"))
plot_h_edu_mixed = plot(
    box=:on,
    size=[800, 500],
    xlim=[18, 65],
    xticks=18:5:65,
    ylim=[0.8, 4.2],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
plot_h_edu_mixed = plot!(
    parameters.age_min:parameters.age_ret,
    parameters.h_grid[1:(parameters.age_ret-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:blue
)
plot_h_edu_mixed = plot!(
    parameters.age_min:parameters.age_ret,
    parameters_h_edu.h_grid[1:(parameters.age_ret-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:red,
    ls=:dashdot
)
savefig(plot_h_edu_mixed, string("plot_h_edu_mixed.pdf"))

# average conception rate
avg_conception_rate = zeros(parameters.age_size)
avg_conception_rate_h_edu = zeros(parameters.age_size)
avg_conception_rate_no_inf_risk = zeros(parameters.age_size)
for t=1:parameters.age_size
    avg_conception_rate[t] = sum(panel_K[:,t].-1) / (num_hh/1000)
    avg_conception_rate_h_edu[t] = sum(panel_K_h_edu[:,t].-1) / (num_hh/1000) 
    avg_conception_rate_no_inf_risk[t] = sum(panel_K_no_inf_risk[:,t].-1) / (num_hh/1000) 
end
# plot(parameters.age_min:parameters.age_inf, avg_conception_rate[1:(parameters.age_inf-parameters.age_min+1)], legend=:topright)
# plot!(parameters.age_min:parameters.age_inf, avg_conception_rate_h_edu[1:(parameters.age_inf-parameters.age_min+1)], legend=:topright)
plot_conception_dist_by_age = plot(
    box=:on,
    size=[800, 500],
    xlim=[18, 45],
    xticks=18:3:45,
    ylim=[-10,160],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
plot_conception_dist_by_age = plot!(
    parameters.age_min:parameters.age_inf,
    avg_conception_rate[1:(parameters.age_inf-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:blue
)
savefig(plot_conception_dist_by_age, string("plot_conception_dist_by_age.pdf"))
plot_conception_dist_by_age_h_edu = plot(
    box=:on,
    size=[800, 500],
    xlim=[18, 45],
    xticks=18:3:45,
    ylim=[-10,160],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
plot_conception_dist_by_age_h_edu = plot!(
    parameters.age_min:parameters.age_inf,
    avg_conception_rate_h_edu[1:(parameters.age_inf-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:red,
    ls=:dashdot
)
savefig(plot_conception_dist_by_age_h_edu, string("plot_conception_dist_by_age_h_edu.pdf"))
plot_conception_dist_by_age_mixed = plot(
    box=:on,
    size=[800, 500],
    xlim=[18, 45],
    xticks=18:3:45,
    ylim=[-10,160],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
plot_conception_dist_by_age_mixed = plot!(
    parameters.age_min:parameters.age_inf,
    avg_conception_rate[1:(parameters.age_inf-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:blue
)
plot_conception_dist_by_age_mixed = plot!(
    parameters.age_min:parameters.age_inf,
    avg_conception_rate_h_edu[1:(parameters.age_inf-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:red,
    ls=:dashdot
)
savefig(plot_conception_dist_by_age_mixed, string("plot_conception_dist_by_age_mixed.pdf"))

# average cumulative number of kids at home
avg_cum_number_of_kids = zeros(parameters.age_size)
avg_cum_number_of_kids_h_edu = zeros(parameters.age_size)
avg_cum_number_of_kids_no_inf_risk = zeros(parameters.age_size)
for t=1:parameters.age_size
    avg_cum_number_of_kids[t] = mean(panel_n[:,t].-1)
    avg_cum_number_of_kids_h_edu[t] = mean(panel_n_h_edu[:,t].-1)
    avg_cum_number_of_kids_no_inf_risk[t] = mean(panel_n_no_inf_risk[:,t].-1)
end
plot(parameters.age_min:parameters.age_inf, avg_cum_number_of_kids[1:(parameters.age_inf-parameters.age_min+1)], legend=:topright)
plot!(parameters.age_min:parameters.age_inf, avg_cum_number_of_kids_h_edu[1:(parameters.age_inf-parameters.age_min+1)], legend=:topright)

# cumulative number of children born
avg_cum_number_of_born = zeros(parameters.age_size)
avg_cum_number_of_born_h_edu = zeros(parameters.age_size)
avg_cum_number_of_born_no_inf_risk = zeros(parameters.age_size)
for t=1:parameters.age_size
    avg_cum_number_of_born[t] = sum(panel_K[:,1:t].-1) / num_hh
    avg_cum_number_of_born_h_edu[t] = sum(panel_K_h_edu[:,1:t].-1) / num_hh
    avg_cum_number_of_born_no_inf_risk[t] = sum(panel_K_no_inf_risk[:,1:t].-1) / num_hh
end
# plot(parameters.age_min:parameters.age_inf, avg_cum_number_of_born[1:(parameters.age_inf-parameters.age_min+1)], legend=:topright)
# plot!(parameters.age_min:parameters.age_inf, avg_cum_number_of_born_h_edu[1:(parameters.age_inf-parameters.age_min+1)], legend=:topright)
avg_cum_number_of_born_by_age = plot(
    box=:on,
    size=[800, 500],
    xlim=[18, 45],
    xticks=18:3:45,
    ylim=[-0.05,1.25],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
avg_cum_number_of_born_by_age = plot!(
    parameters.age_min:parameters.age_inf,
    avg_cum_number_of_born[1:(parameters.age_inf-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:blue
)
savefig(avg_cum_number_of_born_by_age, string("avg_cum_number_of_born_by_age.pdf"))
avg_cum_number_of_born_by_age_h_edu = plot(
    box=:on,
    size=[800, 500],
    xlim=[18, 45],
    xticks=18:3:45,
    ylim=[-0.05,1.25],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
avg_cum_number_of_born_by_age_h_edu = plot!(
    parameters.age_min:parameters.age_inf,
    avg_cum_number_of_born_h_edu[1:(parameters.age_inf-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:red,
    ls=:dashdot
)
savefig(avg_cum_number_of_born_by_age_h_edu, string("avg_cum_number_of_born_by_age_h_edu.pdf"))
avg_cum_number_of_born_by_age_mixed = plot(
    box=:on,
    size=[800, 500],
    xlim=[18, 45],
    xticks=18:3:45,
    ylim=[-0.05,1.25],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
avg_cum_number_of_born_by_age_mixed = plot!(
    parameters.age_min:parameters.age_inf,
    avg_cum_number_of_born[1:(parameters.age_inf-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:blue
)
avg_cum_number_of_born_by_age_mixed = plot!(
    parameters.age_min:parameters.age_inf,
    avg_cum_number_of_born_h_edu[1:(parameters.age_inf-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:red,
    ls=:dashdot
)
savefig(avg_cum_number_of_born_by_age_mixed, string("avg_cum_number_of_born_by_age_mixed.pdf"))

#==========================================#
# simulation results (no infertility risk) #
#==========================================#
# infertility risk
plot_inf_risk_mixed = plot(
    box=:on,
    size=[800, 500],
    xlim=[18, 65],
    xticks=18:5:65,
    ylim=[-0.05, 1.05],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
plot_inf_risk_mixed = plot!(
    parameters.age_min:parameters.age_ret,
    parameters.inf_grid[1:(parameters.age_ret-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:blue
)
plot_inf_risk_mixed = plot!(
    parameters.age_min:parameters.age_ret,
    parameters_no_inf_risk.inf_grid[1:(parameters.age_ret-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:red,
    ls=:dashdot
)
savefig(plot_inf_risk_mixed, string("plot_inf_risk_mixed.pdf"))

# fertility distribution across ages
plot_conception_dist_by_age_mixed_no_inf = plot(
    box=:on,
    size=[800, 500],
    xlim=[18, 45],
    xticks=18:3:45,
    ylim=[-10,50],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
plot_conception_dist_by_age_mixed_no_inf = plot!(
    parameters.age_min:parameters.age_inf,
    avg_conception_rate[1:(parameters.age_inf-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:blue
)
plot_conception_dist_by_age_mixed_no_inf = plot!(
    parameters.age_min:parameters.age_inf,
    avg_conception_rate_no_inf_risk[1:(parameters.age_inf-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:red,
    ls=:dashdot
)
savefig(plot_conception_dist_by_age_mixed_no_inf, string("plot_conception_dist_by_age_mixed_no_inf.pdf"))
# cumulative kids born
avg_cum_number_of_born_by_age_mixed_no_inf_risk = plot(
    box=:on,
    size=[800, 500],
    xlim=[18, 45],
    xticks=18:3:45,
    ylim=[-0.05,6.05],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
avg_cum_number_of_born_by_age_mixed_no_inf_risk = plot!(
    parameters.age_min:parameters.age_inf,
    avg_cum_number_of_born[1:(parameters.age_inf-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:blue
)
avg_cum_number_of_born_by_age_mixed_no_inf_risk = plot!(
    parameters.age_min:parameters.age_inf,
    avg_cum_number_of_born_no_inf_risk[1:(parameters.age_inf-parameters.age_min+1)],
    label="",
    lw=3,
    lc=:red,
    ls=:dashdot
)
savefig(avg_cum_number_of_born_by_age_mixed_no_inf_risk, string("avg_cum_number_of_born_by_age_mixed_no_inf_risk.pdf"))

#======================#
# fertility elasticity #
#======================#
var_h = zeros(Real, num_hh, num_periods)
var_fertility_h = zeros(Real, num_hh, num_periods)
var_h_ = zeros(Real, num_hh, num_periods)
var_fertility_h_ = zeros(Real, num_hh, num_periods)
var_asset = zeros(Real, num_hh, num_periods)
var_fertility_asset = zeros(Real, num_hh, num_periods)
for period_i in 1:num_periods
    for hh_i in 1:num_hh
        # increased h
        h_old = period_i
        h_new = h_old == parameters.h_size ? h_old : period_i + 1
        var_h[hh_i, period_i] = parameters.h_grid[h_new] - parameters.h_grid[h_old]
        K_old = panel_K[hh_i, period_i] 
        K_new = policy_K[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], h_new]
        var_fertility_h[hh_i, period_i] = var_h[hh_i, period_i] == 0 ? 0 : K_new - K_old
        # decreased h
        h_old = period_i
        h_new = h_old == 1 ? h_old : period_i - 1
        var_h_[hh_i, period_i] = parameters.h_grid[h_new] - parameters.h_grid[h_old]
        K_old = panel_K[hh_i, period_i] 
        K_new = policy_K[panel_a[hh_i, period_i], panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], h_new]
        var_fertility_h_[hh_i, period_i] = var_h[hh_i, period_i] == 0 ? 0 : K_new - K_old
        # increased a
        a_old = panel_a[hh_i, period_i]
        a_new = a_old < parameters.a_size ? a_old + 1 : parameters.a_size
        var_asset[hh_i, period_i] = parameters.a_grid[a_new] - parameters.a_grid[a_old]
        # var_asset[hh_i, period_i] = a_old == 1 ? log(parameters.a_grid[a_new]) - log(parameters.a_grid[a_old] + eps()) : log(parameters.a_grid[a_new]) - log(parameters.a_grid[a_old])
        K_old = panel_K[hh_i, period_i] 
        K_new = policy_K[a_new, panel_n[hh_i, period_i], shock_ϵ[hh_i, period_i], shock_ν[hh_i, period_i], shock_f[hh_i, period_i], period_i]
        var_fertility_asset[hh_i, period_i] = var_asset[hh_i, period_i] == 0 ? 0 : K_new - K_old
    end
end
α_estimated_temp = zeros(Real, 2, num_hh, num_periods)
α_estimated = zeros(Real, 2, num_periods)
for period_i in 1:num_periods
    # model(var_h, α) = α[1] .+ α[2] .* var_h
    # α_int = [0.0, 0.1]
    # fit = curve_fit(model, var_h[:,period_i], var_fertility_h[:,period_i], α_int)
    # α_estimated[1, period_i] = fit.param[1] 
    # α_estimated[2, period_i] = fit.param[2]
    for hh_i in 1:num_hh
        α_estimated_temp[1, hh_i, period_i] = var_h[hh_i,period_i] != 0.0 ? var_fertility_h[hh_i,period_i] / var_h[hh_i,period_i] : 0.0
        α_estimated_temp[2, hh_i, period_i] = var_h_[hh_i,period_i] != 0.0 ? var_fertility_h_[hh_i,period_i] / var_h_[hh_i,period_i] : 0.0
    end
    α_estimated[1, period_i] = mean(α_estimated_temp[1,:,period_i])
    α_estimated[2, period_i] = mean(α_estimated_temp[2,:,period_i])
end
β_estimated = zeros(Real, 2, num_periods)
for period_i in 1:num_periods 
    model(var_a, β) = β[1] .+ β[2] .* var_a
    β_int = [0.0, 0.1]
    fit = curve_fit(model, var_asset[:,period_i], var_fertility_asset[:,period_i], β_int)
    β_estimated[1, period_i] = fit.param[1] 
    β_estimated[2, period_i] = fit.param[2] 
end
age_group = [19, 24, 29, 34, 39, 44, 50] 
age_label = ["18-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-50"]
α_estimated_grouped = zeros(Real, length(age_group))
α_estimated_grouped_ = zeros(Real, length(age_group))
β_estimated_grouped = zeros(Real, length(age_group))
for group_i = 1:length(age_group)
    if group_i == 1
        age_min_index = 1
        age_max_index = findall(parameters.age_grid .== age_group[group_i])[]
    else
        age_min_index = findall(parameters.age_grid .== age_group[group_i-1])[] + 1
        age_max_index = findall(parameters.age_grid .== age_group[group_i])[]
    end
    α_estimated_grouped[group_i] = sum(α_estimated[1, age_min_index:age_max_index])
    α_estimated_grouped_[group_i] = sum(α_estimated[2, age_min_index:age_max_index])
    β_estimated_grouped[group_i] = sum(β_estimated[2, age_min_index:age_max_index])
end
plot_fertility_sensitivity_α = plot(
    box=:on,
    size=[800, 500],
    ylims=[-0.05,0.90],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
plot_fertility_sensitivity_α = plot!(
    1:length(age_group),
    α_estimated_grouped,
    label="",
    lw=3,
    xticks=(1:length(age_group),age_label),
    lc=:blue
)
savefig(plot_fertility_sensitivity_α, string("plot_fertility_sensitivity_α.pdf"))
plot_fertility_sensitivity_α_ = plot(
    box=:on,
    size=[800, 500],
    ylims=[-0.05,0.90],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
plot_fertility_sensitivity_α_ = plot!(
    1:length(age_group),
    α_estimated_grouped_,
    label="",
    lw=3,
    xticks=(1:length(age_group),age_label),
    lc=:red,
    ls=:dashdot
)
savefig(plot_fertility_sensitivity_α_, string("plot_fertility_sensitivity_α_.pdf"))
plot_fertility_sensitivity_β = plot(
    box=:on,
    size=[800, 500],
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm
)
plot_fertility_sensitivity_β = plot!(
    1:length(age_group),
    β_estimated_grouped,
    label="",
    lw=3,
    xticks=(1:length(age_group),age_label),
    lc=:blue
)
savefig(plot_fertility_sensitivity_β, string("plot_fertility_sensitivity_β.pdf"))

#================#
# checking plots #
#================#
# infertility risk from 18 to 45
# plot(parameters.data_age, parameters.data_inf, seriestype=:scatter, legend=:topleft, label="data")
# plot!(parameters.age_min:parameters.age_inf, parameters.inf_grid[1:(parameters.age_inf-parameters.age_min+1)], label="exponential fit")
plot_inf_18_45 = plot(
    box=:on,
    size=[800, 500],
    xtickfont=font(12, "Computer Modern", :black),
    ytickfont=font(12, "Computer Modern", :black),
    legendfont=font(12, "Computer Modern", :black),
    guidefont=font(14, "Computer Modern", :black),
    titlefont=font(14, "Computer Modern", :black),
    margin=4mm
)
plot_inf_18_45 = plot!(
    parameters.age_min:parameters.age_inf,
    parameters.inf_grid[1:(parameters.age_inf-parameters.age_min+1)],
    label="model (exponential fit)",
    lw=3,
    lc=:blue
)
plot_inf_18_45 = scatter!(
    parameters.data_age,
    parameters.data_inf,
    ylim=[0.0, 1.05],
    legend=:topleft,
    label="data from Trussell and Wilson (1985, Population Studies)",
    ms=5,
    mc=:black,
    markershape=:rect,
    xlabel="Age",
    title="Infertility Risk (Age 18-45)"
)
savefig(plot_inf_18_45, string("plot_inf_18_45.pdf"))

# infertility risk from 18 to 80
# plot(parameters.age_grid, parameters.inf_grid, legend=:topleft, label="exponential fit")
plot_inf_18_80 = plot(
    box=:on,
    size=[800, 500],
    xtickfont=font(12, "Computer Modern", :black),
    ytickfont=font(12, "Computer Modern", :black),
    legendfont=font(12, "Computer Modern", :black),
    guidefont=font(14, "Computer Modern", :black),
    titlefont=font(14, "Computer Modern", :black),
    margin=4mm
)
plot_inf_18_80 = plot!(
    parameters.age_grid,
    parameters.inf_grid,
    label="model (exponential fit)",
    lw=3,
    lc=:blue
)
plot_inf_18_80 = scatter!(
    parameters.data_age,
    parameters.data_inf,
    ylim=[0.0, 1.05],
    legend=:bottomright,
    label="data from Trussell and Wilson (1985, Population Studies)",
    ms=5,
    mc=:black,
    markershape=:rect,
    xlabel="Age",
    title="Infertility Risk (Age 18-80)"
)
savefig(plot_inf_18_80, string("plot_inf_18_80.pdf"))

# life-cycle earning profile
# plot(parameters.age_min:parameters.age_ret, parameters.data_h[1:(parameters.age_ret-parameters.age_min+1)], legend=:topleft, label="data")
# plot!(parameters.age_min:parameters.age_ret, parameters.h_grid[1:(parameters.age_ret-parameters.age_min+1)], linestyle=:dash)
