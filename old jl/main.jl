using Distributions
using LsqFit
using Plots
using Parameters: @unpack
using NLopt
using QuantEcon: gridmake, rouwenhorst, tauchen, stationary_distributions, MarkovChain
using FLOWMath: Akima

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
    model_inf[findall(model_age .== age_inf)[]:end] .= 1.0
    return model_age, model_inf
end

# plot(data_age, data_inf, seriestype=:scatter, legend=:none)
# plot!(age_min:45, model_inf[1:(45-18+1)])
# plot(model_age,model_inf)

function utility_function(c::Real, n::Real, q::Real, γ::Real, ψ::Real, κ::Real)
    """
    utility function with child qquality
    """
    if c > 0.0
        return 1.0 / ((1.0 - γ) * c^(γ - 1.0)) + ψ / ((1.0 - κ) * (n * q)^(κ - 1.0))
    else
        return -10^9
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
    q_bar::Real=0.34,               # lower bound on children's consumption
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
    a_max::Real=100,                # max of asset holding
    a_size::Integer=50,             # number of asset
    a_degree::Integer=2             # curvature of asset gridpoints
)
    """
    contruct an immutable object containg all paramters
    """
    # infertility parameters: taken from Trussell and Wilson (1985, Population Studies)
    data_inf = [0.07, 0.131, 0.231, 0.345, 0.576, 0.952]
    data_age = [20, 25, 30, 35, 40, 45]
    age_grid, inf_grid = infertility_risk_function(data_age, data_inf, age_min, age_max, age_inf)
    age_size = length(age_grid)

    # transition of child dependence
    n_grid = collect(0:n_max)
    n_size = length(n_grid)
    n_Γ = binomial_matrix_function(n_max, p)

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
        n_max=n_max,
        n_size=n_size,
        n_grid=n_grid,
        n_Γ=n_Γ,
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
        a_degree=a_degree
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
end

function variables_function(parameters::NamedTuple)
    """
    construct a mutable object containing endogenous variables
    """
    # unpack parameters
    @unpack a_size, n_size, ϵ_size, ν_size, age_size = parameters

    # define value and policy functions: (a,n,ϵ,ν,f,t)
    V = zeros(a_size, n_size, ϵ_size, ν_size, 2, age_size)
    policy_a_p = zeros(a_size, n_size, ϵ_size, ν_size, 2, age_size)
    policy_x = zeros(a_size, n_size, ϵ_size, ν_size, 2, age_size)
    policy_l = zeros(a_size, n_size, ϵ_size, ν_size, 2, age_size)

    # return outputs
    variables = Mutable_Variables(V, policy_a_p, policy_x, policy_l)
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
    @unpack n_size, n_grid, n_Γ = parameters
    @unpack a_size, a_grid = parameters
    @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar = parameters

    # loop over all states
    for age_i = age_size:(-1):1
        age = age_grid[age_i]
        if age == age_max # terminal condition
            for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                n = n_grid[n_i]
                a = a_grid[a_i]
                w_bar = exp(ϵ + ν) * b
                variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = utility_function((1.0 + r) * a + w_bar, n, 0.0, γ, ψ, κ)
            end
        elseif age_ret < age < age_max # after retirement
            for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                n = n_grid[n_i]
                a = a_grid[a_i]
                w_bar = exp(ϵ + ν) * b
                EV = zeros(a_size)
                for n_p_i = 1:n_size
                    EV += n_Γ[n_i, n_p_i] * variables.V[:, n_p_i, ϵ_i, ν_i, 2, age_i+1]
                end
                EV_func_itp = Akima(a_grid, EV)
                function V_func(x::Vector, grad::Vector)
                    # x[1] = a_p 
                    return utility_function((1.0 + r) * a + w_bar - x[1], n, 0.0, γ, ψ, κ) + β * EV_func_itp(x[1])
                end
                function c_constraint(x::Vector, grad::Vector)
                    return x[1] + 10^(-16) - (1.0 + r) * a - w_bar
                end
                opt = Opt(:LN_COBYLA, 1)
                inequality_constraint!(opt, c_constraint)
                opt.lower_bounds = [0.0]
                opt.upper_bounds = [(1.0 + r) * a + w_bar]
                opt.xtol_rel = 1e-4
                opt.max_objective = V_func
                (maxf, maxx, ret) = optimize(opt, [a])
                variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = maxf
                variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = maxx[1]
            end
        elseif age == age_ret # at retirement age
            for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
                ν = ν_grid[ν_i]
                ϵ = ϵ_grid[ϵ_i]
                n = n_grid[n_i]
                a = a_grid[a_i]
                w = exp(ϵ + ν)
                if n == 0
                    EV = variables.V[:, n_i, ϵ_i, ν_i, 2, age_i+1]
                    EV_func_itp = Akima(a_grid, EV)
                    function V_func(x::Vector, grad::Vector)
                        # x[1] = a_p 
                        return utility_function((1.0 + r) * a + w - x[1], n, 0.0, γ, ψ, κ) + β * EV_func_itp(x[1])
                    end
                    function c_constraint(x::Vector, grad::Vector)
                        return x[1] + 10^(-16) - (1.0 + r) * a - w
                    end
                    opt = Opt(:LN_COBYLA, 1)
                    inequality_constraint!(opt, c_constraint)
                    opt.lower_bounds = [0.0]
                    opt.upper_bounds = [(1.0 + r) * a + w]
                    opt.xtol_rel = 1e-4
                    opt.max_objective = V_func
                    (maxf, maxx, ret) = optimize(opt, [a])
                    variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = maxf
                    variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = maxx[1]
                else
                    println("($a_i,$n_i,$ϵ_i,$ν_i,$age_i)")
                    EV = zeros(a_size)
                    for n_p_i = 1:n_size
                        EV += n_Γ[n_i, n_p_i] * variables.V[:, n_p_i, ϵ_i, ν_i, 2, age_i+1]
                    end
                    EV_func_itp = Akima(a_grid, EV)
                    function V_func_x(x::Vector, grad::Vector)
                        # x[1] = a_p, x[2] = x, x[3] = l
                        q = (μ * (x[2] / (n^ψ_1))^θ + (1.0 - μ) * (x[3] / (n^ψ_2))^θ)^(1.0 / θ)
                        return utility_function((1.0 + r) * a + (1.0 - x[3]) * w - x[1] - x[2], n, q, γ, ψ, κ) + β * EV_func_itp(x[1])
                    end
                    function c_constraint_x(x::Vector, grad::Vector)
                        return x[1] + x[2] + 10^(-16) - (1.0 + r) * a - (1.0 - x[3]) * w
                    end
                    function q_constraint_x(x::Vector, grad::Vector)
                        return q_bar - (μ * (x[2] / (n^ψ_1))^θ + (1.0 - μ) * (x[3] / (n^ψ_2))^θ)^(1.0 / θ)
                    end
                    opt = Opt(:LN_COBYLA, 3)
                    inequality_constraint!(opt, c_constraint_x)
                    inequality_constraint!(opt, q_constraint_x)
                    opt.lower_bounds = [0.0, 0.0, 0.0]
                    opt.upper_bounds = [(1.0 + r) * a + w, (1.0 + r) * a + w, 1.0]
                    opt.xtol_rel = 1e-4
                    opt.max_objective = V_func_x
                    a_p_initial = a
                    l_initial = 0.5
                    x_initial = (1.0 + r) * a + (1.0 - l_initial) * w - a_p_initial
                    opt.xtol_rel = 1e-4
                    opt.max_objective = V_func_x
                    (maxf, maxx, ret) = optimize(opt, [a_p_initial, x_initial, l_initial])
                    variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = maxf
                    variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = maxx[1]
                    variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = maxx[2]
                    variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = maxx[3]
                end
            end
        elseif age_inf <= age < age_ret # berween infertile age and retirement age
        end
    end
end

parameters = parameters_function()
variables = variables_function(parameters)
solve_value_and_policy_function!(variables, parameters)
