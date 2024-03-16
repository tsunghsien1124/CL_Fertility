#=============================#
# Solve stationary equlibrium #
#=============================#

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

function parameters_function(;
    β::Real=0.940,                # discount factor (households)
    ρ::Real=0.980,                # survival rate
    r_f::Real=0.04,               # risk-free rate # 1.04*ρ-1.0
    β_f::Real=1.0 / (1.0 + r_f),  # discount factor (bank)
    τ::Real=0.00,                 # transaction cost
    σ::Real=2.00,                 # CRRA coefficient
    δ::Real=0.08,                 # depreciation rate
    α::Real=0.36,                 # capital share
    ψ::Real=0.972^4,              # exogenous retention ratio # 1.0 - 1.0 / 20.0
    θ::Real=1.0 / (4.57 * 0.75),  # diverting fraction # 1.0 / 3.0
    p_h::Real=1.0 / 10.0,         # prob. of history erased
    η::Real=0.40,                 # wage garnishment rate
    ξ::Real=0.00,                 # stigma utility filing cost
    κ::Real=0.02,                 # out-of-pocket monetary filing cost
    e_1_σ::Real=0.448,            # s.d. of permanent endowment shock
    e_1_size::Integer=2,          # number of permanent endowment shock
    e_2_ρ::Real=0.957,            # AR(1) of persistent endowment shock
    e_2_σ::Real=0.129,            # s.d. of persistent endowment shock
    e_2_size::Integer=3,          # number of persistent endowment shock
    e_3_σ::Real=0.351,            # s.d. of transitory endowment shock
    e_3_size::Integer=3,          # number oftransitory endowment shock
    ν_size::Integer=3,            # number of expenditure shock
    a_min::Real=-5.0,             # min of asset holding
    a_max::Real=900.0,            # max of asset holding
    a_size_neg::Integer=501,      # number of grid of negative asset holding for VFI
    a_size_pos::Integer=101,      # number of grid of positive asset holding for VFI
    a_degree::Integer=2,          # curvature of the positive asset gridpoints
    μ_scale::Integer=1            # scale for the asset holding gridpoints for distribution
)
    """
    contruct an immutable object containg all paramters
    """

    # permanent endowment shock
    e_1_grid, e_1_Γ = adda_cooper(e_1_size, 0.0, e_1_σ)
    e_1_Γ = e_1_Γ[1,:]
    # e_1_grid = [-e_1_σ, e_1_σ]
    # e_1_Γ = Matrix(1.0I, e_1_size, e_1_size)
    # G_e_1 = [1.0 / e_1_size for i = 1:e_1_size]
    G_e_1 = e_1_Γ

    # persistent endowment shock
    e_2_MC = tauchen(e_2_size, e_2_ρ, e_2_σ, 0.0, 3)
    # e_2_MC = rouwenhorst(e_2_size, e_2_ρ, e_2_σ, 0.0)
    e_2_Γ = e_2_MC.p
    e_2_grid = collect(e_2_MC.state_values)
    # e_2_grid, e_2_Γ = adda_cooper(e_2_size, e_2_ρ, e_2_σ)
    G_e_2 = stationary_distributions(MarkovChain(e_2_Γ, e_2_grid))[1]
    # G_e_2 = [1.0, 0.0, 0.0]

    # transitory endowment shock
    e_3_grid, e_3_Γ = adda_cooper(e_3_size, 0.0, e_3_σ)
    e_3_Γ = e_3_Γ[1,:]
    # e_3_bar = sqrt((3 / 2) * e_3_σ^2)
    # e_3_grid = [-e_3_bar, 0.0, e_3_bar]
    # e_3_Γ = [1.0 / e_3_size for i = 1:e_3_size]
    G_e_3 = e_3_Γ # [0.0, 1.0, 0.0]

    # aggregate labor endowment
    E = 1.0

    # expenditure schock
    ν_grid = zeros(ν_size)
    ν_Γ = [0.0, 0.0, 1.0]
    G_ν = [0.0, 1.0]

    # asset holding grid for VFI
    a_grid_neg = collect(range(a_min, 0.0, length=a_size_neg))
    a_grid_pos = ((range(0.0, stop=a_size_pos - 1, length=a_size_pos) / (a_size_pos - 1)) .^ a_degree) * a_max
    a_grid = cat(a_grid_neg[1:(end-1)], a_grid_pos, dims=1)
    a_size = length(a_grid)
    a_ind_zero = findall(iszero, a_grid)[]

    # asset holding grid for μ
    a_size_neg_μ = a_size_neg * μ_scale
    a_size_pos_μ = a_size_pos * μ_scale
    a_grid_neg_μ = collect(range(a_min, 0.0, length=a_size_neg_μ))
    a_grid_pos_μ = collect(range(0.0, a_max, length = a_size_pos_μ))
    # a_grid_pos_μ = ((range(0.0, stop=a_size_pos_μ - 1, length=a_size_pos_μ) / (a_size_pos_μ - 1)) .^ a_degree) * a_max
    a_grid_μ = cat(a_grid_neg_μ[1:(end-1)], a_grid_pos_μ, dims=1)
    a_size_μ = length(a_grid_μ)
    a_ind_zero_μ = findall(iszero, a_grid_μ)[]

    # return values
    return (
        β=β,
        ρ=ρ,
        r_f=r_f,
        β_f=β_f,
        τ=τ,
        σ=σ,
        δ=δ,
        α=α,
        ψ=ψ,
        θ=θ,
        p_h=p_h,
        η=η,
        ξ=ξ,
        κ=κ,
        e_1_σ=e_1_σ,
        e_1_size=e_1_size,
        e_1_Γ=e_1_Γ,
        e_1_grid=e_1_grid,
        G_e_1=G_e_1,
        e_2_ρ=e_2_ρ,
        e_2_σ=e_2_σ,
        e_2_size=e_2_size,
        e_2_Γ=e_2_Γ,
        e_2_grid=e_2_grid,
        G_e_2=G_e_2,
        e_3_σ=e_3_σ,
        e_3_size=e_3_size,
        e_3_Γ=e_3_Γ,
        e_3_grid=e_3_grid,
        G_e_3=G_e_3,
        E=E,
        ν_size=ν_size,
        ν_Γ=ν_Γ,
        ν_grid=ν_grid,
        G_ν=G_ν,
        a_grid=a_grid,
        a_grid_neg=a_grid_neg,
        a_grid_pos=a_grid_pos,
        a_size=a_size,
        a_size_neg=a_size_neg,
        a_size_pos=a_size_pos,
        a_ind_zero=a_ind_zero,
        a_grid_μ=a_grid_μ,
        a_grid_neg_μ=a_grid_neg_μ,
        a_grid_pos_μ=a_grid_pos_μ,
        a_size_μ=a_size_μ,
        a_size_neg_μ=a_size_neg_μ,
        a_size_pos_μ=a_size_pos_μ,
        a_ind_zero_μ=a_ind_zero_μ,
        a_degree=a_degree,
    )
end

mutable struct Mutable_Aggregate_Prices
    """
    construct a type for mutable aggregate prices
    """
    λ::Real
    ξ_λ::Real
    Λ_λ::Real
    leverage_ratio_λ::Real
    KL_to_D_ratio_λ::Real
    ι_λ::Real
    r_k_λ::Real
    K_λ::Real
    w_λ::Real
end

mutable struct Mutable_Aggregate_Variables
    """
    construct a type for mutable aggregate variables
    """
    K::Real
    L::Real
    L_adj::Real
    D::Real
    N::Real
    profit::Real
    ω::Real
    leverage_ratio::Real
    KL_to_D_ratio::Real
    debt_to_earning_ratio::Real
    share_of_filers::Real
    share_of_involuntary_filers::Real
    share_in_debts::Real
    avg_loan_rate::Real
    avg_loan_rate_pw::Real
end

mutable struct Mutable_Variables
    """
    construct a type for mutable variables
    """
    aggregate_prices::Mutable_Aggregate_Prices
    aggregate_variables::Mutable_Aggregate_Variables
    R::Array{Float64,3}
    q::Array{Float64,3}
    rbl::Array{Float64,3}
    V::Array{Float64,5}
    V_d::Array{Float64,4}
    V_nd::Array{Float64,5}
    V_pos::Array{Float64,5}
    policy_a::Array{Float64,5}
    policy_d::Array{Float64,5}
    policy_pos_a::Array{Float64,5}
    threshold_a::Array{Float64,4}
    threshold_e_2::Array{Float64,4}
    μ::Array{Float64,6}
end

function min_bounds_function(obj::Function, grid_min::Real, grid_max::Real; grid_length::Integer=120, obj_range::Integer=1)
    """
    compute bounds for minimization
    """

    grid = range(grid_min, grid_max, length=grid_length)
    grid_size = length(grid)
    obj_grid = obj.(grid)
    obj_index = argmin(obj_grid)
    if obj_index < (1 + obj_range)
        lb = grid_min
        @inbounds ub = grid[obj_index+obj_range]
    elseif obj_index > (grid_size - obj_range)
        @inbounds lb = grid[obj_index-obj_range]
        ub = grid_max
    else
        @inbounds lb = grid[obj_index-obj_range]
        @inbounds ub = grid[obj_index+obj_range]
    end
    return lb, ub
end

# function min_bounds_function(obj::Function, a_grid::Array{Float64,1}, grid_min::Real, grid_max::Real; obj_range::Integer = 1)
#     """
#     compute bounds for minimization
#     """
#     grid_min_ind = findfirst(grid_min .<= a_grid)[]
#     grid_max_ind = findlast(a_grid .<= grid_max)[]
#     grid_adj = a_grid[grid_min_ind:grid_max_ind]
#     grid_size = length(grid_adj)
#     obj_grid = obj.(grid_adj)
#     obj_index = argmin(obj_grid)
#     if obj_index < (1 + obj_range)
#         lb = grid_min
#         @inbounds ub = grid_adj[obj_index+obj_range]
#     elseif obj_index > (grid_size - obj_range)
#         @inbounds lb = grid_adj[obj_index-obj_range]
#         ub = grid_max
#     else
#         @inbounds lb = grid_adj[obj_index-obj_range]
#         @inbounds ub = grid_adj[obj_index+obj_range]
#     end
#     return lb, ub
# end

function zero_bounds_function(V_d::Real, V_nd::Array{Float64,1}, a_grid::Array{Float64,1})
    """
    compute bounds for (zero) root finding
    """

    @inbounds lb = a_grid[minimum(findall(V_nd .> V_d))]
    @inbounds ub = a_grid[maximum(findall(V_nd .< V_d))]
    return lb, ub
end

function utility_function(c::Real, γ::Real)
    """
    compute utility of CRRA utility function with coefficient γ
    """

    if c > 0.0
        return γ == 1.0 ? log(c) : 1.0 / ((1.0 - γ) * c^(γ - 1.0))
    else
        return -Inf
    end
end

function log_function(threshold_e::Real)
    """
    adjusted log funciton where assigning -Inf to negative domain
    """

    if threshold_e > 0.0
        return log(threshold_e)
    else
        return -Inf
    end
end

function threshold_function(V_d::Array{Float64,4}, V_nd::Array{Float64,5}, w::Real, parameters::NamedTuple)
    """
    update default thresholds
    """

    # unpack parameters
    @unpack a_size_neg, a_grid, e_1_size, e_1_grid, e_2_size, e_2_grid, e_3_size, e_3_grid, ν_size, ν_grid = parameters

    # construct containers
    threshold_a = zeros(e_1_size, e_2_size, e_3_size, ν_size)
    threshold_e_2 = zeros(a_size_neg, e_1_size, e_3_size, ν_size)

    # loop over states
    for ν_i = 1:ν_size, e_3_i = 1:e_3_size, e_1_i = 1:e_1_size

        # println("v_i = $ν_i, e_3_i = $e_3_i, e_1_i = $e_1_i")

        # defaulting thresholds in wealth (a)
        for e_2_i = 1:e_2_size
            @inbounds @views V_nd_Non_Inf = findall(V_nd[:, e_1_i, e_2_i, e_3_i, ν_i] .!= -Inf)
            @inbounds @views a_grid_itp = a_grid[V_nd_Non_Inf]
            @inbounds @views V_nd_grid_itp = V_nd[V_nd_Non_Inf, e_1_i, e_2_i, e_3_i, ν_i]
            V_nd_itp = Akima(a_grid_itp, V_nd_grid_itp)
            @inbounds V_diff_itp(a) = V_nd_itp(a) - V_d[e_1_i, e_2_i, e_3_i, ν_i]

            if minimum(V_nd_grid_itp) > V_d[e_1_i, e_2_i, e_3_i, ν_i]
                @inbounds threshold_a[e_1_i, e_2_i, e_3_i, ν_i] = -Inf
            else
                @inbounds V_diff_lb, V_diff_ub = zero_bounds_function(V_d[e_1_i, e_2_i, e_3_i, ν_i], V_nd[:, e_1_i, e_2_i, e_3_i, ν_i], a_grid)
                @inbounds threshold_a[e_1_i, e_2_i, e_3_i, ν_i] = find_zero(a -> V_diff_itp(a), (V_diff_lb, V_diff_ub), Bisection())
            end
        end

        # defaulting thresholds in endowment (e)
        @inbounds @views thres_a_Non_Inf = findall(threshold_a[e_1_i, :, e_3_i, ν_i] .!= -Inf)
        @inbounds @views thres_a_grid_itp = -threshold_a[e_1_i, thres_a_Non_Inf, e_3_i, ν_i]
        earning_grid_itp = w * exp.(e_1_grid[e_1_i] .+ e_2_grid[thres_a_Non_Inf] .+ e_3_grid[e_3_i]) .- ν_grid[ν_i]
        threshold_earning_itp = Spline1D(thres_a_grid_itp, earning_grid_itp; k=1, bc="extrapolate")
        # threshold_earning_itp = Akima(thres_a_grid_itp, earning_grid_itp)

        Threads.@threads for a_i = 1:a_size_neg
            @inbounds earning_thres = threshold_earning_itp(-a_grid[a_i])
            e_2_thres = log_function(earning_thres / w) - e_1_grid[e_1_i] - e_3_grid[e_3_i]
            @inbounds threshold_e_2[a_i, e_1_i, e_3_i, ν_i] = e_2_thres
        end
    end

    return threshold_a, threshold_e_2
end

function repayment_function(e_1_i::Integer, e_2_i::Integer, e_3_p_i::Integer, a_p::Real, threshold_e_2::Real, w::Real, parameters::NamedTuple; wage_garnishment::Bool=true)
    """
    evaluate repayment analytically with and without wage garnishment
    """

    # unpack parameters
    @unpack e_1_grid, e_2_grid, e_3_grid, e_2_ρ, e_2_σ, η = parameters

    # permanent and transitory components
    e_1 = e_1_grid[e_1_i]
    e_3 = e_3_grid[e_3_p_i]

    # compute expected repayment amount
    @inbounds e_2_μ = e_2_ρ * e_2_grid[e_2_i]

    # (1) not default
    default_prob = cdf(Normal(e_2_μ, e_2_σ), threshold_e_2)
    amount_repay = -a_p * (1.0 - default_prob)

    # (2) default and reclaiming wage garnishment is enabled
    amount_default = 0.0
    if wage_garnishment == true
        default_adjusted_prob = cdf(Normal(e_2_μ + e_2_σ^2.0, e_2_σ), threshold_e_2)
        amount_default = η * w * exp(e_1 + e_3) * exp(e_2_μ + e_2_σ^2.0 / 2.0) * default_adjusted_prob
    end

    # (3) total amount collected by banks
    total_amount = amount_repay + amount_default

    return total_amount
end

function aggregate_prices_λ_funtion(parameters::NamedTuple; λ::Real)
    """
    compute aggregate prices for given incentive multiplier λ
    """
    @unpack ρ, α, ψ, β_f, θ, r_f, δ, E = parameters

    ξ_λ = (1.0 - ψ) / (1.0 - λ - ψ)
    Λ_λ = β_f * (1.0 - ψ + ψ * ξ_λ)
    leverage_ratio_λ = ξ_λ / θ
    KL_to_D_ratio_λ = leverage_ratio_λ / (leverage_ratio_λ - 1.0)
    ι_λ = λ * θ / Λ_λ
    r_k_λ = r_f + ι_λ
    K_λ = E * ((r_k_λ + δ) / α)^(1.0 / (α - 1.0))
    w_λ = (1.0 - α) * (K_λ / E)^α

    return ξ_λ, Λ_λ, leverage_ratio_λ, KL_to_D_ratio_λ, ι_λ, r_k_λ, K_λ, w_λ
end

function variables_function(parameters::NamedTuple; λ::Real)
    """
    construct a mutable object containing endogenous variables
    """

    # unpack parameters
    @unpack a_ind_zero, a_size, a_grid, a_size_pos, a_size_neg, a_grid_neg, a_ind_zero_μ, a_size_μ, a_size_pos_μ = parameters
    @unpack e_1_size, e_1_grid, e_1_Γ, e_2_size, e_2_grid, e_2_Γ, e_2_ρ, e_2_σ, e_3_size, e_3_grid, e_3_Γ = parameters
    @unpack ν_size, ν_Γ = parameters
    @unpack ρ, r_f, τ = parameters

    # define aggregate prices
    ξ_λ, Λ_λ, leverage_ratio_λ, KL_to_D_ratio_λ, ι_λ, r_k_λ, K_λ, w_λ = aggregate_prices_λ_funtion(parameters; λ=λ)
    aggregate_prices = Mutable_Aggregate_Prices(λ, ξ_λ, Λ_λ, leverage_ratio_λ, KL_to_D_ratio_λ, ι_λ, r_k_λ, K_λ, w_λ)

    # define aggregate variables
    K = 0.0
    L = 0.0
    L_adj = 0.0
    D = 0.0
    N = 0.0
    profit = 0.0
    ω = 0.0
    leverage_ratio = 0.0
    KL_to_D_ratio = 0.0
    debt_to_earning_ratio = 0.0
    share_of_filers = 0.0
    share_of_involuntary_filers = 0.0
    share_in_debts = 0.0
    avg_loan_rate = 0.0
    avg_loan_rate_pw = 0.0
    aggregate_variables = Mutable_Aggregate_Variables(K, L, L_adj, D, N, profit, ω, leverage_ratio, KL_to_D_ratio, debt_to_earning_ratio, share_of_filers, share_of_involuntary_filers, share_in_debts, avg_loan_rate, avg_loan_rate_pw)

    # define repayment probability, pricing function, and risky borrowing limit
    R = zeros(a_size_neg, e_1_size, e_2_size)
    q = ones(a_size, e_1_size, e_2_size) .* ρ ./ (1.0 + r_f)
    rbl = zeros(e_1_size, e_2_size, 2)
    for e_2_i = 1:e_2_size, e_1_i = 1:e_1_size
        for a_p_i = 1:(a_size_neg-1)
            @inbounds a_p = a_grid_neg[a_p_i]
            for ν_p_i = 1:ν_size, e_3_p_i = 1:e_3_size
                @inbounds threshold_e_2 = log_function(-a_p / w_λ) - e_1_grid[e_1_i] - e_3_grid[e_3_p_i]
                @inbounds R[a_p_i, e_1_i, e_2_i] += e_3_Γ[e_3_p_i] * ν_Γ[ν_p_i] * repayment_function(e_1_i, e_2_i, e_3_p_i, a_p, threshold_e_2, w_λ, parameters)
            end
            @inbounds q[a_p_i, e_1_i, e_2_i] = ρ * R[a_p_i, e_1_i, e_2_i] / ((-a_p) * (1.0 + r_f + τ + ι_λ))
        end

        qa_funcion_itp = Akima(a_grid_neg, q[1:a_ind_zero, e_1_i, e_2_i] .* a_grid_neg)
        qa_funcion(a_p) = qa_funcion_itp(a_p)
        @inbounds rbl_lb, rbl_ub = min_bounds_function(qa_funcion, a_grid[1], 0.0)
        # @inbounds rbl_lb, rbl_ub = min_bounds_function(qa_funcion, a_grid_neg, a_grid[1], 0.0)
        res_rbl = optimize(qa_funcion, rbl_lb, rbl_ub)
        @inbounds rbl[e_1_i, e_2_i, 1] = Optim.minimizer(res_rbl)
        @inbounds rbl[e_1_i, e_2_i, 2] = Optim.minimum(res_rbl)
    end

    # define value and policy functions
    V = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size)
    V_d = zeros(e_1_size, e_2_size, e_3_size, ν_size)
    V_nd = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size)
    V_pos = zeros(a_size_pos, e_1_size, e_2_size, e_3_size, ν_size)
    policy_a = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size)
    policy_d = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size)
    policy_pos_a = zeros(a_size_pos, e_1_size, e_2_size, e_3_size, ν_size)
    threshold_a = zeros(e_1_size, e_2_size, e_3_size, ν_size)
    threshold_e_2 = zeros(a_size_neg, e_1_size, e_3_size, ν_size)

    # define cross-sectional distribution
    μ = zeros(a_size_μ, e_1_size, e_2_size, e_3_size, ν_size, 2)
    μ_size = (a_size_μ + a_size_pos_μ) * e_1_size * e_2_size * e_3_size * ν_size
    μ[:, :, :, :, :, 1] .= 1.0 ./ μ_size
    μ[a_ind_zero_μ:end, :, :, :, :, 2] .= 1.0 ./ μ_size

    # return outputs
    variables = Mutable_Variables(aggregate_prices, aggregate_variables, R, q, rbl, V, V_d, V_nd, V_pos, policy_a, policy_d, policy_pos_a, threshold_a, threshold_e_2, μ)
    return variables
end

function variables_function_update!(variables::Mutable_Variables, parameters::NamedTuple; λ::Real)
    """
    construct a mutable object containing endogenous variables
    """

    # define aggregate prices
    ξ_λ, Λ_λ, leverage_ratio_λ, KL_to_D_ratio_λ, ι_λ, r_k_λ, K_λ, w_λ = aggregate_prices_λ_funtion(parameters; λ=λ)
    variables.aggregate_prices = Mutable_Aggregate_Prices(λ, ξ_λ, Λ_λ, leverage_ratio_λ, KL_to_D_ratio_λ, ι_λ, r_k_λ, K_λ, w_λ)
end

function EV_function(e_1_i::Integer, e_2_i::Integer, V_p::Array{Float64,5}, parameters::NamedTuple)
    """
    construct expected value function
    """

    # unpack parameters
    @unpack e_2_size, e_2_Γ, e_3_size, e_3_Γ, ν_size, ν_Γ = parameters

    # construct container
    a_size_ = size(V_p)[1]
    EV = zeros(a_size_)

    # update expected value
    for ν_p_i = 1:ν_size, e_3_p_i = 1:e_3_size, e_2_p_i = 1:e_2_size
        @inbounds @views EV += e_2_Γ[e_2_i, e_2_p_i] * e_3_Γ[e_3_p_i] * ν_Γ[ν_p_i] * V_p[:, e_1_i, e_2_p_i, e_3_p_i, ν_p_i]
    end

    # repalce NaN with -Inf
    replace!(EV, NaN => -Inf)

    # return value
    return EV
end

function value_and_policy_function(
    V_p::Array{Float64,5},
    V_d_p::Array{Float64,4},
    V_nd_p::Array{Float64,5},
    V_pos_p::Array{Float64,5},
    q::Array{Float64,3},
    rbl::Array{Float64,3},
    w::Real,
    parameters::NamedTuple;
    slow_updating::Real=1.0
)
    """
    one-step update of value and policy functions
    """

    # unpack parameters
    @unpack a_size, a_grid, a_size_pos, a_grid_pos, a_ind_zero = parameters
    @unpack e_1_size, e_1_grid, e_1_Γ, e_2_size, e_2_grid, e_2_Γ, e_3_size, e_3_grid, e_3_Γ = parameters
    @unpack ν_size, ν_grid, ν_Γ = parameters
    @unpack ρ, β, σ, r_f = parameters
    @unpack p_h, η, κ, ξ = parameters

    # construct containers
    V = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size)
    V_d = zeros(e_1_size, e_2_size, e_3_size, ν_size)
    V_nd = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size)
    V_pos = zeros(a_size_pos, e_1_size, e_2_size, e_3_size, ν_size)
    policy_a = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size)
    policy_d = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size)
    policy_pos_a = zeros(a_size_pos, e_1_size, e_2_size, e_3_size, ν_size)

    # loop over all states
    for ν_i = 1:ν_size, e_3_i = 1:e_3_size, e_2_i = 1:e_2_size, e_1_i = 1:e_1_size

        # extract unexpected expenses
        @inbounds ν = ν_grid[ν_i]

        # construct earning
        @inbounds y = w * exp(e_1_grid[e_1_i] + e_2_grid[e_2_i] + e_3_grid[e_3_i])

        # extract risky borrowing limit and maximum discounted borrowing amount
        @inbounds @views rbl_a, rbl_qa = rbl[e_1_i, e_2_i, :]

        # construct interpolated discounted borrowing amount functions
        @inbounds @views qa = q[:, e_1_i, e_2_i] .* a_grid
        qa_function_itp = Akima(a_grid, qa)

        # compute the next-period discounted expected value funtions and interpolated functions
        V_hat = ρ * β * EV_function(e_1_i, e_2_i, V_p, parameters)
        V_hat_pos = ρ * β * EV_function(e_1_i, e_2_i, V_pos_p, parameters)
        V_hat_itp = Akima(a_grid, V_hat)
        V_hat_pos_itp = Akima(a_grid_pos, (p_h * V_hat[a_ind_zero:end] + (1.0 - p_h) * V_hat_pos))

        # compute defaulting value
        @inbounds V_d[e_1_i, e_2_i, e_3_i, ν_i] = utility_function((1.0 - η) * y - κ, σ) - ξ + V_hat_pos[1]

        # compute non-defaulting value
        Threads.@threads for a_i = 1:a_size

            # cash on hand
            @inbounds a = a_grid[a_i]
            CoH = y + a - ν

            # compute non-defaulting value
            if (CoH - rbl_qa) >= 0.0
                object_nd(a_p) = -(utility_function(CoH - qa_function_itp(a_p), σ) + V_hat_itp(a_p))
                lb, ub = min_bounds_function(object_nd, rbl_a - eps(), CoH + eps())
                res_nd = optimize(object_nd, lb, ub)
                @inbounds V_nd[a_i, e_1_i, e_2_i, e_3_i, ν_i] = -Optim.minimum(res_nd)
                @inbounds policy_a[a_i, e_1_i, e_2_i, e_3_i, ν_i] = Optim.minimizer(res_nd)
            else
                # involuntary default
                @inbounds V_nd[a_i, e_1_i, e_2_i, e_3_i, ν_i] = -Inf
            end

            # compute value with good credit history
            V_max = max(V_nd[a_i, e_1_i, e_2_i, e_3_i, ν_i], V_d[e_1_i, e_2_i, e_3_i, ν_i])
            if V_max == V_nd[a_i, e_1_i, e_2_i, e_3_i, ν_i]
                V[a_i, e_1_i, e_2_i, e_3_i, ν_i] = V_nd[a_i, e_1_i, e_2_i, e_3_i, ν_i]
                policy_d[a_i, e_1_i, e_2_i, e_3_i, ν_i] = 0.0
            else
                V[a_i, e_1_i, e_2_i, e_3_i, ν_i] = V_d[e_1_i, e_2_i, e_3_i, ν_i]
                policy_d[a_i, e_1_i, e_2_i, e_3_i, ν_i] = 1.0
            end

            # bad credit history
            if a_i >= a_ind_zero
                a_pos_i = a_i - a_ind_zero + 1
                object_pos(a_p) = -(utility_function(CoH - qa_function_itp(a_p), σ) + V_hat_pos_itp(a_p))
                lb, ub = min_bounds_function(object_pos, 0.0, CoH + eps())
                res_pos = optimize(object_pos, lb, ub)
                @inbounds V_pos[a_pos_i, e_1_i, e_2_i, e_3_i, ν_i] = -Optim.minimum(res_pos)
                @inbounds policy_pos_a[a_pos_i, e_1_i, e_2_i, e_3_i, ν_i] = Optim.minimizer(res_pos)
            end
        end
    end

    # slow updating
    if slow_updating != 1.0
        V = slow_updating * V + (1.0 - slow_updating) * V_p
        V_d = slow_updating * V_d + (1.0 - slow_updating) * V_d_p
        V_nd = slow_updating * V_nd + (1.0 - slow_updating) * V_nd_p
        V_pos = slow_updating * V_pos + (1.0 - slow_updating) * V_pos_p
    end

    # return results
    return V, V_d, V_nd, V_pos, policy_a, policy_d, policy_pos_a
end

function pricing_and_rbl_function(threshold_e_2::Array{Float64,4}, w::Real, ι::Real, parameters::NamedTuple)
    """
    update pricing function and borrowing risky limit
    """

    # unpack parameters
    @unpack ρ, r_f, τ, η = parameters
    @unpack a_ind_zero, a_size, a_grid, a_size_neg, a_grid_neg = parameters
    @unpack e_1_size, e_1_grid, e_1_Γ, e_2_size, e_2_grid, e_2_Γ, e_3_size, e_3_grid, e_3_Γ = parameters
    @unpack ν_size, ν_Γ = parameters

    # contruct containers
    R = zeros(a_size_neg, e_1_size, e_2_size)
    q = ones(a_size, e_1_size, e_2_size) .* ρ ./ (1.0 + r_f)
    rbl = zeros(e_1_size, e_2_size, 2)

    # loop over states
    for e_2_i = 1:e_2_size, e_1_i = 1:e_1_size
        for a_p_i = 1:(a_size_neg-1)
            @inbounds a_p = a_grid[a_p_i]
            for ν_p_i = 1:ν_size, e_3_p_i = 1:e_3_size
                @inbounds R[a_p_i, e_1_i, e_2_i] += e_3_Γ[e_3_p_i] * ν_Γ[ν_p_i] * repayment_function(e_1_i, e_2_i, e_3_p_i, a_p, threshold_e_2[a_p_i, e_1_i, e_3_p_i, ν_p_i], w, parameters)
            end
            @inbounds q[a_p_i, e_1_i, e_2_i] = ρ * R[a_p_i, e_1_i, e_2_i] / ((-a_p) * (1.0 + r_f + τ + ι))
        end

        # risky borrowing limit and maximum discounted borrwoing amount
        qa_funcion_itp = Akima(a_grid_neg, q[1:a_ind_zero, e_1_i, e_2_i] .* a_grid_neg)
        # qa_funcion_itp = Spline1D(a_grid_neg, q[1:a_ind_zero, e_1_i, e_2_i] .* a_grid_neg; k = 1, bc = "extrapolate")
        qa_funcion(a_p) = qa_funcion_itp(a_p)
        @inbounds rbl_lb, rbl_ub = min_bounds_function(qa_funcion, a_grid[1], 0.0)
        # @inbounds rbl_lb, rbl_ub = min_bounds_function(qa_funcion, a_grid_neg, a_grid[1], 0.0)
        res_rbl = optimize(qa_funcion, rbl_lb, rbl_ub)
        @inbounds rbl[e_1_i, e_2_i, 1] = Optim.minimizer(res_rbl)
        @inbounds rbl[e_1_i, e_2_i, 2] = Optim.minimum(res_rbl)
    end

    # return results
    return R, q, rbl
end

function solve_value_and_pricing_function!(variables::Mutable_Variables, parameters::NamedTuple; tol::Real=1E-8, iter_max::Integer=1000, slow_updating::Real=1.0)
    """
    solve household and banking problems using one-loop algorithm
    """

    # initialize the iteration number and criterion
    search_iter = 0
    crit = Inf
    prog = ProgressThresh(tol, "Solving household and banking problems (one-loop): ")

    # construct containers
    V_p = similar(variables.V)
    V_d_p = similar(variables.V_d)
    V_nd_p = similar(variables.V_nd)
    V_pos_p = similar(variables.V_pos)
    q_p = similar(variables.q)

    while crit > tol && search_iter < iter_max

        # copy previous values
        copyto!(V_p, variables.V)
        copyto!(V_d_p, variables.V_d)
        copyto!(V_nd_p, variables.V_nd)
        copyto!(V_pos_p, variables.V_pos)
        copyto!(q_p, variables.q)

        # value and policy functions
        variables.V, variables.V_d, variables.V_nd, variables.V_pos, variables.policy_a, variables.policy_d, variables.policy_pos_a =
            value_and_policy_function(V_p, V_d_p, V_nd_p, V_pos_p, variables.q, variables.rbl, variables.aggregate_prices.w_λ, parameters; slow_updating=slow_updating)

        # default thresholds
        variables.threshold_a, variables.threshold_e_2 = threshold_function(variables.V_d, variables.V_nd, variables.aggregate_prices.w_λ, parameters)

        # pricing function and borrowing risky limit
        variables.R, variables.q, variables.rbl = pricing_and_rbl_function(variables.threshold_e_2, variables.aggregate_prices.w_λ, variables.aggregate_prices.ι_λ, parameters)

        # check convergence
        V_crit = norm(variables.V .- V_p, Inf)
        V_pos_crit = norm(variables.V_pos .- V_pos_p, Inf)
        q_crit = norm(variables.q .- q_p, Inf)
        crit = max(V_crit, V_pos_crit, q_crit)

        # update the iteration number
        search_iter += 1

        # manually report convergence progress
        # println("|V| = $V_crit, |V_pos| = $V_pos_crit, |q| = $q_crit")
        # println("Solving household and banking problems (one-loop): search_iter = $search_iter and crit = $crit > tol = $tol")
        ProgressMeter.update!(prog, crit)
    end

    return crit
end

function stationary_distribution_function(μ_p::Array{Float64,6}, policy_a::Array{Float64,5}, threshold_a::Array{Float64,4}, policy_pos_a::Array{Float64,5}, parameters::NamedTuple)
    """
    update stationary distribution
    """

    # unpack parameters
    @unpack e_1_size, e_1_Γ, G_e_1, e_2_size, e_2_Γ, G_e_2, e_3_size, e_3_Γ, G_e_3, ν_size, ν_Γ, G_ν, a_grid, a_grid_pos, a_size_μ, a_grid_μ, a_ind_zero_μ, ρ, p_h = parameters

    # construct container
    μ = zeros(a_size_μ, e_1_size, e_2_size, e_3_size, ν_size, 2)

    for e_1_i = 1:e_1_size, e_2_i = 1:e_2_size, e_3_i = 1:e_3_size, ν_i = 1:ν_size

        # interpolated decision rules
        @inbounds @views policy_a_Non_Inf = findall(policy_a[:, e_1_i, e_2_i, e_3_i, ν_i] .!= -Inf)
        @inbounds policy_a_itp = Akima(a_grid[policy_a_Non_Inf], policy_a[policy_a_Non_Inf, e_1_i, e_2_i, e_3_i, ν_i])
        # @inbounds policy_d_itp = Akima(a_grid, policy_d[:, e_1_i, e_2_i, e_3_i, ν_i])
        @inbounds policy_d_itp(x) = x < threshold_a[e_1_i, e_2_i, e_3_i, ν_i] ? 1.0 : 0.0
        @inbounds policy_pos_a_itp = Akima(a_grid_pos, policy_pos_a[:, e_1_i, e_2_i, e_3_i, ν_i])

        # loop over the dimension of asset holding
        for a_μ_i = 1:a_size_μ

            # extract wealth and compute asset choice
            @inbounds a_μ = a_grid_μ[a_μ_i]
            @inbounds a_p = clamp(policy_a_itp(a_μ), a_grid[1], a_grid[end])

            # locate it on the original grid
            a_p_lb = findall(a_grid_μ .<= a_p)[end]
            a_p_ub = findall(a_p .<= a_grid_μ)[1]

            # compute weights
            if a_p_lb != a_p_ub
                @inbounds a_p_lower = a_grid_μ[a_p_lb]
                @inbounds a_p_upper = a_grid_μ[a_p_ub]
                weight_lower = (a_p_upper - a_p) / (a_p_upper - a_p_lower)
                weight_upper = (a_p - a_p_lower) / (a_p_upper - a_p_lower)
            else
                weight_lower = 0.5
                weight_upper = 0.5
            end

            # loop over the dimension of exogenous individual states
            for e_1_p_i = 1:e_1_size, e_2_p_i = 1:e_2_size, e_3_p_i = 1:e_3_size, ν_p_i = 1:ν_size
                if e_1_p_i == e_1_i
                    if policy_d_itp(a_μ) == 0.0
                        @inbounds μ[a_p_lb, e_1_i, e_2_p_i, e_3_p_i, ν_p_i, 1] += ρ * e_2_Γ[e_2_i, e_2_p_i] * e_3_Γ[e_3_p_i] * ν_Γ[ν_p_i] * weight_lower * μ_p[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1]
                        @inbounds μ[a_p_ub, e_1_i, e_2_p_i, e_3_p_i, ν_p_i, 1] += ρ * e_2_Γ[e_2_i, e_2_p_i] * e_3_Γ[e_3_p_i] * ν_Γ[ν_p_i] * weight_upper * μ_p[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1]
                    else
                        @inbounds μ[a_ind_zero_μ, e_1_i, e_2_p_i, e_3_p_i, ν_p_i, 2] += ρ * e_2_Γ[e_2_i, e_2_p_i] * e_3_Γ[e_3_p_i] * ν_Γ[ν_p_i] * μ_p[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1]
                    end
                end
                @inbounds μ[a_ind_zero_μ, e_1_p_i, e_2_p_i, e_3_p_i, ν_p_i, 1] += (1.0 - ρ) * G_e_1[e_1_p_i] * G_e_2[e_2_p_i] * G_e_3[e_3_p_i] * G_ν[ν_p_i] * μ_p[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1]
            end

            if a_μ >= 0.0
                @inbounds a_p = clamp(policy_pos_a_itp(a_μ), 0.0, a_grid[end])
                a_p_lb = findall(a_grid_μ .<= a_p)[end]
                a_p_ub = findall(a_p .<= a_grid_μ)[1]
                if a_p_lb != a_p_ub
                    @inbounds a_p_lower = a_grid_μ[a_p_lb]
                    @inbounds a_p_upper = a_grid_μ[a_p_ub]
                    weight_lower = (a_p_upper - a_p) / (a_p_upper - a_p_lower)
                    weight_upper = (a_p - a_p_lower) / (a_p_upper - a_p_lower)
                else
                    weight_lower = 0.5
                    weight_upper = 0.5
                end
                for e_1_p_i = 1:e_1_size, e_2_p_i = 1:e_2_size, e_3_p_i = 1:e_3_size, ν_p_i = 1:ν_size
                    if e_1_p_i == e_1_i
                        @inbounds μ[a_p_lb, e_1_i, e_2_p_i, e_3_p_i, ν_p_i, 1] += ρ * p_h * e_2_Γ[e_2_i, e_2_p_i] * e_3_Γ[e_3_p_i] * ν_Γ[ν_p_i] * weight_lower * μ_p[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 2]
                        @inbounds μ[a_p_ub, e_1_i, e_2_p_i, e_3_p_i, ν_p_i, 1] += ρ * p_h * e_2_Γ[e_2_i, e_2_p_i] * e_3_Γ[e_3_p_i] * ν_Γ[ν_p_i] * weight_upper * μ_p[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 2]
                        @inbounds μ[a_p_lb, e_1_i, e_2_p_i, e_3_p_i, ν_p_i, 2] += ρ * (1.0 - p_h) * e_2_Γ[e_2_i, e_2_p_i] * e_3_Γ[e_3_p_i] * ν_Γ[ν_p_i] * weight_lower * μ_p[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 2]
                        @inbounds μ[a_p_ub, e_1_i, e_2_p_i, e_3_p_i, ν_p_i, 2] += ρ * (1.0 - p_h) * e_2_Γ[e_2_i, e_2_p_i] * e_3_Γ[e_3_p_i] * ν_Γ[ν_p_i] * weight_upper * μ_p[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 2]
                    end
                    @inbounds μ[a_ind_zero_μ, e_1_p_i, e_2_p_i, e_3_p_i, ν_p_i, 1] += (1.0 - ρ) * G_e_1[e_1_p_i] * G_e_2[e_2_p_i] * G_e_3[e_3_p_i] * G_ν[ν_p_i] * μ_p[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 2]
                end
            end
        end
    end

    # standardize distribution
    # sum_μ = sum(μ)
    # println("sum_μ = $sum_μ")
    μ = μ ./ sum(μ)

    # return result
    return μ
end

function solve_stationary_distribution_function!(variables::Mutable_Variables, parameters::NamedTuple; tol::Real=1E-8, iter_max::Integer=2000)
    """
    solve stationary distribution
    """

    # initialize the iteration number and criterion
    search_iter = 0
    crit = Inf
    prog = ProgressThresh(tol, "Solving stationary distribution: ")

    # construct container
    μ_p = similar(variables.μ)

    while crit > tol && search_iter < iter_max

        # copy previous value
        copyto!(μ_p, variables.μ)

        # update stationary distribution
        variables.μ = stationary_distribution_function(μ_p, variables.policy_a, variables.threshold_a, variables.policy_pos_a, parameters)

        # check convergence
        crit = norm(variables.μ .- μ_p, Inf)

        # update the iteration number
        search_iter += 1

        # manually report convergence progress
        # println("Solving stationary distribution: search_iter = $search_iter and crit = $crit > tol = $tol")
        ProgressMeter.update!(prog, crit)
    end

    return crit
end

function solve_aggregate_variable_function(
    policy_a::Array{Float64,5},
    threshold_a::Array{Float64,4},
    policy_pos_a::Array{Float64,5},
    q::Array{Float64,3},
    rbl::Array{Float64,3},
    μ::Array{Float64,6},
    K::Real,
    w::Real,
    ι::Real,
    parameters::NamedTuple,
)
    """
    compute equlibrium aggregate variables
    """

    # unpack parameters
    @unpack e_1_size, e_1_grid, e_2_size, e_2_grid, e_3_size, e_3_grid, ν_size, a_grid, a_grid_neg, a_grid_pos, a_ind_zero_μ, a_grid_pos_μ, a_grid_neg_μ, a_size_neg_μ, a_grid_μ, a_size_μ, r_f, τ, ψ, η = parameters

    # initialize container
    K = K
    L = 0.0
    L_adj = 0.0
    D = 0.0
    N = 0.0
    profit = 0.0
    ω = 0.0
    leverage_ratio = 0.0
    KL_to_D_ratio = 0.0
    debt_to_earning_ratio = 0.0
    debt_to_earning_ratio_num = 0.0
    debt_to_earning_ratio_den = 0.0
    share_of_filers = 0.0
    share_of_involuntary_filers = 0.0
    share_in_debts = 0.0
    avg_loan_rate = 0.0
    avg_loan_rate_num = 0.0
    avg_loan_rate_den = 0.0
    avg_loan_rate_pw = 0.0
    avg_loan_rate_pw_num = 0.0
    avg_loan_rate_pw_den = 0.0

    # total loans, deposits, share of filers, nad debt-to-earning ratio
    for e_1_i = 1:e_1_size, e_2_i = 1:e_2_size, e_3_i = 1:e_3_size, ν_i = 1:ν_size

        # interpolated decision rules
        @inbounds @views policy_a_Non_Inf = findall(policy_a[:, e_1_i, e_2_i, e_3_i, ν_i] .!= -Inf)
        @inbounds policy_a_itp = Akima(a_grid[policy_a_Non_Inf], policy_a[policy_a_Non_Inf, e_1_i, e_2_i, e_3_i, ν_i])
        # @inbounds policy_d_itp = Akima(a_grid, policy_d[:, e_1_i, e_2_i, e_3_i, ν_i])
        @inbounds policy_d_itp(x) = x < threshold_a[e_1_i, e_2_i, e_3_i, ν_i] ? 1.0 : 0.0
        @inbounds policy_pos_a_itp = Akima(a_grid_pos, policy_pos_a[:, e_1_i, e_2_i, e_3_i, ν_i])

        # interpolated discounted borrowing amount
        @inbounds @views q_e = q[:, e_1_i, e_2_i]
        q_function_itp = Akima(a_grid, q_e)
        qa_function_itp = Akima(a_grid, q_e .* a_grid)

        # loop over the dimension of asset holding
        for a_μ_i = 1:a_size_μ

            # extract wealth and compute asset choice
            @inbounds a_μ = a_grid_μ[a_μ_i]
            @inbounds a_p = clamp(policy_a_itp(a_μ), a_grid[1], a_grid[end])

            if a_p < 0.0
                # total loans
                @inbounds L += -(μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (1.0 - policy_d_itp(a_μ)) * qa_function_itp(a_p))

                # average loan rate
                avg_loan_rate_num += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (1.0 - policy_d_itp(a_μ)) * (1.0 / q_function_itp(a_p) - 1.0)
                avg_loan_rate_den += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (1.0 - policy_d_itp(a_μ))

                # average loan rate (persons-weighted)
                avg_loan_rate_pw_num += (1.0 - policy_d_itp(a_μ)) * (1.0 / q_function_itp(a_p) - 1.0)
                avg_loan_rate_pw_den += 1
            else
                # total deposits
                if a_p > 0.0
                    @inbounds D += (μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (1.0 - policy_d_itp(a_μ)) * qa_function_itp(a_p))
                    # @inbounds D += (μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * qa_function_itp(a_p))

                end
            end

            if a_μ >= 0.0
                @inbounds a_pos_p = clamp(policy_pos_a_itp(a_μ), 0.0, a_grid[end])
                if a_pos_p > 0.0
                    @inbounds D += (μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 2] * qa_function_itp(a_pos_p))
                end
            end

            if a_μ < 0.0
                # share of filers
                @inbounds share_of_filers += (μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * policy_d_itp(a_μ))

                # share of involuntary filers
                if w * exp(e_1_grid[e_1_i] + e_2_grid[e_2_i] + e_3_grid[e_3_i]) + a_μ - rbl[e_1_i, e_2_i, 2] < 0.0
                    @inbounds share_of_involuntary_filers += (μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * policy_d_itp(a_μ))
                end

                # debt-to-earning ratio
                # @inbounds debt_to_earning_ratio += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (-a_μ / (w * exp(e_1_grid[e_1_i] + e_2_grid[e_2_i] + e_3_grid[e_3_i])))
                @inbounds debt_to_earning_ratio_num += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (-a_μ)
                # @inbounds debt_to_earning_ratio_den += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (w * exp(e_1_grid[e_1_i] + e_2_grid[e_2_i] + e_3_grid[e_3_i]))

                # loans returned
                L_adj += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * ((-a_μ) * (1.0 - policy_d_itp(a_μ)) + policy_d_itp(a_μ) * η * w * exp(e_1_grid[e_1_i] + e_2_grid[e_2_i] + e_3_grid[e_3_i]))
            end

            @inbounds debt_to_earning_ratio_den += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (w * exp(e_1_grid[e_1_i] + e_2_grid[e_2_i] + e_3_grid[e_3_i]))
            @inbounds debt_to_earning_ratio_den += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 2] * (w * exp(e_1_grid[e_1_i] + e_2_grid[e_2_i] + e_3_grid[e_3_i]))
        end
    end

    # net worth
    N = (K + L) - D

    # exogenous dividend policy
    # profit = (1.0 + r_f + ι) * K + (1.0 + τ + ι) * L - (1.0 + r_f) * D
    profit = ι * (K + L) + (1.0 + r_f) * N
    # ω = (N - ψ * profit) / ((1.0 - ψ) * profit)
    # ω = N / (ψ * profit)
    # ω = (N - ψ * profit) / ((1.0 - ψ) * (K + L))
    ω = (N - ψ * profit) / (K + L)
    # ω = N - ψ * profit

    # leverage ratio
    leverage_ratio = (K + L) / N

    # capital-loan-to-deposit ratio
    KL_to_D_ratio = (K + L) / D

    # debt-to-earning ratio
    # debt_to_earning_ratio = debt_to_earning_ratio_num / debt_to_earning_ratio_den
    # debt_to_earning_ratio = L / w
    debt_to_earning_ratio = debt_to_earning_ratio_num / w

    # average loan rate
    avg_loan_rate = avg_loan_rate_num / avg_loan_rate_den
    avg_loan_rate_pw = avg_loan_rate_pw_num / avg_loan_rate_pw_den

    # share in debt
    share_in_debts = sum(μ[1:(a_ind_zero_μ-1), :, :, :, :, 1])

    # return results
    aggregate_variables = Mutable_Aggregate_Variables(K, L, L_adj, D, N, profit, ω, leverage_ratio, KL_to_D_ratio, debt_to_earning_ratio, share_of_filers, share_of_involuntary_filers, share_in_debts, avg_loan_rate, avg_loan_rate_pw)
    return aggregate_variables
end

function solve_aggregate_variable_across_HH_function(
    policy_a::Array{Float64,5},
    policy_d::Array{Float64,5},
    policy_pos_a::Array{Float64,5},
    q::Array{Float64,3},
    μ::Array{Float64,6},
    w::Real,
    parameters::NamedTuple,
)
    """
    compute equlibrium aggregate variables
    """

    # unpack parameters
    @unpack e_1_size, e_1_grid, e_2_size, e_2_grid, e_3_size, e_3_grid, ν_size, a_grid, a_grid_neg, a_grid_pos, a_ind_zero_μ, a_grid_pos_μ, a_grid_neg_μ, a_size_neg_μ, a_grid_μ, a_size_μ, r_f, τ, ψ, η = parameters

    # initialize container
    debt_to_earning_ratio = 0.0
    debt_to_earning_ratio_permanent_low = 0.0
    debt_to_earning_ratio_permanent_high = 0.0

    debt_to_earning_ratio_num = 0.0
    debt_to_earning_ratio_num_permanent_low = 0.0
    debt_to_earning_ratio_num_permanent_high = 0.0

    share_of_filers = 0.0
    share_of_filers_permanent_low = 0.0
    share_of_filers_permanent_high = 0.0

    share_in_debts = 0.0
    share_in_debts_permanent_low = 0.0
    share_in_debts_permanent_high = 0.0

    avg_loan_rate = 0.0
    avg_loan_rate_num = 0.0
    avg_loan_rate_den = 0.0

    avg_loan_rate_permanent_low = 0.0
    avg_loan_rate_num_permanent_low = 0.0
    avg_loan_rate_den_permanent_low = 0.0

    avg_loan_rate_permanent_high = 0.0
    avg_loan_rate_num_permanent_high = 0.0
    avg_loan_rate_den_permanent_high = 0.0

    # total loans, deposits, share of filers, nad debt-to-earning ratio
    for e_1_i = 1:e_1_size, e_2_i = 1:e_2_size, e_3_i = 1:e_3_size, ν_i = 1:ν_size

        # interpolated decision rules
        @inbounds @views policy_a_Non_Inf = findall(policy_a[:, e_1_i, e_2_i, e_3_i, ν_i] .!= -Inf)
        @inbounds policy_a_itp = Akima(a_grid[policy_a_Non_Inf], policy_a[policy_a_Non_Inf, e_1_i, e_2_i, e_3_i, ν_i])
        @inbounds policy_d_itp = Akima(a_grid, policy_d[:, e_1_i, e_2_i, e_3_i, ν_i])
        @inbounds policy_pos_a_itp = Akima(a_grid_pos, policy_pos_a[:, e_1_i, e_2_i, e_3_i, ν_i])

        # interpolated discounted borrowing amount
        @inbounds @views q_e = q[:, e_1_i, e_2_i]
        q_function_itp = Akima(a_grid, q_e)

        # loop over the dimension of asset holding
        for a_μ_i = 1:a_size_μ

            # extract wealth and compute asset choice
            @inbounds a_μ = a_grid_μ[a_μ_i]
            @inbounds a_p = clamp(policy_a_itp(a_μ), a_grid[1], a_grid[end])

            if a_p < 0.0
                # average loan rate
                avg_loan_rate_num += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (1.0 - policy_d_itp(a_μ)) * (1.0 / q_function_itp(a_p) - 1.0)
                avg_loan_rate_den += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (1.0 - policy_d_itp(a_μ))
                if (e_1_i == 1) && (e_2_i == 2)
                    avg_loan_rate_num_permanent_low += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (1.0 - policy_d_itp(a_μ)) * (1.0 / q_function_itp(a_p) - 1.0)
                    avg_loan_rate_den_permanent_low += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (1.0 - policy_d_itp(a_μ))
                end
                if (e_1_i == 2) && (e_2_i == 2)
                    avg_loan_rate_num_permanent_high += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (1.0 - policy_d_itp(a_μ)) * (1.0 / q_function_itp(a_p) - 1.0)
                    avg_loan_rate_den_permanent_high += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (1.0 - policy_d_itp(a_μ))
                end
            end

            if a_μ < 0.0
                # share of filers
                @inbounds share_of_filers += (μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * policy_d_itp(a_μ))
                if (e_1_i == 1) && (e_2_i == 2)
                    share_of_filers_permanent_low += (μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * policy_d_itp(a_μ)) / sum(μ[:, e_1_i, e_2_i, :, :, :])
                end
                if (e_1_i == 2) && (e_2_i == 2)
                    share_of_filers_permanent_high += (μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * policy_d_itp(a_μ)) / sum(μ[:, e_1_i, e_2_i, :, :, :])
                end

                # debt-to-earning ratio
                @inbounds debt_to_earning_ratio_num += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (-a_μ)
                if (e_1_i == 1) && (e_2_i == 2)
                    @inbounds debt_to_earning_ratio_num_permanent_low += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (-a_μ) / sum(μ[:, e_1_i, e_2_i, :, :, :])
                end
                if (e_1_i == 2) && (e_2_i == 2)
                    @inbounds debt_to_earning_ratio_num_permanent_high += μ[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1] * (-a_μ) / sum(μ[:, e_1_i, e_2_i, :, :, :])
                end
            end
        end
    end

    # debt-to-earning ratio
    debt_to_earning_ratio = debt_to_earning_ratio_num / w
    debt_to_earning_ratio_permanent_low = debt_to_earning_ratio_num_permanent_low / (w * exp(e_1_grid[1]))
    debt_to_earning_ratio_permanent_high = debt_to_earning_ratio_num_permanent_high / (w * exp(e_1_grid[2]))

    # average loan rate
    avg_loan_rate = avg_loan_rate_num / avg_loan_rate_den
    avg_loan_rate_permanent_low = avg_loan_rate_num_permanent_low / avg_loan_rate_den_permanent_low
    avg_loan_rate_permanent_high = avg_loan_rate_num_permanent_high / avg_loan_rate_den_permanent_high

    # share in debt
    share_in_debts = sum(μ[1:(a_ind_zero_μ-1), :, :, :, :, 1])
    share_in_debts_permanent_low = sum(μ[1:(a_ind_zero_μ-1), 1, 2, :, :, 1]) ./ sum(μ[:, 1, 2, :, :, :])
    share_in_debts_permanent_high = sum(μ[1:(a_ind_zero_μ-1), 2, 2, :, :, 1]) ./ sum(μ[:, 2, 2, :, :, :])

    # return results
    return debt_to_earning_ratio, debt_to_earning_ratio_permanent_low, debt_to_earning_ratio_permanent_high, share_of_filers, share_of_filers_permanent_low, share_of_filers_permanent_high, share_in_debts, share_in_debts_permanent_low, share_in_debts_permanent_high, avg_loan_rate, avg_loan_rate_permanent_low, avg_loan_rate_permanent_high
end

function solve_economy_function!(variables::Mutable_Variables, parameters::NamedTuple; tol_h::Real=1E-6, tol_μ::Real=1E-8, slow_updating::Real=1.0)
    """
    solve the economy with given liquidity multiplier ι
    """

    # solve household and banking problems
    crit_V = solve_value_and_pricing_function!(variables, parameters; tol=tol_h, iter_max=500, slow_updating=slow_updating)

    # solve the cross-sectional distribution
    crit_μ = solve_stationary_distribution_function!(variables, parameters; tol=tol_μ, iter_max=1000)

    # compute aggregate variables
    variables.aggregate_variables = solve_aggregate_variable_function(variables.policy_a, variables.threshold_a, variables.policy_pos_a, variables.q, variables.rbl, variables.μ, variables.aggregate_prices.K_λ, variables.aggregate_prices.w_λ, variables.aggregate_prices.ι_λ, parameters)

    # compute the difference between demand and supply sides
    ED_KL_to_D_ratio = variables.aggregate_variables.KL_to_D_ratio - variables.aggregate_prices.KL_to_D_ratio_λ
    ED_leverage_ratio = variables.aggregate_variables.leverage_ratio - variables.aggregate_prices.leverage_ratio_λ

    # printout results
    data_spec = Any[
        "Wage Garnishment Rate" parameters.η #=1=#
        "Probability of Bad History Removal" parameters.p_h #=2=#
        "Liquidity Multiplier" variables.aggregate_prices.λ #=3=#
        "Asset-to-Debt Ratio (Demand)" variables.aggregate_variables.KL_to_D_ratio #=4=#
        "Asset-to-Debt Ratio (Supply)" variables.aggregate_prices.KL_to_D_ratio_λ #=5=#
        "Difference" ED_KL_to_D_ratio #=6=#
        "Leverage Ratio (Demand)" variables.aggregate_variables.leverage_ratio #=7=#
        "Leverage Ratio (Supply)" variables.aggregate_prices.leverage_ratio_λ #=8=#
        "Difference" ED_leverage_ratio #=9=#
    ]
    pretty_table(data_spec; header=["Name", "Value"], alignment=[:l, :r], formatters=ft_round(8), body_hlines=[3, 5, 6, 8])

    # return excess demand
    return ED_KL_to_D_ratio, ED_leverage_ratio, crit_V, crit_μ
end

function optimal_multiplier_function(parameters::NamedTuple; λ_min_adhoc::Real=-Inf, λ_max_adhoc::Real=Inf, tol::Real=1E-5, iter_max::Real=200, slow_updating::Real=1.0)
    """
    solve for optimal liquidity multiplier
    """

    # check the case of λ_min = 0.0
    λ_min = 0.0
    variables_λ_min = variables_function(parameters; λ=λ_min)
    ED_KL_to_D_ratio_λ_min, ED_leverage_ratio_λ_min, crit_V_min, crit_μ_min = solve_economy_function!(variables_λ_min, parameters; slow_updating=slow_updating)
    # if ED_KL_to_D_ratio_λ_min > 0.0
    #     return variables_λ_min, variables_λ_min, 1
    # end
    if ED_leverage_ratio_λ_min < 0.0
        return variables_λ_min, variables_λ_min, 1, crit_V_min, crit_μ_min
    end

    # check the case of λ_max = 1-ψ^(1/2)
    λ_max = 1.0 - sqrt(parameters.ψ)
    variables_λ_max = variables_function(parameters; λ=λ_max)
    ED_KL_to_D_ratio_λ_max, ED_leverage_ratio_λ_max, crit_V_max, crit_μ_max = solve_economy_function!(variables_λ_max, parameters; slow_updating=slow_updating)
    # if ED_KL_to_D_ratio_λ_max < 0.0
    #     return variables_λ_min, variables_λ_max, 2 # meaning solution doesn't exist!
    # end
    if ED_leverage_ratio_λ_max > 0.0
        return variables_λ_min, variables_λ_max, 2, crit_V_max, crit_μ_max # meaning solution doesn't exist!
    end

    # initialization
    search_iter = 0
    crit = Inf
    λ_optimal = 0.0
    crit_V_optimal = 0.0
    crit_μ_optimal = 0.0
    variables_λ_optimal = []
    λ_lower = max(λ_min_adhoc, λ_min)
    λ_upper = min(λ_max_adhoc, λ_max)

    # solve equlibrium multiplier by bisection
    while crit > tol && search_iter < iter_max

        # update the multiplier
        λ_optimal = (λ_lower + λ_upper) / 2

        # compute the associated results
        # if search_iter == 0
        #     variables_λ_optimal = variables_function(parameters; λ = λ_optimal)
        # else
        #     variables_function_update!(variables_λ_optimal, parameters; λ = λ_optimal)
        # end
        variables_λ_optimal = variables_function(parameters; λ=λ_optimal)
        ED_KL_to_D_ratio_λ_optimal, ED_leverage_ratio_λ_optimal, crit_V_optimal, crit_μ_optimal = solve_economy_function!(variables_λ_optimal, parameters; slow_updating=slow_updating)

        # update search region
        # if ED_KL_to_D_ratio_λ_optimal > 0.0
        #     λ_upper = λ_optimal
        # else
        #     λ_lower = λ_optimal
        # end
        if ED_leverage_ratio_λ_optimal < 0.0
            λ_upper = λ_optimal
        else
            λ_lower = λ_optimal
        end

        # check convergence
        # crit = abs(ED_KL_to_D_ratio_λ_optimal)
        crit = abs(ED_leverage_ratio_λ_optimal)

        # update the iteration number
        search_iter += 1

    end

    # return results
    return variables_λ_min, variables_λ_optimal, 3, crit_V_optimal, crit_μ_optimal
end

function results_η_function(; η_min::Real, η_max::Real, η_step::Real, θ::Real=1.0 / (4.57 * 0.75))
    """
    compute stationary equilibrium with various η
    """

    # initialize η grid
    η_grid = collect(η_max:-η_step:η_min)
    η_size = length(η_grid)

    # initialize pparameters
    parameters = parameters_function(θ=θ)
    @unpack a_size, a_size_pos, a_size_μ, e_1_size, e_2_size, e_3_size, ν_size = parameters

    # initialize variables that will be saved
    var_names = [
        "Wage Garnishment Rate", #=1=#
        "Rental Rate", #=2=#
        "Liquidity Multiplier", #=3=#
        "Liquidity Premium", #=4=#
        "Capital", #=5=#
        "Loans", #=6=#
        "Deposits", #=7=#
        "Net Worth", #=8=#
        "Leverage Ratio", #=9=#
        "Share of Filers", #=10=#
        "Sahre in Debt", #=11=#
        "Debt-to-Earning Ratio", #=12=#
        "Average Loan Rate", #=13=#
        "Flag",  #=14=#
    ]
    var_size = length(var_names)

    # initialize containers
    results_A_NFF = zeros(var_size, η_size)
    results_V_NFF = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size, η_size)
    results_V_pos_NFF = zeros(a_size_pos, e_1_size, e_2_size, e_3_size, ν_size, η_size)
    results_μ_NFF = zeros(a_size_μ, e_1_size, e_2_size, e_3_size, ν_size, 2, η_size)
    results_A_FF = zeros(var_size, η_size)
    results_V_FF = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size, η_size)
    results_V_pos_FF = zeros(a_size_pos, e_1_size, e_2_size, e_3_size, ν_size, η_size)
    results_μ_FF = zeros(a_size_μ, e_1_size, e_2_size, e_3_size, ν_size, 2, η_size)

    # compute the optimal multipliers with different η
    for η_i = 1:η_size

        η = η_grid[η_i]
        parameters_η = parameters_function(η=η, θ=θ)
        # λ_min_adhoc_η = η_i > 1 ? results_A_FF[3,η_i-1] : -Inf
        λ_min_adhoc_η = -Inf
        variables_NFF, variables_FF, flag, crit_V, crit_μ = optimal_multiplier_function(parameters_η; λ_min_adhoc=λ_min_adhoc_η, slow_updating=slow_updating)

        # save results
        results_A_NFF[1, η_i] = parameters_η.η
        results_A_NFF[2, η_i] = variables_NFF.aggregate_prices.r_k_λ
        results_A_NFF[3, η_i] = variables_NFF.aggregate_prices.λ
        results_A_NFF[4, η_i] = variables_NFF.aggregate_prices.ι_λ
        results_A_NFF[5, η_i] = variables_NFF.aggregate_variables.K
        results_A_NFF[6, η_i] = variables_NFF.aggregate_variables.L
        results_A_NFF[7, η_i] = variables_NFF.aggregate_variables.D
        results_A_NFF[8, η_i] = variables_NFF.aggregate_variables.N
        results_A_NFF[9, η_i] = variables_NFF.aggregate_variables.leverage_ratio
        results_A_NFF[10, η_i] = variables_NFF.aggregate_variables.share_of_filers
        results_A_NFF[11, η_i] = variables_NFF.aggregate_variables.share_in_debts
        results_A_NFF[12, η_i] = variables_NFF.aggregate_variables.debt_to_earning_ratio
        results_A_NFF[13, η_i] = variables_NFF.aggregate_variables.avg_loan_rate
        results_A_NFF[14, η_i] = 1
        results_V_NFF[:, :, :, :, :, η_i] = variables_NFF.V
        results_V_pos_NFF[:, :, :, :, :, η_i] = variables_NFF.V_pos
        results_μ_NFF[:, :, :, :, :, :, η_i] = variables_NFF.μ

        results_A_FF[1, η_i] = parameters_η.η
        results_A_FF[2, η_i] = variables_FF.aggregate_prices.r_k_λ
        results_A_FF[3, η_i] = variables_FF.aggregate_prices.λ
        results_A_FF[4, η_i] = variables_FF.aggregate_prices.ι_λ
        results_A_FF[5, η_i] = variables_FF.aggregate_variables.K
        results_A_FF[6, η_i] = variables_FF.aggregate_variables.L
        results_A_FF[7, η_i] = variables_FF.aggregate_variables.D
        results_A_FF[8, η_i] = variables_FF.aggregate_variables.N
        results_A_FF[9, η_i] = variables_FF.aggregate_variables.leverage_ratio
        results_A_FF[10, η_i] = variables_FF.aggregate_variables.share_of_filers
        results_A_FF[11, η_i] = variables_FF.aggregate_variables.share_in_debts
        results_A_FF[12, η_i] = variables_FF.aggregate_variables.debt_to_earning_ratio
        results_A_FF[13, η_i] = variables_FF.aggregate_variables.avg_loan_rate
        results_A_FF[14, η_i] = flag
        results_V_FF[:, :, :, :, :, η_i] = variables_FF.V
        results_V_pos_FF[:, :, :, :, :, η_i] = variables_FF.V_pos
        results_μ_FF[:, :, :, :, :, :, η_i] = variables_FF.μ
    end

    # return results
    return var_names, results_A_NFF, results_V_NFF, results_V_pos_NFF, results_μ_NFF, results_A_FF, results_V_FF, results_V_pos_FF, results_μ_FF
end

function results_p_h_function(; p_h_min::Real, p_h_max::Real, p_h_step::Real, θ::Real=1.0 / (4.57 * 0.75))
    """
    compute stationary equilibrium with various p_h
    """

    # initialize p_h grid
    p_h_grid_ = collect(p_h_max:-p_h_step:p_h_min)
    p_h_grid = 1.0 ./ p_h_grid_
    p_h_size = length(p_h_grid)

    # initialize pparameters
    parameters = parameters_function(θ=θ)
    @unpack a_size, a_size_pos, a_size_μ, e_1_size, e_2_size, e_3_size, ν_size = parameters

    # initialize variables that will be saved
    var_names = [
        "Prob. of Bad Record Removal", #=1=#
        "Rental Rate", #=2=#
        "Liquidity Multiplier", #=3=#
        "Liquidity Premium", #=4=#
        "Capital", #=5=#
        "Loans", #=6=#
        "Deposits", #=7=#
        "Net Worth", #=8=#
        "Leverage Ratio", #=9=#
        "Share of Filers", #=10=#
        "Sahre in Debt", #=11=#
        "Debt-to-Earning Ratio", #=12=#
        "Average Loan Rate", #=13=#
        "Flag",  #=14=#
    ]
    var_size = length(var_names)

    # initialize containers
    results_A_NFF = zeros(var_size, p_h_size)
    results_V_NFF = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size, p_h_size)
    results_V_pos_NFF = zeros(a_size_pos, e_1_size, e_2_size, e_3_size, ν_size, p_h_size)
    results_μ_NFF = zeros(a_size_μ, e_1_size, e_2_size, e_3_size, ν_size, 2, p_h_size)
    results_A_FF = zeros(var_size, p_h_size)
    results_V_FF = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size, p_h_size)
    results_V_pos_FF = zeros(a_size_pos, e_1_size, e_2_size, e_3_size, ν_size, p_h_size)
    results_μ_FF = zeros(a_size_μ, e_1_size, e_2_size, e_3_size, ν_size, 2, p_h_size)

    # compute the optimal multipliers with different p_h
    for p_h_i = 1:p_h_size

        p_h = p_h_grid[p_h_i]
        parameters_p_h = parameters_function(p_h=p_h, θ=θ)
        # λ_min_adhoc_p_h = p_h_i > 1 ? results_A_FF[3,p_h_i-1] : -Inf
        λ_min_adhoc_p_h = -Inf
        variables_NFF, variables_FF, flag, crit_V, crit_μ = optimal_multiplier_function(parameters_p_h; λ_min_adhoc=λ_min_adhoc_p_h, slow_updating=slow_updating)

        # save results
        results_A_NFF[1, p_h_i] = parameters_p_h.p_h
        results_A_NFF[2, p_h_i] = variables_NFF.aggregate_prices.r_k_λ
        results_A_NFF[3, p_h_i] = variables_NFF.aggregate_prices.λ
        results_A_NFF[4, p_h_i] = variables_NFF.aggregate_prices.ι_λ
        results_A_NFF[5, p_h_i] = variables_NFF.aggregate_variables.K
        results_A_NFF[6, p_h_i] = variables_NFF.aggregate_variables.L
        results_A_NFF[7, p_h_i] = variables_NFF.aggregate_variables.D
        results_A_NFF[8, p_h_i] = variables_NFF.aggregate_variables.N
        results_A_NFF[9, p_h_i] = variables_NFF.aggregate_variables.leverage_ratio
        results_A_NFF[10, p_h_i] = variables_NFF.aggregate_variables.share_of_filers
        results_A_NFF[11, p_h_i] = variables_NFF.aggregate_variables.share_in_debts
        results_A_NFF[12, p_h_i] = variables_NFF.aggregate_variables.debt_to_earning_ratio
        results_A_NFF[13, p_h_i] = variables_NFF.aggregate_variables.avg_loan_rate
        results_A_NFF[14, p_h_i] = 1
        results_V_NFF[:, :, :, :, :, p_h_i] = variables_NFF.V
        results_V_pos_NFF[:, :, :, :, :, p_h_i] = variables_NFF.V_pos
        results_μ_NFF[:, :, :, :, :, :, p_h_i] = variables_NFF.μ

        results_A_FF[1, p_h_i] = parameters_p_h.p_h
        results_A_FF[2, p_h_i] = variables_FF.aggregate_prices.r_k_λ
        results_A_FF[3, p_h_i] = variables_FF.aggregate_prices.λ
        results_A_FF[4, p_h_i] = variables_FF.aggregate_prices.ι_λ
        results_A_FF[5, p_h_i] = variables_FF.aggregate_variables.K
        results_A_FF[6, p_h_i] = variables_FF.aggregate_variables.L
        results_A_FF[7, p_h_i] = variables_FF.aggregate_variables.D
        results_A_FF[8, p_h_i] = variables_FF.aggregate_variables.N
        results_A_FF[9, p_h_i] = variables_FF.aggregate_variables.leverage_ratio
        results_A_FF[10, p_h_i] = variables_FF.aggregate_variables.share_of_filers
        results_A_FF[11, p_h_i] = variables_FF.aggregate_variables.share_in_debts
        results_A_FF[12, p_h_i] = variables_FF.aggregate_variables.debt_to_earning_ratio
        results_A_FF[13, p_h_i] = variables_FF.aggregate_variables.avg_loan_rate
        results_A_FF[14, p_h_i] = flag
        results_V_FF[:, :, :, :, :, p_h_i] = variables_FF.V
        results_V_pos_FF[:, :, :, :, :, p_h_i] = variables_FF.V_pos
        results_μ_FF[:, :, :, :, :, :, p_h_i] = variables_FF.μ
    end

    # return results
    return var_names, results_A_NFF, results_V_NFF, results_V_pos_NFF, results_μ_NFF, results_A_FF, results_V_FF, results_V_pos_FF, results_μ_FF
end

function results_θ_function(; θ_min::Real, θ_max::Real, θ_step::Real)
    """
    compute stationary equilibrium with various θ
    """

    # initialize θ grid
    θ_grid = collect(θ_max:-θ_step:θ_min)
    θ_size = length(θ_grid)

    # initialize pparameters
    parameters = parameters_function()
    @unpack a_size, a_size_pos, a_size_μ, e_1_size, e_2_size, e_3_size, ν_size = parameters

    # initialize variables that will be saved
    var_names = [
        "Fraction of asset diversion", #=1=#
        "Rental Rate", #=2=#
        "Liquidity Multiplier", #=3=#
        "Liquidity Premium", #=4=#
        "Capital", #=5=#
        "Loans", #=6=#
        "Deposits", #=7=#
        "Net Worth", #=8=#
        "Leverage Ratio", #=9=#
        "Share of Filers", #=10=#
        "Sahre in Debt", #=11=#
        "Debt-to-Earning Ratio", #=12=#
        "Average Loan Rate", #=13=#
        "Flag",  #=14=#
    ]
    var_size = length(var_names)

    # initialize containers
    results_A_NFF = zeros(var_size, θ_size)
    results_V_NFF = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size, θ_size)
    results_V_pos_NFF = zeros(a_size_pos, e_1_size, e_2_size, e_3_size, ν_size, θ_size)
    results_μ_NFF = zeros(a_size_μ, e_1_size, e_2_size, e_3_size, ν_size, 2, θ_size)
    results_A_FF = zeros(var_size, θ_size)
    results_V_FF = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size, θ_size)
    results_V_pos_FF = zeros(a_size_pos, e_1_size, e_2_size, e_3_size, ν_size, θ_size)
    results_μ_FF = zeros(a_size_μ, e_1_size, e_2_size, e_3_size, ν_size, 2, θ_size)

    # compute the optimal multipliers with different θ
    for θ_i = 1:θ_size

        θ = θ_grid[θ_i]
        parameters_θ = parameters_function(θ=θ)
        # λ_min_adhoc_θ = θ_i > 1 ? results_A_FF[3,θ_i-1] : -Inf
        λ_min_adhoc_θ = -Inf
        variables_NFF, variables_FF, flag, crit_V, crit_μ = optimal_multiplier_function(parameters_θ; λ_min_adhoc=λ_min_adhoc_θ, slow_updating=slow_updating)

        # save results
        results_A_NFF[1, θ_i] = parameters_θ.θ
        results_A_NFF[2, θ_i] = variables_NFF.aggregate_prices.r_k_λ
        results_A_NFF[3, θ_i] = variables_NFF.aggregate_prices.λ
        results_A_NFF[4, θ_i] = variables_NFF.aggregate_prices.ι_λ
        results_A_NFF[5, θ_i] = variables_NFF.aggregate_variables.K
        results_A_NFF[6, θ_i] = variables_NFF.aggregate_variables.L
        results_A_NFF[7, θ_i] = variables_NFF.aggregate_variables.D
        results_A_NFF[8, θ_i] = variables_NFF.aggregate_variables.N
        results_A_NFF[9, θ_i] = variables_NFF.aggregate_variables.leverage_ratio
        results_A_NFF[10, θ_i] = variables_NFF.aggregate_variables.share_of_filers
        results_A_NFF[11, θ_i] = variables_NFF.aggregate_variables.share_in_debts
        results_A_NFF[12, θ_i] = variables_NFF.aggregate_variables.debt_to_earning_ratio
        results_A_NFF[13, θ_i] = variables_NFF.aggregate_variables.avg_loan_rate
        results_A_NFF[14, θ_i] = 1
        results_V_NFF[:, :, :, :, :, θ_i] = variables_NFF.V
        results_V_pos_NFF[:, :, :, :, :, θ_i] = variables_NFF.V_pos
        results_μ_NFF[:, :, :, :, :, :, θ_i] = variables_NFF.μ

        results_A_FF[1, θ_i] = parameters_θ.θ
        results_A_FF[2, θ_i] = variables_FF.aggregate_prices.r_k_λ
        results_A_FF[3, θ_i] = variables_FF.aggregate_prices.λ
        results_A_FF[4, θ_i] = variables_FF.aggregate_prices.ι_λ
        results_A_FF[5, θ_i] = variables_FF.aggregate_variables.K
        results_A_FF[6, θ_i] = variables_FF.aggregate_variables.L
        results_A_FF[7, θ_i] = variables_FF.aggregate_variables.D
        results_A_FF[8, θ_i] = variables_FF.aggregate_variables.N
        results_A_FF[9, θ_i] = variables_FF.aggregate_variables.leverage_ratio
        results_A_FF[10, θ_i] = variables_FF.aggregate_variables.share_of_filers
        results_A_FF[11, θ_i] = variables_FF.aggregate_variables.share_in_debts
        results_A_FF[12, θ_i] = variables_FF.aggregate_variables.debt_to_earning_ratio
        results_A_FF[13, θ_i] = variables_FF.aggregate_variables.avg_loan_rate
        results_A_FF[14, θ_i] = flag
        results_V_FF[:, :, :, :, :, θ_i] = variables_FF.V
        results_V_pos_FF[:, :, :, :, :, θ_i] = variables_FF.V_pos
        results_μ_FF[:, :, :, :, :, :, θ_i] = variables_FF.μ
    end

    # return results
    return var_names, results_A_NFF, results_V_NFF, results_V_pos_NFF, results_μ_NFF, results_A_FF, results_V_FF, results_V_pos_FF, results_μ_FF
end

function results_ψ_function(; ψ_min::Real, ψ_max::Real, ψ_step::Real)
    """
    compute stationary equilibrium with various ψ
    """

    # initialize ψ grid
    ψ_grid_ = collect(ψ_max:-ψ_step:ψ_min)
    ψ_grid = 1.0 .- 1.0 ./ ψ_grid_
    ψ_size = length(ψ_grid)

    # initialize pparameters
    parameters = parameters_function()
    @unpack a_size, a_size_pos, a_size_μ, e_1_size, e_2_size, e_3_size, ν_size = parameters

    # initialize variables that will be saved
    var_names = [
        "Planning horizon of banks", #=1=#
        "Rental Rate", #=2=#
        "Liquidity Multiplier", #=3=#
        "Liquidity Premium", #=4=#
        "Capital", #=5=#
        "Loans", #=6=#
        "Deposits", #=7=#
        "Net Worth", #=8=#
        "Leverage Ratio", #=9=#
        "Share of Filers", #=10=#
        "Sahre in Debt", #=11=#
        "Debt-to-Earning Ratio", #=12=#
        "Average Loan Rate", #=13=#
        "Flag",  #=14=#
    ]
    var_size = length(var_names)

    # initialize containers
    results_A_NFF = zeros(var_size, ψ_size)
    results_V_NFF = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size, ψ_size)
    results_V_pos_NFF = zeros(a_size_pos, e_1_size, e_2_size, e_3_size, ν_size, ψ_size)
    results_μ_NFF = zeros(a_size_μ, e_1_size, e_2_size, e_3_size, ν_size, 2, ψ_size)
    results_A_FF = zeros(var_size, ψ_size)
    results_V_FF = zeros(a_size, e_1_size, e_2_size, e_3_size, ν_size, ψ_size)
    results_V_pos_FF = zeros(a_size_pos, e_1_size, e_2_size, e_3_size, ν_size, ψ_size)
    results_μ_FF = zeros(a_size_μ, e_1_size, e_2_size, e_3_size, ν_size, 2, ψ_size)

    # compute the optimal multipliers with different ψ
    for ψ_i = 1:ψ_size

        ψ = ψ_grid[ψ_i]
        parameters_ψ = parameters_function(ψ=ψ)
        # λ_min_adhoc_ψ = ψ_i > 1 ? results_A_FF[3,ψ_i-1] : -Inf
        λ_min_adhoc_ψ = -Inf
        variables_NFF, variables_FF, flag, crit_V, crit_μ = optimal_multiplier_function(parameters_ψ; λ_min_adhoc=λ_min_adhoc_ψ, slow_updating=slow_updating)

        # save results
        results_A_NFF[1, ψ_i] = parameters_ψ.ψ
        results_A_NFF[2, ψ_i] = variables_NFF.aggregate_prices.r_k_λ
        results_A_NFF[3, ψ_i] = variables_NFF.aggregate_prices.λ
        results_A_NFF[4, ψ_i] = variables_NFF.aggregate_prices.ι_λ
        results_A_NFF[5, ψ_i] = variables_NFF.aggregate_variables.K
        results_A_NFF[6, ψ_i] = variables_NFF.aggregate_variables.L
        results_A_NFF[7, ψ_i] = variables_NFF.aggregate_variables.D
        results_A_NFF[8, ψ_i] = variables_NFF.aggregate_variables.N
        results_A_NFF[9, ψ_i] = variables_NFF.aggregate_variables.leverage_ratio
        results_A_NFF[10, ψ_i] = variables_NFF.aggregate_variables.share_of_filers
        results_A_NFF[11, ψ_i] = variables_NFF.aggregate_variables.share_in_debts
        results_A_NFF[12, ψ_i] = variables_NFF.aggregate_variables.debt_to_earning_ratio
        results_A_NFF[13, ψ_i] = variables_NFF.aggregate_variables.avg_loan_rate
        results_A_NFF[14, ψ_i] = 1
        results_V_NFF[:, :, :, :, :, ψ_i] = variables_NFF.V
        results_V_pos_NFF[:, :, :, :, :, ψ_i] = variables_NFF.V_pos
        results_μ_NFF[:, :, :, :, :, :, ψ_i] = variables_NFF.μ

        results_A_FF[1, ψ_i] = parameters_ψ.ψ
        results_A_FF[2, ψ_i] = variables_FF.aggregate_prices.r_k_λ
        results_A_FF[3, ψ_i] = variables_FF.aggregate_prices.λ
        results_A_FF[4, ψ_i] = variables_FF.aggregate_prices.ι_λ
        results_A_FF[5, ψ_i] = variables_FF.aggregate_variables.K
        results_A_FF[6, ψ_i] = variables_FF.aggregate_variables.L
        results_A_FF[7, ψ_i] = variables_FF.aggregate_variables.D
        results_A_FF[8, ψ_i] = variables_FF.aggregate_variables.N
        results_A_FF[9, ψ_i] = variables_FF.aggregate_variables.leverage_ratio
        results_A_FF[10, ψ_i] = variables_FF.aggregate_variables.share_of_filers
        results_A_FF[11, ψ_i] = variables_FF.aggregate_variables.share_in_debts
        results_A_FF[12, ψ_i] = variables_FF.aggregate_variables.debt_to_earning_ratio
        results_A_FF[13, ψ_i] = variables_FF.aggregate_variables.avg_loan_rate
        results_A_FF[14, ψ_i] = flag
        results_V_FF[:, :, :, :, :, ψ_i] = variables_FF.V
        results_V_pos_FF[:, :, :, :, :, ψ_i] = variables_FF.V_pos
        results_μ_FF[:, :, :, :, :, :, ψ_i] = variables_FF.μ
    end

    # return results
    return var_names, results_A_NFF, results_V_NFF, results_V_pos_NFF, results_μ_NFF, results_A_FF, results_V_FF, results_V_pos_FF, results_μ_FF
end

# function results_CEV_function(parameters::NamedTuple, results_A::Array{Float64,2}, results_V::Array{Float64,6}, results_V_pos::Array{Float64,6})
#     """
#     compute consumption equivalent variation (CEV) with various η compared to the smallest η (most lenient policy)
#     """
#
#     # initialize pparameters
#     @unpack a_grid, a_size, a_grid_pos, a_size_pos, a_grid_μ, a_size_μ, a_ind_zero_μ, e_1_size, e_2_size, e_3_size, ν_size, σ = parameters
#
#     # initialize result matrix
#     η_size = size(results_A)[2]
#     results_CEV = zeros(a_size_μ, e_1_size, e_2_size, e_3_size, ν_size, 2, η_size)
#
#     # compute CEV for different η compared to the smallest η
#     for η_i = 1:η_size, e_1_i = 1:e_1_size, e_2_i = 1:e_2_size, e_3_i = 1:e_3_size, ν_i = 1:ν_size
#         @inbounds @views V_itp_new = Akima(a_grid, results_V[:, e_1_i, e_2_i, e_3_i, ν_i, η_i])
#         @inbounds @views V_itp_old = Akima(a_grid, results_V[:, e_1_i, e_2_i, e_3_i, ν_i, end])
#         @inbounds @views V_pos_itp_new = Akima(a_grid_pos, results_V_pos[:, e_1_i, e_2_i, e_3_i, ν_i, η_i])
#         @inbounds @views V_pos_itp_old = Akima(a_grid_pos, results_V_pos[:, e_1_i, e_2_i, e_3_i, ν_i, end])
#         for a_μ_i = 1:a_size_μ
#             @inbounds a_μ = a_grid_μ[a_μ_i]
#             @inbounds results_CEV[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1, η_i] = (V_itp_new(a_μ) / V_itp_old(a_μ))^(1.0 / (1.0 - σ)) - 1.0
#             if a_μ >= 0.0
#                 # a_pos_μ_i = a_μ_i - a_ind_zero_μ + 1
#                 @inbounds results_CEV[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 2, η_i] = (V_pos_itp_new(a_μ) / V_pos_itp_old(a_μ))^(1.0 / (1.0 - σ)) - 1.0
#             end
#         end
#     end
#
#     # return results
#     return results_CEV
# end

# function results_HHs_favor_function(parameters::NamedTuple, results_A::Array{Float64,2}, results_V::Array{Float64,6}, results_V_pos::Array{Float64,6})
#     """
#     compute consumption share of HHs in favor of the given η compared to the smallest η (most lenient policy)
#     """
#
#     # initialize pparameters
#     @unpack a_grid, a_size, a_grid_pos, a_size_pos, a_grid_μ, a_size_μ, a_ind_zero_μ, e_1_size, e_2_size, e_3_size, ν_size, σ = parameters
#
#     # initialize result matrix
#     η_size = size(results_A)[2]
#     results_HHs_favor = zeros(a_size_μ, e_1_size, e_2_size, e_3_size, ν_size, 2, η_size)
#
#     # compute CEV for different η compared to the smallest η
#     for η_i = 1:η_size, ν_i = 1:ν_size, e_3_i = 1:e_3_size, e_2_i = 1:e_2_size, e_1_i = 1:e_1_size
#         @inbounds @views V_itp_new = Akima(a_grid, results_V[:, e_1_i, e_2_i, e_3_i, ν_i, η_i])
#         @inbounds @views V_itp_old = Akima(a_grid, results_V[:, e_1_i, e_2_i, e_3_i, ν_i, end])
#         @inbounds @views V_pos_itp_new = Akima(a_grid_pos, results_V_pos[:, e_1_i, e_2_i, e_3_i, ν_i, η_i])
#         @inbounds @views V_pos_itp_old = Akima(a_grid_pos, results_V_pos[:, e_1_i, e_2_i, e_3_i, ν_i, end])
#         for a_μ_i = 1:a_size_μ
#             @inbounds a_μ = a_grid_μ[a_μ_i]
#             @inbounds results_HHs_favor[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 1, η_i] = V_itp_new(a_μ) >= V_itp_old(a_μ)
#             if a_μ >= 0.0
#                 # a_pos_μ_i = a_μ_i - a_ind_zero_μ + 1
#                 @inbounds results_HHs_favor[a_μ_i, e_1_i, e_2_i, e_3_i, ν_i, 2, η_i] = V_pos_itp_new(a_μ) >= V_pos_itp_old(a_μ)
#             end
#         end
#     end
#
#     # return results
#     return results_HHs_favor
# end
