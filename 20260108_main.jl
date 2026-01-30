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

function adda_cooper(N::Integer, ρ::Real, σ::Real; μ::Real = 0.0)
	"""
	Approximation of an autoregression process with a Markov chain proposed by Adda and Cooper (2003)
	"""
	σ_ϵ = σ / sqrt(1.0 - ρ^2.0)
	q = quantile(Normal(), range(0, 1; length = N + 1))
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
	Π = Π ./ sum(Π, dims = 2)
	return z, Π
end

function tauchen_grid(N::Int, ρ::Float64, σ::Float64; μ::Float64 = 0.0, m::Float64 = 3.0)
	σ_z = σ / sqrt(1 - ρ^2)
	z_max = μ + m * σ_z
	z_min = μ - m * σ_z
	z = collect(range(z_min, z_max; length = N))
	return z
end

function tauchen_transition_matrix(z::Vector{Float64}, ρ::Float64, σ::Float64; μ::Float64 = 0.0)
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

function infertility_risk_function(data_age::Array{Int64, 1}, data_inf::Array{Float64, 1}, age_min::Integer, age_max::Integer, age_inf::Integer)
	"""
	Exponential fit of infertility probability, intrapolated on ages up to age_inf
	"""
	model(t, ω) = ω[1] * exp.(ω[2] * t)
	ω_int = [0.5, 0.5]
	fit = curve_fit(model, data_age, data_inf, ω_int)
	model_age = collect(age_min:age_max)
	model_inf = fit.param[1] .* exp.(fit.param[2] .* model_age)
	model_inf[findall(model_age .== age_inf)[]:end] .= 1.0
	return model_age, model_inf
end

function infertility_risk_low_function(data_age::Array{Int64, 1}, data_inf::Array{Float64, 1}, age_min::Integer, age_max::Integer)
	"""
	Exponential fit of infertility probability, intrapolated on ages up to age_inf
	"""
	model(t, ω) = ω[1] * exp.(ω[2] * t)
	ω_int = [0.5, 0.5]
	fit = curve_fit(model, data_age, data_inf, ω_int)
	model_age = collect(age_min:age_max)
	model_inf = fit.param[1] .* exp.(fit.param[2] .* model_age)
	age_inf = model_age[findlast(model_inf .< 1.0)[]]
	model_inf[findall(model_age .== (age_inf+1))[]:end] .= 1.0
	return model_age, model_inf, age_inf
end

function h_function(data_h::Array{Float64, 1}, age_min::Integer, age_ret::Integer)
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

@inline function u_CRRA_e(c::Float64, e::Float64, one_m_γ::Float64, inv_one_m_γ::Float64, one_m_κ::Float64, ψ_inv_one_m_κ::Float64, scale::Float64; c_floor::Float64 = 1e-12, V_penalty::Float64 = -1.0e16)
	c <= c_floor && return V_penalty
	return c^one_m_γ * inv_one_m_γ + scale * e^one_m_κ * ψ_inv_one_m_κ
end

@inline function u_CRRA(c::Float64, one_m_γ::Float64, inv_one_m_γ::Float64; c_floor::Float64 = 1e-12, V_penalty::Float64 = -1.0e16)
	c <= c_floor && return V_penalty
	return c^one_m_γ * inv_one_m_γ
end

@inline function u_log(c::Float64; c_floor::Float64 = 1e-12, V_penalty::Float64 = -1.0e16)
	c <= c_floor && return V_penalty
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
	r::Real = 0.04,                   # interest rate
	β::Real = 1.0 / (1.0 + r),        # discount factor
	γ::Real = 1.5,                    # risk aversion
	ρ::Real = 0.95,                   # persistance coefficient
	σ_ϵ::Real = 0.21,                 # std of persistent shock
	σ_ν::Real = 0.17,                 # std of transitory shock
	b::Real = 0.40,                   # replacement rate
	#----------------------#
	# estimated parameters #
	#----------------------#
	κ::Real = 0.14,                   # preference curvature
	ψ::Real = 3.50,                   # preference scale
	μ::Real = 0.35,                   # production share
	θ::Real = 0.70,                   # elasticity of substitution in production
	q_bar::Real = 1.50,               # lower bound on children's consumption
	ψ_1::Real = 0.91,                 # HH economies to money input to production
	ψ_2::Real = 0.54,                 # HH economies to time input to production
	p::Real = 0.02, #====================##====================#                   # prob that a child becomes independent

	# numerical solution #

	age_min::Integer = 18,            # min age
	age_max::Integer = 80,            # max age
	age_edu::Integer = 22,            # education age
	age_inf::Integer = 45,            # infertile age
	age_ret::Integer = 65,            # retirement age
	n_max::Integer = 4,               # max number of kids
	ϵ_size::Integer = 7,              # number of persistent shock
	ν_size::Integer = 5,              # number of transitory shock
	a_max::Real = 120,                # max of asset holding
	a_size_neg::Integer = 11,         # number of negative asset
	a_size::Integer = 50,             # number of asset
	a_degree::Integer = 2,            # curvature of asset gridpoints
	q_x::Real = 1.0, #=================##=================#                  # price of monetary input $x$

	# case indicators #

	inf_scale::Real = 1.0,            # scale of infertility risk
	edu_h_ind::Real = 0.0,            # education indicator
	edu_h_scale::Real = 1.0,          # scale of edu-dependent life-cycle income
	edu_σ_ϵ_scale::Real = 1.0,        # scale of persistent income uncertainty    
	edu_σ_ν_scale::Real = 1.0,        # scale of transitory income uncertainty
	lr_exp::Integer = 0,               # long-run response experiment
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
	n_grid = collect(0.0:n_max)
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
		2.509402,
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
	ϵ_grid = tauchen_grid(ϵ_size, ρ, σ_ϵ; m = 2.0)
	ϵ_Γ = tauchen_transition_matrix(ϵ_grid, ρ, σ_ϵ * edu_σ_ϵ_scale)
	# ϵ_grid, ϵ_Γ = adda_cooper(ϵ_size, ρ, σ_ϵ * edu_σ_ϵ_scale)
	ϵ_G = stationary_distributions(MarkovChain(ϵ_Γ, ϵ_grid))[1]

	# transitory income shock
	ν_grid = tauchen_grid(ν_size, 0.0, σ_ν; m = 2.0)
	ν_Γ = tauchen_transition_matrix(ν_grid, 0.0, σ_ν * edu_σ_ν_scale)
	# ν_grid, ν_Γ = adda_cooper(ν_size, 0.0, σ_ν * edu_σ_ν_scale)
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

	Γ = zeros(n_size, ϵ_size, ν_size, inf_size, n_size, ϵ_size, inf_size, age_inf-age_min)
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
	if edu_h_ind == 0.0
		a_grid = ((range(0.0, stop = a_size - 1, length = a_size) / (a_size - 1)) .^ a_degree) * a_max
		a_ind_zero = 1
	else
		a_grid_neg = collect(range(a_min, 0.0, length = a_size_neg))
		a_grid_pos = ((range(0.0, stop = a_size - 1, length = a_size) / (a_size - 1)) .^ a_degree) * a_max
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

	# numerical guards
	c_floor   = 1.0e-12
	V_penalty = -1.0e16

	# return values
	return (
		r = r,
		R = R,
		inv_R = inv_R,
		β = β,
		γ = γ,
		m_inv_γ = m_inv_γ,
		one_m_γ = one_m_γ,
		inv_one_m_γ = inv_one_m_γ,
		βR = βR,
		EGM_fac = EGM_fac,
		one_m_κ = one_m_κ,
		inv_κ = inv_κ,
		inv_one_m_κ = inv_one_m_κ,
		ψ_inv_one_m_κ = ψ_inv_one_m_κ,
		γ_by_κ = γ_by_κ,
		ρ = ρ,
		σ_ϵ = σ_ϵ,
		σ_ν = σ_ν,
		b = b,
		κ = κ,
		ψ = ψ,
		μ = μ,
		θ = θ,
		q_bar = q_bar,
		ψ_1 = ψ_1,
		ψ_2 = ψ_2,
		p = p,
		age_min = age_min,
		age_max = age_max,
		age_edu = age_edu,
		age_inf = age_inf,
		age_ret = age_ret,
		age_size = age_size,
		age_grid = age_grid,
		inf_grid = inf_grid,
		data_age = data_age,
		data_inf = data_inf,
		inf_size = inf_size,
		d_κ = d_κ,
		d_ι = d_ι,
		n_max = n_max,
		n_size = n_size,
		n_grid = n_grid,
		n_Γ = n_Γ,
		h_size = h_size,
		h_grid = h_grid,
		data_h = data_h,
		ϵ_size = ϵ_size,
		ϵ_grid = ϵ_grid,
		ϵ_Γ = ϵ_Γ,
		ϵ_G = ϵ_G,
		ν_size = ν_size,
		ν_grid = ν_grid,
		ν_Γ = ν_Γ,
		ν_G = ν_G,
		Γ_ret = Γ_ret,
		Γ_inf = Γ_inf,
		Γ = Γ,
		w_grid = w_grid,
		P_grid = P_grid,
		q_bar_P_grid = q_bar_P_grid,
		a_min = a_min,
		a_max = a_max,
		a_ind_zero = a_ind_zero,
		a_size = a_size,
		a_size_neg = a_size_neg,
		a_grid = a_grid,
		aR_grid = aR_grid,
		a_degree = a_degree,
		l_size = l_size,
		l_grid = l_grid,
		x_size = x_size,
		x_grid = x_grid,
		q_x = q_x,
		inf_scale = inf_scale,
		edu_h_ind = edu_h_ind,
		edu_h_scale = edu_h_scale,
		edu_σ_ϵ_scale = edu_σ_ϵ_scale,
		edu_σ_ν_scale = edu_σ_ν_scale,
		c_floor = c_floor,
		V_penalty = V_penalty,
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
	EV_next::Array{Float64, 4}
	c_next_CE::Array{Float64, 4}
end

mutable struct Mutable_Variables
	"""
	Construct a type for mutable variables
	"""
	V::Array{Float64, 6}
	policy_c::Array{Float64, 6}
	policy_a_p::Array{Float64, 6}
	policy_e::Array{Float64, 6}
	policy_K::Array{Float64, 6}
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
			j  = searchsortedlast(x, xq)
			x0 = x[j];
			x1 = x[j+1]
			y0 = y[j];
			y1 = y[j+1]
			w  = (xq - x0) / (x1 - x0)
			return muladd(w, (y1 - y0), y0)
		end
	end
end

# @inline function terminal_step!(
# 	V_end::AbstractArray{Float64, 3},
# 	policy_c_end::AbstractArray{Float64, 3},
# 	parameters::NamedTuple,
# 	ϵ_i::Integer,
# 	ν_i::Integer,
# )
# 	V_end_a = @views V_end[:, ϵ_i, ν_i]
# 	policy_c_end_a = @views policy_c_end[:, ϵ_i, ν_i]

# 	@unpack w_grid, aR_grid, one_m_γ, inv_one_m_γ = parameters
# 	@inbounds w_bar = w_grid[ϵ_i, ν_i, end]

# 	@inbounds for i in eachindex(aR_grid)
# 		c = aR_grid[i] + w_bar
# 		policy_c_end_a[i] = c
# 		V_end_a[i] = u_CRRA(c, one_m_γ, inv_one_m_γ)
# 	end

# 	return nothing
# end

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
		V_end_a[i] = u_CRRA(c, one_m_γ, inv_one_m_γ; c_floor = c_floor, V_penalty = V_penalty)
	end
	return nothing
end

function retired_step!(
	V_endo::Vector{Float64},
	a_endo::Vector{Float64},
	ap_endo::Vector{Float64},
	V_next_a::AbstractVector{Float64},
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
		Vp = V_next_a[ap_i]
		if !isfinite(Vp) || Vp <= V_penalty
			continue
		end
		cnext = c_next_a[ap_i]
		if !isfinite(cnext) || cnext <= c_floor
			continue
		end
		c = EGM_fac * cnext
		if !isfinite(c) || c <= c_floor
			continue
		end
		ap      = a_grid[ap_i]
		m       = c + ap
		a_today = (m - w_bar) * inv_R
		if !isfinite(a_today)
			continue
		end
		nv          += 1
		V_endo[nv]  = Vp
		ap_endo[nv] = ap
		a_endo[nv]  = a_today
	end

	if nv == 0
		@inbounds for a_i in 1:a_size
			a_p_current_a[a_i] = a_min
			c_current_a[a_i] = c_floor
			V_current_a[a_i] = V_penalty
		end
		return nothing
	end

	V_endo_v  = @view V_endo[1:nv]
	a_endo_v  = @view a_endo[1:nv]
	ap_endo_v = @view ap_endo[1:nv]

	@inbounds for j in 2:nv
		if a_endo_v[j] <= a_endo_v[j-1]
			a_endo_v[j] = nextfloat(a_endo_v[j-1])
		end
	end

	ibind = searchsortedlast(a_grid, a_endo_v[1])

	if ibind > 0
		Vp_amin = V_next_a[1]
	
		@inbounds for a_i in 1:ibind
			aR = aR_grid[a_i]
			ap = a_min
			c  = aR + w_bar - ap
	
			if !isfinite(c) || c <= c_floor || !isfinite(Vp_amin) || Vp_amin <= V_penalty
				a_p_current_a[a_i] = ap
				c_current_a[a_i]   = c_floor
				V_current_a[a_i]   = V_penalty
			else
				a_p_current_a[a_i] = ap
				c_current_a[a_i]   = c
				V_current_a[a_i]   = u_CRRA(c, one_m_γ, inv_one_m_γ; c_floor=c_floor, V_penalty=V_penalty) + β * Vp_amin
			end
		end
	end

	@inbounds for a_i in (ibind+1):a_size
		a  = a_grid[a_i]
		aR = aR_grid[a_i]
		ap = lininterp(a_endo_v, ap_endo_v, a)
		c  = aR + w_bar - ap
		if !isfinite(c) || c <= c_floor
			a_p_current_a[a_i] = ap
			c_current_a[a_i] = c_floor
			V_current_a[a_i] = V_penalty
		else
			a_p_current_a[a_i] = ap
			c_current_a[a_i] = c
			Vp = lininterp(ap_endo_v, V_endo_v, ap)
			if !isfinite(Vp) || Vp <= V_penalty
				c_current_a[a_i] = c_floor
				V_current_a[a_i] = V_penalty
				continue
			end
			V_current_a[a_i] = u_CRRA(c, one_m_γ, inv_one_m_γ; c_floor = c_floor, V_penalty = V_penalty) + β * Vp
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
	c_floor::Float64 = 1e-12,
	maxit::Int64 = 80,
	tol::Float64 = 1e-12,
)::Float64

	# auxiliary parameter
	A = ψ * (n / P)^one_m_κ

	# check feasibility
	e_hi = m - c_floor
	if e_hi < e_bar
		return e_hi
	end

	# Check lower bound (constraint) corner: if g(e_bar) <= 0, optimum at e = e_bar
	gl = g_eval(e_bar, m, A, γ, κ)
	if gl <= 0.0
		return e_bar
	end

	# Check upper bound corner: if g(e_hi) >= 0, utility still wants more e -> choose e_hi
	gh = g_eval(e_hi, m, A, γ, κ)
	if gh >= 0.0
		return e_hi
	end

	# Now we have g(e_min) > 0 and g(e_hi) < 0 -> unique root inside (e_min, e_hi)
	lo = e_bar
	hi = e_hi
	mid = 0.5 * (lo + hi)
	@inbounds for _ in 1:maxit
		mid = 0.5 * (lo + hi)
		gm = g_eval(mid, m, A, γ, κ)
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
	V_endo::Vector{Float64},
	a_endo::Vector{Float64},
	ap_endo::Vector{Float64},
	V_next_a::AbstractVector{Float64},
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
			V_current_a[a_i]   = V_penalty
			c_current_a[a_i]   = c_floor
			a_p_current_a[a_i] = a_min
			e_current_a[a_i]   = e_bar
		end
		return nothing
	end

	nv = 0
	@inbounds for ap_i in 1:a_size
		Vp = V_next_a[ap_i]
		if !isfinite(Vp) || Vp <= V_penalty
			continue
		end
		cnext = c_next_a[ap_i]
		if !isfinite(cnext) || cnext <= c_floor
			continue
		end
		c = EGM_fac * cnext
		if !isfinite(c) || c <= c_floor
			continue
		end
		e_foc = e_para * c^γ_by_κ
		e     = max(e_foc, e_bar)
		if !isfinite(e)
			continue
		end
		ap      = a_grid[ap_i]
		m       = c + ap + e
		a_today = (m - w_bar) * inv_R
		if !isfinite(a_today)
			continue
		end
		nv          += 1
		V_endo[nv]  = Vp
		ap_endo[nv] = ap
		a_endo[nv]  = a_today
	end

	if nv == 0
		@inbounds for a_i in 1:a_size
			V_current_a[a_i]   = V_penalty
			c_current_a[a_i]   = c_floor
			a_p_current_a[a_i] = a_min
			e_current_a[a_i]   = e_bar
		end
		return nothing
	end

	V_endo_v  = @view V_endo[1:nv]
	a_endo_v  = @view a_endo[1:nv]
	ap_endo_v = @view ap_endo[1:nv]

	@inbounds for j in 2:nv
		if a_endo_v[j] <= a_endo_v[j-1]
			a_endo_v[j] = nextfloat(a_endo_v[j-1])
		end
	end

	ibind = searchsortedlast(a_grid, a_endo_v[1])

	if ibind > 0
		Vp_amin = V_next_a[1]
	
		@inbounds for a_i in 1:ibind
			aR = aR_grid[a_i]
			M  = aR + w_bar
			ap = a_min
			m  = M - ap
	
			# static feasibility: must afford e_bar + c_floor
			if m <= e_bar + c_floor || !isfinite(Vp_amin) || Vp_amin <= V_penalty
				V_current_a[a_i]   = V_penalty
				c_current_a[a_i]   = c_floor
				a_p_current_a[a_i] = ap
				e_current_a[a_i]   = e_bar
				continue
			end
	
			e = solve_e_bisect(m, n, P, ψ, γ, κ, one_m_κ, e_bar)
			if !isfinite(e) || e < e_bar
				V_current_a[a_i]   = V_penalty
				c_current_a[a_i]   = c_floor
				a_p_current_a[a_i] = ap
				e_current_a[a_i]   = e_bar
				continue
			end
	
			c = m - e
			if !isfinite(c) || c <= c_floor
				V_current_a[a_i]   = V_penalty
				c_current_a[a_i]   = c_floor
				a_p_current_a[a_i] = ap
				e_current_a[a_i]   = e
				continue
			end
	
			V_current_a[a_i] = u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ, scale;
										c_floor=c_floor, V_penalty=V_penalty) + β * Vp_amin
			c_current_a[a_i]   = c
			a_p_current_a[a_i] = ap
			e_current_a[a_i]   = e
		end
	end	

	@inbounds for a_i in (ibind+1):a_size
		a = a_grid[a_i]
		aR = aR_grid[a_i]
		ap = lininterp(a_endo_v, ap_endo_v, a)
		M = aR + w_bar
		m = M - ap
		if m <= e_bar + c_floor
			V_current_a[a_i]   = V_penalty
			c_current_a[a_i]   = c_floor
			a_p_current_a[a_i] = ap
			e_current_a[a_i]   = e_bar
		else
			Vp = lininterp(ap_endo_v, V_endo_v, ap)
			if !isfinite(Vp) || Vp <= V_penalty
				V_current_a[a_i]   = V_penalty
				c_current_a[a_i]   = c_floor
				a_p_current_a[a_i] = ap
				e_current_a[a_i]   = e_bar
				continue
			end
			e = solve_e_bisect(m, n, P, ψ, γ, κ, one_m_κ, e_bar)
			if !isfinite(e) || e < e_bar
				V_current_a[a_i]   = V_penalty
				c_current_a[a_i]   = c_floor
				a_p_current_a[a_i] = ap
				e_current_a[a_i]   = e_bar
				continue
			end
			c = m - e
			V_current_a[a_i] = u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ, scale; c_floor = c_floor, V_penalty = V_penalty) + β * Vp
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
				V_current_a[i]   = V1[i]
				c_current_a[i]   = c1[i]
				a_p_current_a[i] = ap1[i]
				K_current_a[i]   = 1.0
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
				V_current_a[i]   = V1[i]
				c_current_a[i]   = c1[i]
				a_p_current_a[i] = ap1[i]
				e_current_a[i]   = e1[i]
				K_current_a[i]   = 1.0
			end
		end
	end
	return nothing
end

function fill_EV_Euc!(
	EV_next::AbstractArray{Float64, 4},
	c_next_CE::AbstractArray{Float64, 4},
	V_next::AbstractArray{Float64, N},
	policy_c_next::AbstractArray{Float64, N},
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
			EV  = 0.0
			Euc = 0.0
			bad = false
			for k in eachindex(Γt)
				p = Γt[k]
				if p == 0.0
					continue
				end
				Vp = Vt[k]
				c  = ct[k]
				if !isfinite(Vp) || Vp <= V_penalty || !isfinite(c) || c <= c_floor
					bad = true
					break
				end
				EV  += p * Vp
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
			EV  = 0.0
			Euc = 0.0
			bad = false
			for k in eachindex(Γt)
				p = Γt[k]
				if p == 0.0
					continue
				end
				Vp = Vt[k]
				c  = ct[k]
				if !isfinite(Vp) || Vp <= V_penalty || !isfinite(c) || c <= c_floor
					bad = true
					break
				end
				EV  += p * Vp
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
			EV  = 0.0
			Euc = 0.0
			bad = false
			for k in eachindex(Γt)
				p = Γt[k]
				if p == 0.0
					continue
				end
				Vp = Vt[k]
				c  = ct[k]
				if !isfinite(Vp) || Vp <= V_penalty || !isfinite(c) || c <= c_floor
					bad = true
					break
				end
				EV  += p * Vp
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

parameters = parameters_function()
variables = variables_function(parameters)
solve_value_and_policy_function!(variables, parameters)
save_JLD_function!(variables, parameters, filename = "workspace_benchmark_new.jld2")

# parameters_lr_exp_1 = parameters_function(lr_exp = 1)
# variables_lr_exp_1 = variables_function(parameters_lr_exp_1)
# solve_value_and_policy_function!(variables_lr_exp_1, parameters_lr_exp_1)
# save_JLD_function!(variables_lr_exp_1, parameters_lr_exp_1, filename = "workspace_lr_exp_1.jld2")

# parameters_lr_exp_2 = parameters_function(lr_exp = 2)
# variables_lr_exp_2 = variables_function(parameters_lr_exp_2)
# solve_value_and_policy_function!(variables_lr_exp_2, parameters_lr_exp_2)
# save_JLD_function!(variables_lr_exp_2, parameters_lr_exp_2, filename = "workspace_lr_exp_2.jld2")

# parameters_lr_exp_3 = parameters_function(lr_exp = 3)
# variables_lr_exp_3 = variables_function(parameters_lr_exp_3)
# solve_value_and_policy_function!(variables_lr_exp_3, parameters_lr_exp_3)
# save_JLD_function!(variables_lr_exp_3, parameters_lr_exp_3, filename = "workspace_lr_exp_3.jld2")

# parameters_lr_exp_4 = parameters_function(lr_exp = 4)
# variables_lr_exp_4 = variables_function(parameters_lr_exp_4)
# solve_value_and_policy_function!(variables_lr_exp_4, parameters_lr_exp_4)
# save_JLD_function!(variables_lr_exp_4, parameters_lr_exp_4, filename = "workspace_lr_exp_4.jld2")

# parameters_lr_exp_5 = parameters_function(lr_exp = 5)
# variables_lr_exp_5 = variables_function(parameters_lr_exp_5)
# solve_value_and_policy_function!(variables_lr_exp_5, parameters_lr_exp_5)
# save_JLD_function!(variables_lr_exp_5, parameters_lr_exp_5, filename = "workspace_lr_exp_5.jld2")

# parameters_lr_exp_6 = parameters_function(lr_exp = 6)
# variables_lr_exp_6 = variables_function(parameters_lr_exp_6)
# solve_value_and_policy_function!(variables_lr_exp_6, parameters_lr_exp_6)
# save_JLD_function!(variables_lr_exp_6, parameters_lr_exp_6, filename = "workspace_lr_exp_6.jld2")

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

# checking plots #

ϵ_ind = 1
ν_ind = 1
age_ind = 6

plot_K_a = plot(
    box=:on,
    size=[800, 600],
    # xlim=[0, 80],
    # xticks=0:10:80,
    ylim=[-0.1, 1.1],
    yticks=0:1:1,
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm,
    xlabel="Asset",
    ylabel="New Child",
	legend=:right
)

plot_K_a = plot!(
    parameters.a_grid,
    variables.policy_K[:,1,ϵ_ind,ν_ind,1,age_ind],
    label="Age $(parameters.age_grid[age_ind])",
    lw=3,
    lc=:blue
)

plot_K_a = plot!(
    parameters.a_grid,
    variables.policy_K[:,1,ϵ_ind,ν_ind,1,age_ind+1],
    label="Age $(parameters.age_grid[age_ind+1])",
    lw=3,
    lc=:red,
	ls=:dash
)

plot_K_a = plot!(
    parameters.a_grid,
    variables.policy_K[:,1,ϵ_ind,ν_ind,1,age_ind+2],
    label="Age $(parameters.age_grid[age_ind+2])",
    lw=3,
    lc=:black,
	ls=:dot	
)

savefig(plot_K_a, string("plot_K_a.pdf"))


a_ind = 10
ν_ind = 1
age_ind = 7

plot_K_ϵ = plot(
    box=:on,
    size=[800, 600],
    ylim=[-0.1, 1.1],
    yticks=0:1:1,
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm,
    xlabel="Persistent Income",
    ylabel="New Child",
	legend=:left
)

plot_K_ϵ = plot!(
    parameters.ϵ_grid,
    variables.policy_K[a_ind,1,:,ν_ind,1,age_ind],
    label="Age $(parameters.age_grid[age_ind])",
    lw=3,
    lc=:blue
)

plot_K_ϵ = plot!(
    parameters.ϵ_grid,
    variables.policy_K[a_ind,1,:,ν_ind,1,age_ind+1],
    label="Age $(parameters.age_grid[age_ind+1])",
    lw=3,
    lc=:red,
	ls=:dash
)

plot_K_ϵ = plot!(
    parameters.ϵ_grid,
    variables.policy_K[a_ind,1,:,ν_ind,1,age_ind+2],
    label="Age $(parameters.age_grid[age_ind+2])",
    lw=3,
    lc=:black,
	ls=:dot	
)

savefig(plot_K_ϵ, string("plot_K_ϵ.pdf"))


a_ind = 12
ϵ_ind = 4
age_ind = 7

plot_K_ν = plot(
    box=:on,
    size=[800, 600],
    ylim=[-0.1, 1.1],
    yticks=0:1:1,
    xtickfont=font(16, "Computer Modern", :black),
    ytickfont=font(16, "Computer Modern", :black),
    legendfont=font(16, "Computer Modern", :black),
    guidefont=font(18, "Computer Modern", :black),
    titlefont=font(18, "Computer Modern", :black),
    margin=4mm,
    xlabel="Transitory Income",
    ylabel="New Child",
	legend=:left
)

plot_K_ν = plot!(
    parameters.ν_grid,
    variables.policy_K[a_ind,1,ϵ_ind,:,1,age_ind],
    label="Age $(parameters.age_grid[age_ind])",
    lw=3,
    lc=:blue
)

plot_K_ν = plot!(
    parameters.ν_grid,
    variables.policy_K[a_ind,1,ϵ_ind,:,1,age_ind+1],
    label="Age $(parameters.age_grid[age_ind+1])",
    lw=3,
    lc=:red,
	ls=:dash
)

plot_K_ν = plot!(
    parameters.ν_grid,
    variables.policy_K[a_ind,1,ϵ_ind,:,1,age_ind+2],
    label="Age $(parameters.age_grid[age_ind+2])",
    lw=3,
    lc=:black,
	ls=:dot	
)

savefig(plot_K_ν, string("plot_K_ν.pdf"))


# simuation #

using Random
using Random123
using QuantEcon: Categorical
using Polyester

@inline function make_thread_rngs(; seed::UInt64 = UInt64(1124))
    nt = Threads.nthreads()
    key = (seed, UInt64(0))
    return [Philox4x(UInt64, key) for _ in 1:nt]
end

@inline base_counter(h_id::UInt64, t_id::UInt64) = (h_id << 44) | (t_id << 24)

@inline function ap_policy_to_index(
    rng,
    a_grid::AbstractVector{Float64},
    ap::Float64,
	a_size::Int
)::Int
    if ap ≤ a_grid[1]
        return 1
    elseif ap ≥ a_grid[end]
        return a_size
    else
        j  = searchsortedlast(a_grid, ap)
        a0 = a_grid[j]
        a1 = a_grid[j+1]
        w  = (ap - a0) / (a1 - a0)
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
    seed::UInt64 = UInt64(1124),
)
    @unpack a_grid, a_ind_zero, a_size, ϵ_size, n_size, n_Γ, ϵ_Γ, ϵ_G, ν_Γ, inf_grid, age_ret, age_inf, age_min = parameters

    T, N = size(panel.a_state)

    ϵ_init = Categorical(ϵ_G)
    ν_iid  = Categorical(ν_Γ)
    ϵ_cat  = [Categorical(ϵ_Γ[i, :]) for i in 1:ϵ_size]
    n_cat  = [Categorical(n_Γ[i, :]) for i in 1:n_size]

    rngs = make_thread_rngs(seed = seed)

    @batch for h in 1:N

        rng = rngs[Threads.threadid()]
        h_id = UInt64(h)

        set_counter!(rng, base_counter(h_id, UInt64(0)))

        panel.a_state[1, h] = a_ind_zero
        panel.n_state[1, h] = 1
        panel.f_state[1, h] = 1
        panel.ϵ_state[1, h] = rand(rng, ϵ_init)
        panel.ν_state[1, h] = rand(rng, ν_iid)

        @inbounds for t in 1:(T-1)

			set_counter!(rng, base_counter(h_id, UInt64(t)))

            a_i = panel.a_state[t, h]
            n_i = panel.n_state[t, h]
            ϵ_i = panel.ϵ_state[t, h]
            ν_i = panel.ν_state[t, h]
            f_i = panel.f_state[t, h]

            c_raw  = variables.policy_c[a_i, n_i, ϵ_i, ν_i, f_i, t]
            e_raw  = variables.policy_e[a_i, n_i, ϵ_i, ν_i, f_i, t]
            panel.c_choice[t, h] = c_raw
            panel.e_choice[t, h] = e_raw

			K_raw  = variables.policy_K[a_i, n_i, ϵ_i, ν_i, f_i, t]
			K01 = (K_raw >= 0.5) ? Int8(1) : Int8(0)
            panel.K_choice[t, h] = K01

			Δa_i = a_i != a_size ? a_i + 1 : a_size
            ΔK_raw  = variables.policy_K[Δa_i, n_i, ϵ_i, ν_i, f_i, t]
			ΔK01 = (ΔK_raw >= 0.5) ? Int8(1) : Int8(0)
            panel.ΔK_choice[t, h] = ΔK01

			ap_raw = variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, f_i, t]
            ap_idx = ap_policy_to_index(rng, a_grid, ap_raw, a_size)
            panel.ap_choice[t, h] = ap_idx
            panel.a_state[t+1, h] = ap_idx

			if t < age_ret - age_min
				n_eff = min(n_size, n_i + Int(K_raw))
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

T = parameters.age_size
N = 300_000
panel = initialize_panel(num_households=N, num_periods=T)
simulate_household_panel!(panel, variables, parameters; seed=UInt64(20260118))

# panel_a_low_inf, panel_a_p_low_inf, panel_x_low_inf, panel_l_low_inf, panel_n_low_inf, panel_K_low_inf, shock_ϵ_low_inf, shock_ν_low_inf, shock_f_low_inf = simulation_function(num_hh = num_hh, filename = "workspace_low_inf.jld2")

# panel_a_no_inf, panel_a_p_no_inf, panel_x_no_inf, panel_l_no_inf, panel_n_no_inf, panel_K_no_inf, shock_ϵ_no_inf, shock_ν_no_inf, shock_f_no_inf = simulation_function(num_hh = num_hh, filename = "workspace_no_inf.jld2")

# panel_a_edu_h, panel_a_p_edu_h, panel_x_edu_h, panel_l_edu_h, panel_n_edu_h, panel_K_edu_h, shock_ϵ_edu_h, shock_ν_edu_h, shock_f_edu_h = simulation_function(num_hh = num_hh, filename = "workspace_edu_h.jld2")

# panel_a_edu_h_low_σ, panel_a_p_edu_h_low_σ, panel_x_edu_h_low_σ, panel_l_edu_h_low_σ, panel_n_edu_h_low_σ, panel_K_edu_h_low_σ, shock_ϵ_edu_h_low_σ, shock_ν_edu_h_low_σ, shock_f_edu_h_low_σ = simulation_function(num_hh = num_hh, filename = "workspace_edu_h_low_σ.jld2")


# simulation results #

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


# Short-Run Response #

using DataFrames
using CategoricalArrays
using FixedEffectModels
using Plots

function _term_group_index(term::AbstractString, levels::Vector{String})
    for (k, lev) in pairs(levels)
        occursin(lev, term) && return k
    end
    return missing
end

function short_run_regression(panel::SimulatedPanel, parameters; savepath::String="plot_short_run.pdf")

    @unpack a_size, a_grid, age_grid = parameters

    panel_K  = permutedims(Float64.(panel.K_choice))
    panel_dK = permutedims(Float64.(panel.ΔK_choice)) .- panel_K
    panel_a  = permutedims(panel.a_state)
    panel_a_adj = clamp.(panel_a, 1, a_size-1)
    panel_da = a_grid[panel_a_adj .+ 1] .- a_grid[panel_a_adj]

    I, J = size(panel_dK)

    df = DataFrame(
        i   = repeat(1:I, J),
        age = repeat(age_grid, inner=I),
        dK  = vec(panel_dK),
        ai  = vec(panel_a),
        da  = vec(panel_da),
    )

    df = df[df.ai .!= a_size, :]

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
        elseif 45 <= j <= 49
            "45-49"
        else
            missing
        end
    end
    dropmissing!(df, [:age_group])

    levels = ["15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"]
    df.age_group = CategoricalArray(df.age_group; ordered=true, levels=levels)

    m = reg(df, @formula(dK ~ 0 + age_group + da & age_group), Vcov.cluster(:i))
    println(m)

    out = DataFrame(term = coefnames(m), b = coef(m), se = stderror(m))

    alpha = out[contains.(out.term, "da") .&& contains.(out.term, "age_group"), :]
    alpha.grp = [_term_group_index(t, levels) for t in alpha.term]
    dropmissing!(alpha, [:grp])
    sort!(alpha, :grp)

    # plot
    x = 1:nrow(alpha)
    y = alpha.b
    yerr = 1.96 .* alpha.se
    plot_sr = plot(
        x, y;
        yerror = yerr,
        xlabel = "Age group",
        ylabel = "α_J",
        xticks = (x, levels[alpha.grp]),
        legend = false,
    )
    savefig(plot_sr, savepath)

    return (m=m, out=out, alpha=alpha, df=df, plot=plot_sr)
end


using DataFrames
using CategoricalArrays
using FixedEffectModels

I, J = size(panel_ΔK)
panel_dK = panel_ΔK .- panel_K
panel_a_adj = copy(panel_a)
panel_a_adj[panel_a .== parameters.a_size] .= parameters.a_size - 1
panel_da = parameters.a_grid[panel_a_adj .+ 1] .- parameters.a_grid[panel_a_adj]
df = DataFrame(
	i = repeat(1:I, J),
	age = repeat(parameters.age_grid, inner = I),
	dK = vec(panel_dK),
	ai = vec(panel_a),
	da = vec(panel_da),
)
df = df[df.ai .!= parameters.a_size, :]
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
df.age_group = CategoricalArray(df.age_group, ordered = true)

# run the age-group fixed effects model
m = reg(df, @formula(dK ~ 0 + age_group + da & age_group), Vcov.cluster(:i))
# m = reg(df, @formula(dK ~ 0 + da & age_group ), Vcov.cluster(:i))
println(m)
out = DataFrame(term = coefnames(m), b = coef(m), se = stderror(m))
alpha = out[contains.(out.term, "da").&&contains.(out.term, "age_group"), :]

# plot 
x = 1:nrow(alpha)
y = alpha.b
yerr = 1.96 .* alpha.se
age_group = ["15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-50"]
plot_sr = plot(
	x, y,
	# seriestype = :scatter,
	yerror = yerr,
	xlabel = "Age group",
	ylabel = "α_J",
	xticks = (x, age_group),
	legend = false,
	# xrotation = 45,
)
savefig(plot_sr, string("plot_short_run.pdf"))#===================##===================#


# Long-Run Response #

plot_lr_h = plot(parameters.age_grid, parameters.h_grid, labels = "Benchmark", legend = :bottomright)
plot_lr_h = plot!(parameters_lr_exp_1.age_grid, parameters_lr_exp_1.h_grid, labels = "20-24")
plot_lr_h = plot!(parameters_lr_exp_2.age_grid, parameters_lr_exp_2.h_grid, labels = "25-29")
plot_lr_h = plot!(parameters_lr_exp_3.age_grid, parameters_lr_exp_3.h_grid, labels = "30-34")
plot_lr_h = plot!(parameters_lr_exp_4.age_grid, parameters_lr_exp_4.h_grid, labels = "35-39")
plot_lr_h = plot!(parameters_lr_exp_5.age_grid, parameters_lr_exp_5.h_grid, labels = "40-44")
plot_lr_h = plot!(parameters_lr_exp_6.age_grid, parameters_lr_exp_6.h_grid, labels = "45-49")
savefig(plot_lr_h, string("plot_long_run_h.pdf"))

using DataFrames
using CategoricalArrays
using FixedEffectModels

panel_n_max = maximum(parameters.n_grid[panel_n], dims = 2)
panel_n_1_max = maximum(parameters.n_grid[panel_n_lr_exp_1], dims = 2)
panel_n_2_max = maximum(parameters.n_grid[panel_n_lr_exp_2], dims = 2)
panel_n_3_max = maximum(parameters.n_grid[panel_n_lr_exp_3], dims = 2)
panel_n_4_max = maximum(parameters.n_grid[panel_n_lr_exp_4], dims = 2)
panel_n_5_max = maximum(parameters.n_grid[panel_n_lr_exp_5], dims = 2)
panel_n_6_max = maximum(parameters.n_grid[panel_n_lr_exp_6], dims = 2)
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
	i = repeat(1:I, 6),
	exp = repeat(age_group_lr, inner = I),
	dN = vec(panel_dN),
)

# run the long-run model
m_lr = reg(df_lr, @formula(dN ~ 0 + exp), Vcov.cluster(:i))
println(m_lr)
beta = DataFrame(term = coefnames(m_lr), b = coef(m_lr), se = stderror(m_lr))

# plot 
x = 1:nrow(beta)
y = beta.b .* 1000
yerr = 1.96 .* beta.se .* 1000
plot_lr = plot(
	x, y,
	# seriestype = :scatter,
	yerror = yerr,
	xlabel = "Age group",
	ylabel = "β_J",
	xticks = (x, age_group_lr),
	legend = false,
	# xrotation = 45,
)
savefig(plot_lr, string("plot_long_run.pdf"))
