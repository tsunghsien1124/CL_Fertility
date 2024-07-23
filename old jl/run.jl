#=================#
# Import packages #
#=================#
using Dierckx
using FLOWMath
using Distributions
using QuadGK
using JLD2: @save, @load
using LinearAlgebra
using Optim
using Parameters: @unpack
using PrettyTables
using ProgressMeter
using QuantEcon: gridmake, rouwenhorst, tauchen, stationary_distributions, MarkovChain
using Roots
using CSV
using Tables
using Plots
using Random
using GLM
using DataFrames
using Measures

#==================#
# Import functions #
#==================#
include("solving_stationary_equilibrium.jl")
# include("solving_transitional_dynamics.jl")
# include("simulation.jl")

#==============================#
# Solve stationary equilibrium #
#==============================#
parameters = parameters_function()
variables = variables_function(parameters; λ=0.0)
slow_updating = 1.0

# crit_V = solve_value_and_pricing_function!(variables, parameters; tol = tol_h, iter_max = 500, slow_updating = slow_updating)
# ED_KL_to_D_ratio_min, ED_leverage_ratio_min, crit_V_min, crit_μ_min = solve_economy_function!(variables, parameters; slow_updating=slow_updating)

# e_1_i = 2
# plot(parameters.a_grid_neg,variables.q[1:parameters.a_ind_zero,e_1_i,:],legend=:none)
# plot(parameters.a_grid_neg,-parameters.a_grid_neg.*variables.q[1:parameters.a_ind_zero,e_1_i,:],legend=:none)
# plot!(variables.rbl[e_1_i,:,1],-variables.rbl[e_1_i,:,2],seriestype=:scatter,legend=:none)
# vline!(variables.rbl[e_1_i,:,1],color=:black)

# plot(parameters.a_grid_neg,variables.V[1:parameters.a_ind_zero,e_1_i,:,2,1])
# plot(parameters.a_grid_neg,variables.V_nd[1:parameters.a_ind_zero,e_1_i,:,2,1])
# plot(parameters.e_2_grid, variables.threshold_a[2,:,:,3])
# plot(parameters.a_grid_neg, variables.threshold_e_2[:,1,:,2])

parameters = parameters_function()

variables_min = variables_function(parameters; λ=0.0)
ED_KL_to_D_ratio_min, ED_leverage_ratio_min, crit_V_min, crit_μ_min = solve_economy_function!(variables_min, parameters; slow_updating=slow_updating)

variables_max = variables_function(parameters; λ=1.0 - sqrt(parameters.ψ))
ED_KL_to_D_ratio_max, ED_leverage_ratio_max, crit_V_max, crit_μ_max = solve_economy_function!(variables_max, parameters; slow_updating=slow_updating)

data_spec = Any[
    "λ" variables_min.aggregate_prices.λ variables_max.aggregate_prices.λ
    "(K'+L')/D' (supply)" variables_min.aggregate_prices.KL_to_D_ratio_λ variables_max.aggregate_prices.KL_to_D_ratio_λ
    "(K'+L')/D' (demand)" variables_min.aggregate_variables.KL_to_D_ratio variables_max.aggregate_variables.KL_to_D_ratio
    "(K'+L')/N (supply)" variables_min.aggregate_prices.leverage_ratio_λ variables_max.aggregate_prices.leverage_ratio_λ
    "(K'+L')/N (demand)" variables_min.aggregate_variables.leverage_ratio variables_max.aggregate_variables.leverage_ratio
    "share of filers (%)" variables_min.aggregate_variables.share_of_filers*100 variables_max.aggregate_variables.share_of_filers*100
    "share in debts (%)" variables_min.aggregate_variables.share_in_debts*100 variables_max.aggregate_variables.share_in_debts*100
    "debt-to-earnings ratio (%)" variables_min.aggregate_variables.debt_to_earning_ratio*100 variables_max.aggregate_variables.debt_to_earning_ratio*100
    "avg interest rate (%)" variables_min.aggregate_variables.avg_loan_rate*100 variables_max.aggregate_variables.avg_loan_rate*100
    "policy upper bound" variables_min.policy_a[end, end, end, end, 2]<parameters.a_grid[end] variables_max.policy_a[end, end, end, end, 2]<parameters.a_grid[end]
    "convergence of V" crit_V_min crit_V_max
    "convergence of μ" crit_μ_min crit_μ_max
]
pretty_table(data_spec; header=["Name", "λ minimum", "λ maximum"], alignment=[:c, :c, :c], formatters=ft_round(8), body_hlines=[3, 5, 9, 10])
