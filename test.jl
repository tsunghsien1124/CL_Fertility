using JLD2
using Plots

old = load("workspace_benchmark.jld2")
new = load("workspace_benchmark_new.jld2")

age_ret_ind = old["parameters"].age_max - old["parameters"].age_ret + 1
@assert all(old["V"][:,1,:,:,2,end] .≈ new["V"][:,1,:,:,2,end])

ϵ_ind = 4
ν_ind = 4
plot(old["parameters"].a_grid, old["parameters"].a_grid[old["policy_a_p"][:,1,ϵ_ind,ν_ind,2,end-1]])
plot!(new["parameters"].a_grid, new["policy_a_p"][:,1,ϵ_ind,ν_ind,2,end-1])