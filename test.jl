using Random
using BenchmarkTools

a_size = 20
b_size = 30
c_size = 40
d_size = 50
e_size = 60
a_grid = rand(a_size)
b_grid = rand(b_size)
c_grid = rand(c_size)
d_grid = rand(d_size)
e_grid = rand(e_size)

function sol_reshape(a_size, a_grid, b_size, b_grid, c_size, c_grid, d_size, d_grid, e_size, e_grid)
    res_reshape = Array{Float64}(undef, a_size, b_size, c_size, d_size, e_size)
    a = reshape(a_grid, (a_size, 1, 1, 1, 1))
    b = reshape(b_grid, (1, b_size, 1, 1, 1))
    c = reshape(c_grid, (1, 1, c_size, 1, 1))
    d = reshape(d_grid, (1, 1, 1, d_size, 1))
    e = reshape(e_grid, (1, 1, 1, 1, e_size))
    res_reshape .= a .+ b .- c .* d .- e
    return res_reshape
end

function sol_loop(a_size, a_grid, b_size, b_grid, c_size, c_grid, d_size, d_grid, e_size, e_grid)
    res_loop = Array{Float64}(undef, a_size, b_size, c_size, d_size, e_size)
    # Threads.@threads for (a_i, b_i, c_i) in collect(Iterators.product(1:a_size, 1:b_size, 1:c_size))
    # for a_i = 1:a_size, b_i = 1:b_size, c_i = 1:c_size, d_i = 1:d_size, e_i = 1:e_size
    # Threads.@threads for (e_i, d_i, c_i, b_i, a_i) in collect(Iterators.product(1:e_size, 1:d_size, 1:c_size, 1:b_size, 1:a_size))
    Threads.@threads for ijklm in CartesianIndices((e_size, d_size, c_size, b_size, a_size))
    # for e_i = 1:e_size, d_i = 1:d_size, c_i = 1:c_size, b_i = 1:b_size, a_i = 1:a_size
        e_i, d_i, c_i, b_i, a_i = Tuple(ijklm)
        @inbounds res_loop[a_i,b_i,c_i,d_i,e_i] = a_grid[a_i] + b_grid[b_i] - c_grid[c_i] * d_grid[d_i] - e_grid[e_i]
    end
    return res_loop
end

function sol_loop_inline(a_size, a_grid, b_size, b_grid, c_size, c_grid, d_size, d_grid, e_size, e_grid)
    # res_1oop_inline = Array{Float64}(undef, a_size, b_size, c_size, d_size, e_size)
    # res_1oop_inline .= [a_grid[a_i] + b_grid[b_i] - c_grid[c_i] * d_grid[d_i] - e_grid[e_i] for a_i = 1:a_size, b_i = 1:b_size, c_i = 1:c_size, d_i = 1:d_size, e_i = 1:e_size]
    return [a_grid[a_i] + b_grid[b_i] - c_grid[c_i] * d_grid[d_i] - e_grid[e_i] for a_i = 1:a_size, b_i = 1:b_size, c_i = 1:c_size, d_i = 1:d_size, e_i = 1:e_size]
    # return res_1oop_inline
end

res_reshape = sol_reshape(a_size, a_grid, b_size, b_grid, c_size, c_grid, d_size, d_grid, e_size, e_grid);
res_loop = sol_loop(a_size, a_grid, b_size, b_grid, c_size, c_grid, d_size, d_grid, e_size, e_grid);
res_1oop_inline = sol_loop_inline(a_size, a_grid, b_size, b_grid, c_size, c_grid, d_size, d_grid, e_size, e_grid);

@assert res_reshape == res_loop == res_1oop_inline

@btime sol_reshape($a_size, $a_grid, $b_size, $b_grid, $c_size, $c_grid, $d_size, $d_grid, $e_size, $e_grid);
@btime sol_loop($a_size, $a_grid, $b_size, $b_grid, $c_size, $c_grid, $d_size, $d_grid, $e_size, $e_grid);
@btime sol_loop_inline($a_size, $a_grid, $b_size, $b_grid, $c_size, $c_grid, $d_size, $d_grid, $e_size, $e_grid);

# @benchmark sol_reshape($a_size, $a_grid, $b_size, $b_grid, $c_size, $c_grid, $d_size, $d_grid, $e_size, $e_grid)

# @benchmark sol_loop($a_size, $a_grid, $b_size, $b_grid, $c_size, $c_grid, $d_size, $d_grid, $e_size, $e_grid)

# @benchmark sol_loop_inline($a_size, $a_grid, $b_size, $b_grid, $c_size, $c_grid, $d_size, $d_grid, $e_size, $e_grid)

# # # using BenchmarkTools

# # # parameters = parameters_function()

# # # function sol_1!(variables::Mutable_Variables, parameters::NamedTuple)

# # #     # unpack parameters
# # #     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
# # #     @unpack ν_size, ν_grid, ν_Γ = parameters
# # #     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
# # #     @unpack n_size, n_grid, n_Γ, n_max = parameters
# # #     @unpack a_size, a_grid = parameters
# # #     @unpack l_size, l_grid, x_size, x_grid = parameters
# # #     @unpack inf_size, inf_grid = parameters
# # #     @unpack h_grid = parameters
# # #     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

# # #     age_i = age_size
# # #     age = age_grid[age_i]
# # #     age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
# # #     age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
# # #     for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, a_i = 1:a_size
# # #         ν = ν_grid[ν_i]
# # #         ϵ = ϵ_grid[ϵ_i]
# # #         a = a_grid[a_i]
# # #         w_bar = exp(h + ϵ + ν) * b
# # #         @inbounds variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = utility_function((1.0 + r) * a + w_bar, 0.0, 0.0, γ, ψ, κ, q_bar)
# # #     end
# # #     nothing
# # # end

# # # function sol_2!(variables::Mutable_Variables, parameters::NamedTuple)

# # #     # unpack parameters
# # #     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
# # #     @unpack ν_size, ν_grid, ν_Γ = parameters
# # #     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
# # #     @unpack n_size, n_grid, n_Γ, n_max = parameters
# # #     @unpack a_size, a_grid = parameters
# # #     @unpack l_size, l_grid, x_size, x_grid = parameters
# # #     @unpack inf_size, inf_grid = parameters
# # #     @unpack h_grid = parameters
# # #     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

# # #     age_i = age_size
# # #     age = age_grid[age_i]
# # #     age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
# # #     age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
# # #     Threads.@threads for (ν_i, ϵ_i, a_i) in collect(Iterators.product(1:ν_size, 1:ϵ_size, 1:a_size))
# # #         ν = ν_grid[ν_i]
# # #         ϵ = ϵ_grid[ϵ_i]
# # #         a = a_grid[a_i]
# # #         w_bar = exp(h + ϵ + ν) * b
# # #         @inbounds variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = utility_function((1.0 + r) * a + w_bar, 0.0, 0.0, γ, ψ, κ, q_bar)
# # #     end
# # #     nothing
# # # end

# # # function sol_3!(variables::Mutable_Variables, parameters::NamedTuple)

# # #     # unpack parameters
# # #     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
# # #     @unpack ν_size, ν_grid, ν_Γ = parameters
# # #     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
# # #     @unpack n_size, n_grid, n_Γ, n_max = parameters
# # #     @unpack a_size, a_grid = parameters
# # #     @unpack l_size, l_grid, x_size, x_grid = parameters
# # #     @unpack inf_size, inf_grid = parameters
# # #     @unpack h_grid = parameters
# # #     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

# # #     age_i = age_size
# # #     age = age_grid[age_i]
# # #     age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
# # #     age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
# # #     res = [utility_function((1.0 + r) * a_grid[a_i] + exp(h + ϵ_grid[ϵ_i] + ν_grid[ν_i]) * b, 0.0, 0.0, γ, ψ, κ, q_bar) for a_i = 1:a_size, ϵ_i = 1:ϵ_size, ν_i = 1:ν_size];
# # #     @views @inbounds variables.V[:, 1, :, :, 2, age_i] .= res
# # #     nothing
# # # end

# # # function sol_4!(variables::Mutable_Variables, parameters::NamedTuple)

# # #     # unpack parameters
# # #     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
# # #     @unpack ν_size, ν_grid, ν_Γ = parameters
# # #     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
# # #     @unpack n_size, n_grid, n_Γ, n_max = parameters
# # #     @unpack a_size, a_grid = parameters
# # #     @unpack l_size, l_grid, x_size, x_grid = parameters
# # #     @unpack inf_size, inf_grid = parameters
# # #     @unpack h_grid = parameters
# # #     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

# # #     age_i = age_size
# # #     age = age_grid[age_i]
# # #     age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
# # #     age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
# # #     a = reshape(a_grid, (a_size, 1, 1))
# # #     ϵ = reshape(ϵ_grid, (1, ϵ_size, 1))
# # #     ν = reshape(ν_grid, (1, 1, ν_size))
# # #     c = (1.0 + r) .* a .+ exp.(h .+ ϵ .+ ν) .* b
# # #     @inbounds variables.V[:, 1, :, :, 2, age_i] .= utility_function.(c, Ref(0.0), Ref(0.0), Ref(γ), Ref(ψ), Ref(κ), Ref(q_bar))
# # #     nothing
# # # end

# # # variables_1 = variables_function(parameters)
# # # variables_2 = variables_function(parameters)
# # # variables_3 = variables_function(parameters)
# # # variables_4 = variables_function(parameters)

# # # sol_1!(variables_1, parameters)
# # # sol_2!(variables_2, parameters)
# # # sol_3!(variables_3, parameters)
# # # sol_4!(variables_4, parameters)

# # # @assert variables_1.V[:, 1, :, :, 2, end] == variables_2.V[:, 1, :, :, 2, end] == variables_3.V[:, 1, :, :, 2, end] == variables_4.V[:, 1, :, :, 2, end]

# # # @btime sol_1!($variables_1, $parameters)
# # # @btime sol_2!($variables_2, $parameters)
# # # @btime sol_3!($variables_3, $parameters)
# # # @btime sol_4!($variables_4, $parameters)

# # # @benchmark sol_4!($variables_4, $parameters)

# # using BenchmarkTools
# # using LinearAlgebra

# # parameters = parameters_function()

# # function sol_1!(variables::Mutable_Variables, parameters::NamedTuple)

# #     # unpack parameters
# #     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
# #     @unpack ν_size, ν_grid, ν_Γ = parameters
# #     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
# #     @unpack n_size, n_grid, n_Γ, n_max = parameters
# #     @unpack a_size, a_grid = parameters
# #     @unpack l_size, l_grid, x_size, x_grid = parameters
# #     @unpack inf_size, inf_grid = parameters
# #     @unpack h_grid = parameters
# #     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

# #     for age_i = age_size:(-1):1
# #         age = age_grid[age_i]
# #         age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
# #         age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
# #         if age == age_max # terminal condition
# #             for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, a_i = 1:a_size
# #                 ν = ν_grid[ν_i]
# #                 ϵ = ϵ_grid[ϵ_i]
# #                 a = a_grid[a_i]
# #                 w_bar = exp(h + ϵ + ν) * b
# #                 @inbounds variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = utility_function((1.0 + r) * a + w_bar, 0.0, 0.0, γ, ψ, κ, q_bar)
# #             end
# #         elseif age_ret < age < age_max # after retirement
# #             for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, a_i = 1:a_size
# #                 ν = ν_grid[ν_i]
# #                 ϵ = ϵ_grid[ϵ_i]
# #                 a = a_grid[a_i]
# #                 w_bar = exp(h + ϵ + ν) * b
# #                 EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
# #                 V_all = utility_function.((1.0 + r) * a + w_bar .- a_grid, 0.0, 0.0, γ, ψ, κ, q_bar) .+ β * EV
# #                 V_max_i = argmax(V_all)
# #                 @inbounds variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
# #                 @inbounds variables.policy_a_p[a_i, 1, ϵ_i, ν_i, 2, age_i] = V_max_i # a_grid[V_max_i]
# #             end
# #         end
# #     end
# #     nothing
# # end

# # function sol_2!(variables::Mutable_Variables, parameters::NamedTuple)

# #     # unpack parameters
# #     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
# #     @unpack ν_size, ν_grid, ν_Γ = parameters
# #     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
# #     @unpack n_size, n_grid, n_Γ, n_max = parameters
# #     @unpack a_size, a_grid = parameters
# #     @unpack l_size, l_grid, x_size, x_grid = parameters
# #     @unpack inf_size, inf_grid = parameters
# #     @unpack h_grid = parameters
# #     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

# #     for age_i = age_size:(-1):1
# #         age = age_grid[age_i]
# #         age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
# #         age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
# #         if age == age_max # terminal condition
# #             a = reshape(a_grid, (a_size, 1, 1))
# #             ϵ = reshape(ϵ_grid, (1, ϵ_size, 1))
# #             ν = reshape(ν_grid, (1, 1, ν_size))
# #             c = (1.0 + r) .* a .+ exp.(h .+ ϵ .+ ν) .* b
# #             @inbounds variables.V[:, 1, :, :, 2, age_i] .= utility_function.(c, Ref(0.0), Ref(0.0), Ref(γ), Ref(ψ), Ref(κ), Ref(q_bar))
# #         elseif age_ret < age < age_max # after retirement
# #             for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, a_i = 1:a_size
# #                 ν = ν_grid[ν_i]
# #                 ϵ = ϵ_grid[ϵ_i]
# #                 a = a_grid[a_i]
# #                 w_bar = exp(h + ϵ + ν) * b
# #                 EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
# #                 V_all = utility_function.((1.0 + r) * a + w_bar .- a_grid, 0.0, 0.0, γ, ψ, κ, q_bar) .+ β * EV
# #                 V_max_i = argmax(V_all)
# #                 @inbounds variables.V[a_i, 1, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
# #                 @inbounds variables.policy_a_p[a_i, 1, ϵ_i, ν_i, 2, age_i] = V_max_i # a_grid[V_max_i]
# #             end
# #         end
# #     end
# #     nothing
# # end

# # function sol_3!(variables::Mutable_Variables, parameters::NamedTuple)

# #     # unpack parameters
# #     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
# #     @unpack ν_size, ν_grid, ν_Γ = parameters
# #     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
# #     @unpack n_size, n_grid, n_Γ, n_max = parameters
# #     @unpack a_size, a_grid = parameters
# #     @unpack l_size, l_grid, x_size, x_grid = parameters
# #     @unpack inf_size, inf_grid = parameters
# #     @unpack h_grid = parameters
# #     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

# #     for age_i = age_size:(-1):1
# #         age = age_grid[age_i]
# #         age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
# #         age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
# #         if age == age_max # terminal condition
# #             a = reshape(a_grid, (a_size, 1, 1))
# #             ϵ = reshape(ϵ_grid, (1, ϵ_size, 1))
# #             ν = reshape(ν_grid, (1, 1, ν_size))
# #             c = (1.0 + r) .* a .+ exp.(h .+ ϵ .+ ν) .* b
# #             @inbounds variables.V[:, 1, :, :, 2, age_i] .= utility_function.(c, Ref(0.0), Ref(0.0), Ref(γ), Ref(ψ), Ref(κ), Ref(q_bar))
# #         elseif age_ret < age < age_max # after retirement
# #             a = reshape(a_grid, (a_size, 1, 1, 1))
# #             ϵ = reshape(ϵ_grid, (1, ϵ_size, 1, 1))
# #             ν = reshape(ν_grid, (1, 1, ν_size, 1))
# #             ap = reshape(a_grid, (1, 1, 1, a_size))
# #             c = (1.0 + r) .* a .+ exp.(h .+ ϵ .+ ν) .* b .- ap
# #             ϵ_I = reshape(Matrix{Float64}(I, ϵ_size, ϵ_size), (1, ϵ_size, 1, 1, ϵ_size, 1))
# #             ν_I = reshape(Matrix{Float64}(I, ν_size, ν_size), (1, 1, ν_size, 1, 1, ν_size))
# #             VP = reshape(variables.V[:, 1, :, :, 2, age_i+1], (1, 1, 1, a_size, ϵ_size, ν_size))
# #             EV = dropdims(reduce(+, VP .* ϵ_I .* ν_I, dims=5:6),dims=(5,6))
# #             V_all = utility_function.(c, Ref(0.0), Ref(0.0), Ref(γ), Ref(ψ), Ref(κ), Ref(q_bar)) .+ β .* EV
# #             @inbounds variables.V[:, 1, :, :, 2, age_i] .= maximum(V_all, dims=4)
# #             @inbounds variables.policy_a_p[:, 1, :, :, 2, age_i] .= dropdims(mapslices(argmax,V_all,dims=4),dims=4)
# #         end
# #     end
# #     nothing
# # end

# # variables_1 = variables_function(parameters)
# # variables_2 = variables_function(parameters)
# # variables_3 = variables_function(parameters)

# # sol_1!(variables_1, parameters)
# # sol_2!(variables_2, parameters)
# # sol_3!(variables_3, parameters)

# # @assert variables_1.V == variables_2.V == variables_3.V
# # @assert variables_1.policy_a_p == variables_2.policy_a_p == variables_3.policy_a_p

# # @btime sol_1!($variables_1, $parameters)
# # @btime sol_2!($variables_2, $parameters)
# # @btime sol_3!($variables_3, $parameters)

# using BenchmarkTools
# using LinearAlgebra

# parameters = parameters_function()

# function sol_1!(variables::Mutable_Variables, parameters::NamedTuple)

#     # unpack parameters
#     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
#     @unpack ν_size, ν_grid, ν_Γ = parameters
#     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
#     @unpack n_size, n_grid, n_Γ, n_max = parameters
#     @unpack a_size, a_grid = parameters
#     @unpack l_size, l_grid, x_size, x_grid = parameters
#     @unpack inf_size, inf_grid = parameters
#     @unpack h_grid = parameters
#     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

#     for age_i = age_size:(-1):1
#         age = age_grid[age_i]
#         age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
#         age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
#         if age == age_ret # at retirement age
#             for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
#                 ν = ν_grid[ν_i]
#                 ϵ = ϵ_grid[ϵ_i]
#                 n = n_grid[n_i]
#                 a = a_grid[a_i]
#                 w = exp(h + ϵ + ν)
#                 if n == 0
#                     EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
#                     V_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV
#                     V_max_i = argmax(V_all)
#                     @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
#                     @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_max_i # a_grid[V_max_i]
#                 else
#                     EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
#                     Sol_all = zeros(a_size * x_size * l_size, 4)
#                     Sol_all_i = 0
#                     for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
#                         Sol_all_i += 1
#                         @inbounds Sol_all[Sol_all_i, 1] = a_p_i
#                         @inbounds Sol_all[Sol_all_i, 2] = x_i
#                         @inbounds Sol_all[Sol_all_i, 3] = l_i
#                         a_p = a_grid[a_p_i]
#                         x = x_grid[x_i]
#                         l = l_grid[l_i]
#                         q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
#                         @inbounds Sol_all[Sol_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV[a_p_i]
#                     end
#                     V_max_i = argmax(Sol_all[:, 4])
#                     @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Sol_all[V_max_i, 4]
#                     @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 1]) # a_grid[Int(Sol_all[V_max_i, 1])]
#                     @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 2]) # x_grid[Int(Sol_all[V_max_i, 2])]
#                     @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 3]) # l_grid[Int(Sol_all[V_max_i, 3])]
#                 end
#             end
#         end
#     end
#     nothing
# end

# function sol_2!(variables::Mutable_Variables, parameters::NamedTuple)

#     # unpack parameters
#     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
#     @unpack ν_size, ν_grid, ν_Γ = parameters
#     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
#     @unpack n_size, n_grid, n_Γ, n_max = parameters
#     @unpack a_size, a_grid = parameters
#     @unpack l_size, l_grid, x_size, x_grid = parameters
#     @unpack inf_size, inf_grid = parameters
#     @unpack h_grid = parameters
#     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

#     for age_i = age_size:(-1):1
#         age = age_grid[age_i]
#         age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
#         age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
#         if age == age_ret # at retirement age
#             Threads.@threads for (ν_i, ϵ_i, n_i, a_i) in collect(Iterators.product(1:ν_size, 1:ϵ_size, 1:n_size, 1:a_size))
#                 ν = ν_grid[ν_i]
#                 ϵ = ϵ_grid[ϵ_i]
#                 n = n_grid[n_i]
#                 a = a_grid[a_i]
#                 w = exp(h + ϵ + ν)
#                 if n == 0
#                     EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
#                     V_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV
#                     V_max_i = argmax(V_all)
#                     @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
#                     @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_max_i # a_grid[V_max_i]
#                 else
#                     EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
#                     Sol_all = zeros(a_size * x_size * l_size, 4)
#                     Sol_all_i = 0
#                     for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
#                         Sol_all_i += 1
#                         @inbounds Sol_all[Sol_all_i, 1] = a_p_i
#                         @inbounds Sol_all[Sol_all_i, 2] = x_i
#                         @inbounds Sol_all[Sol_all_i, 3] = l_i
#                         a_p = a_grid[a_p_i]
#                         x = x_grid[x_i]
#                         l = l_grid[l_i]
#                         q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
#                         @inbounds Sol_all[Sol_all_i, 4] = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * EV[a_p_i]
#                     end
#                     V_max_i = argmax(Sol_all[:, 4])
#                     @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Sol_all[V_max_i, 4]
#                     @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 1]) # a_grid[Int(Sol_all[V_max_i, 1])]
#                     @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 2]) # x_grid[Int(Sol_all[V_max_i, 2])]
#                     @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = Int(Sol_all[V_max_i, 3]) # l_grid[Int(Sol_all[V_max_i, 3])]
#                 end
#             end
#         end
#     end
#     nothing
# end

# function sol_3!(variables::Mutable_Variables, parameters::NamedTuple)

#     # unpack parameters
#     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
#     @unpack ν_size, ν_grid, ν_Γ = parameters
#     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
#     @unpack n_size, n_grid, n_Γ, n_max = parameters
#     @unpack a_size, a_grid = parameters
#     @unpack l_size, l_grid, x_size, x_grid = parameters
#     @unpack inf_size, inf_grid = parameters
#     @unpack h_grid = parameters
#     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

#     for age_i = age_size:(-1):1
#         age = age_grid[age_i]
#         age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
#         age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
#         if age == age_ret # at retirement age
#             for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
#                 ν = ν_grid[ν_i]
#                 ϵ = ϵ_grid[ϵ_i]
#                 n = n_grid[n_i]
#                 a = a_grid[a_i]
#                 w = exp(h + ϵ + ν)
#                 if n == 0
#                     EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
#                     V_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV
#                     V_max_i = argmax(V_all)
#                     @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
#                     @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_max_i # a_grid[V_max_i]
#                 else
#                     a_p = reshape(a_grid, (a_size, 1, 1))
#                     x = reshape(x_grid, (1, x_size, 1))
#                     l = reshape(l_grid, (1, 1, l_size))
#                     q = quality_function.(x, l, Ref(n), Ref(μ), Ref(θ), Ref(ψ_1), Ref(ψ_2))
#                     c = (1.0 + r) * a .+ (1.0 .- l) * w .- a_p .- q_x .* x
#                     # V_all = [c[i,j,k] > 0 && q[1,j,k] >= q_bar ? c[i,j,k]^(1.0 - γ) / (1.0 - γ) + ψ * (n * q[1,j,k])^(1.0 - κ)/ (1.0 - κ) + β * variables.V[i, 1, ϵ_i, ν_i, 2, age_i+1] : -10^16 for i = 1:a_size, j = 1:x_size, k = 1:l_size]
#                     V_all = [c[i,j,k] > 0 && q[1,j,k] >= q_bar ? 1.0 / ((1.0 - γ) * c[i,j,k]^(γ - 1.0)) + ψ / ((1.0 - κ) * (n * q[1,j,k])^(κ - 1.0)) + β * variables.V[i, 1, ϵ_i, ν_i, 2, age_i+1] : -10^16 for i = 1:a_size, j = 1:x_size, k = 1:l_size]
#                     V_all_ind = argmax(V_all)
#                     @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = maximum(V_all)
#                     @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all_ind[1]
#                     @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all_ind[2]
#                     @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all_ind[3]
#                 end
#             end
#         end
#     end
#     nothing
# end

# function sol_4!(variables::Mutable_Variables, parameters::NamedTuple)

#     # unpack parameters
#     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
#     @unpack ν_size, ν_grid, ν_Γ = parameters
#     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
#     @unpack n_size, n_grid, n_Γ, n_max = parameters
#     @unpack a_size, a_grid = parameters
#     @unpack l_size, l_grid, x_size, x_grid = parameters
#     @unpack inf_size, inf_grid = parameters
#     @unpack h_grid = parameters
#     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

#     for age_i = age_size:(-1):1
#         age = age_grid[age_i]
#         age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
#         age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
#         if age == age_ret # at retirement age
#             for n_i = 1:n_size
#                 n = n_grid[n_i]
#                 for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, a_i = 1:a_size
#                     ν = ν_grid[ν_i]
#                     ϵ = ϵ_grid[ϵ_i]
#                     a = a_grid[a_i]
#                     w = exp(h + ϵ + ν)
#                     if n == 0
#                         EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
#                         V_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV
#                         V_max_i = argmax(V_all)
#                         @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
#                         @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_max_i # a_grid[V_max_i]
#                     else
#                         a_p = reshape(a_grid, (a_size, 1, 1))
#                         x = reshape(x_grid, (1, x_size, 1))
#                         l = reshape(l_grid, (1, 1, l_size))
#                         q = quality_function.(x, l, Ref(n), Ref(μ), Ref(θ), Ref(ψ_1), Ref(ψ_2))
#                         c = (1.0 + r) * a .+ (1.0 .- l) * w .- a_p .- q_x .* x
#                         # V_all = [c[i,j,k] > 0 && q[1,j,k] >= q_bar ? c[i,j,k]^(1.0 - γ) / (1.0 - γ) + ψ * (n * q[1,j,k])^(1.0 - κ)/ (1.0 - κ) + β * variables.V[i, 1, ϵ_i, ν_i, 2, age_i+1] : -10^16 for i = 1:a_size, j = 1:x_size, k = 1:l_size]
#                         V_all = [c[i,j,k] > 0 && q[1,j,k] >= q_bar ? 1.0 / ((1.0 - γ) * c[i,j,k]^(γ - 1.0)) + ψ / ((1.0 - κ) * (n * q[1,j,k])^(κ - 1.0)) + β * variables.V[i, 1, ϵ_i, ν_i, 2, age_i+1] : -10^16 for i = 1:a_size, j = 1:x_size, k = 1:l_size]
#                         V_all_ind = argmax(V_all)
#                         @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = maximum(V_all)
#                         @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all_ind[1]
#                         @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all_ind[2]
#                         @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all_ind[3]
#                     end
#                 end
#             end
#         end
#     end
#     nothing
# end

# function sol_5!(variables::Mutable_Variables, parameters::NamedTuple)

#     # unpack parameters
#     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
#     @unpack ν_size, ν_grid, ν_Γ = parameters
#     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
#     @unpack n_size, n_grid, n_Γ, n_max = parameters
#     @unpack a_size, a_grid = parameters
#     @unpack l_size, l_grid, x_size, x_grid = parameters
#     @unpack inf_size, inf_grid = parameters
#     @unpack h_grid = parameters
#     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

#     for age_i = age_size:(-1):1
#         age = age_grid[age_i]
#         age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
#         age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
#         if age == age_ret # at retirement age
#             for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
#                 ν = ν_grid[ν_i]
#                 ϵ = ϵ_grid[ϵ_i]
#                 n = n_grid[n_i]
#                 a = a_grid[a_i]
#                 w = exp(h + ϵ + ν)
#                 if n == 0
#                     EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
#                     V_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV
#                     V_max_i = argmax(V_all)
#                     @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
#                     @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_max_i # a_grid[V_max_i]
#                 else
#                     Sol_all_i = 0
#                     for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
#                         a_p = a_grid[a_p_i]
#                         x = x_grid[x_i]
#                         l = l_grid[l_i]
#                         q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
#                         temp = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * variables.V[a_p_i, 1, ϵ_i, ν_i, 2, age_i+1]
#                         if Sol_all_i == 0 
#                             @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = temp
#                             @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = a_p_i
#                             @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = x_i
#                             @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = l_i
#                         elseif temp > variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i]
#                             @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = temp
#                             @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = a_p_i
#                             @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = x_i
#                             @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = l_i
#                         end
#                         Sol_all_i += 1
#                     end
#                 end
#             end
#         end
#     end
#     nothing
# end

# function sol_6!(variables::Mutable_Variables, parameters::NamedTuple)

#     # unpack parameters
#     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
#     @unpack ν_size, ν_grid, ν_Γ = parameters
#     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
#     @unpack n_size, n_grid, n_Γ, n_max = parameters
#     @unpack a_size, a_grid = parameters
#     @unpack l_size, l_grid, x_size, x_grid = parameters
#     @unpack inf_size, inf_grid = parameters
#     @unpack h_grid = parameters
#     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

#     for age_i = age_size:(-1):1
#         age = age_grid[age_i]
#         age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
#         age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
#         if age == age_ret # at retirement age
#             for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
#                 ν = ν_grid[ν_i]
#                 ϵ = ϵ_grid[ϵ_i]
#                 n = n_grid[n_i]
#                 a = a_grid[a_i]
#                 w = exp(h + ϵ + ν)
#                 if n == 0
#                     EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
#                     V_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV
#                     V_max_i = argmax(V_all)
#                     @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
#                     @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_max_i # a_grid[V_max_i]
#                 else
#                     function Sol_all!(variables::Mutable_Variables, parameters::NamedTuple)
#                         # unpack parameters
#                         @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
#                         @unpack ν_size, ν_grid, ν_Γ = parameters
#                         @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
#                         @unpack n_size, n_grid, n_Γ, n_max = parameters
#                         @unpack a_size, a_grid = parameters
#                         @unpack l_size, l_grid, x_size, x_grid = parameters
#                         @unpack inf_size, inf_grid = parameters
#                         @unpack h_grid = parameters
#                         @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters
#                         Sol_all_i = 0
#                         for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
#                             a_p = a_grid[a_p_i]
#                             x = x_grid[x_i]
#                             l = l_grid[l_i]
#                             q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
#                             temp = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * variables.V[a_p_i, 1, ϵ_i, ν_i, 2, age_i+1]
#                             if Sol_all_i == 0 
#                                 @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = temp
#                                 @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = a_p_i
#                                 @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = x_i
#                                 @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = l_i
#                             elseif temp > variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i]
#                                 @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = temp
#                                 @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = a_p_i
#                                 @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = x_i
#                                 @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = l_i
#                             end
#                             Sol_all_i += 1
#                         end
#                         nothing
#                     end
#                     Sol_all!(variables, parameters)
#                 end
#             end
#         end
#     end
#     nothing
# end

# function sol_7!(variables::Mutable_Variables, parameters::NamedTuple)

#     # unpack parameters
#     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
#     @unpack ν_size, ν_grid, ν_Γ = parameters
#     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
#     @unpack n_size, n_grid, n_Γ, n_max = parameters
#     @unpack a_size, a_grid = parameters
#     @unpack l_size, l_grid, x_size, x_grid = parameters
#     @unpack inf_size, inf_grid = parameters
#     @unpack h_grid = parameters
#     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

#     age_i = 48
#         age = age_grid[age_i]
#         age_ret_i = findall(parameters.age_grid .== parameters.age_ret)[]
#         age > age_ret ? h = h_grid[age_ret_i] : h = h_grid[age_i]
#         if age == age_ret # at retirement age
#             for ν_i = 1:ν_size, ϵ_i = 1:ϵ_size, n_i = 1:n_size, a_i = 1:a_size
#                 @views ν = ν_grid[ν_i]
#                 @views ϵ = ϵ_grid[ϵ_i]
#                 @views n = n_grid[n_i]
#                 @views a = a_grid[a_i]
#                 w = exp(h + ϵ + ν)
#                 if n == 0
#                     @views EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
#                     @views V_all = utility_function.((1.0 + r) * a + w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β * EV
#                     @views V_max_i = argmax(V_all)
#                     @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
#                     @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_max_i
#                 else
#                     Sol_all_i = 0
#                     for a_p_i = 1:a_size, x_i = 1:x_size, l_i = 1:l_size
#                         @views a_p = a_grid[a_p_i]
#                         @views x = x_grid[x_i]
#                         @views l = l_grid[l_i]
#                         q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
#                         @views temp = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) + β * variables.V[a_p_i, 1, ϵ_i, ν_i, 2, age_i+1]
#                         if Sol_all_i == 0 
#                             @views @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = temp
#                             @views @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = a_p_i
#                             @views @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = x_i
#                             @views @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = l_i
#                         elseif temp > variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i]
#                             @views @inbounds variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = temp
#                             @views @inbounds variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = a_p_i
#                             @views @inbounds variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = x_i
#                             @views @inbounds variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = l_i
#                         end
#                         Sol_all_i += 1
#                     end
#                 end
#             end
#         end
#     nothing
# end

# function sol_8!(variables::Mutable_Variables, parameters::NamedTuple)
#     # Unpack parameters
#     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
#     @unpack ν_size, ν_grid, ν_Γ = parameters
#     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
#     @unpack n_size, n_grid, n_Γ, n_max = parameters
#     @unpack a_size, a_grid = parameters
#     @unpack l_size, l_grid, x_size, x_grid = parameters
#     @unpack inf_size, inf_grid = parameters
#     @unpack h_grid = parameters
#     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

#     age_i = 48
#     age = age_grid[age_i]
#     age_ret_i = findfirst(parameters.age_grid .== parameters.age_ret)
#     h = age > age_ret ? h_grid[age_ret_i] : h_grid[age_i]

#     if age == age_ret # At retirement age
#         for ν_i in 1:ν_size
#             @views ν = ν_grid[ν_i]
#             for ϵ_i in 1:ϵ_size
#                 @views ϵ = ϵ_grid[ϵ_i]
#                 w = exp(h + ϵ + ν)
#                 for n_i in 1:n_size
#                     @views n = n_grid[n_i]
#                     for a_i in 1:a_size
#                         @views a = a_grid[a_i]
#                         if n == 0
#                             @views EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
#                             V_all = utility_function.((1.0 + r) * a .+ w .- a_grid, n, 0.0, γ, ψ, κ, q_bar) .+ β .* EV
#                             V_max_i = argmax(V_all)
#                             @inbounds begin
#                                 variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
#                                 variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_max_i
#                             end
#                         else
#                             Sol_all_i = 0
#                             V_best = -Inf
#                             best_a_p_i, best_x_i, best_l_i = 0, 0, 0
#                             for a_p_i in 1:a_size
#                                 @views a_p = a_grid[a_p_i]
#                                 for x_i in 1:x_size
#                                     @views x = x_grid[x_i]
#                                     for l_i in 1:l_size
#                                         @views l = l_grid[l_i]
#                                         q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
#                                         temp = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) +
#                                                β * variables.V[a_p_i, 1, ϵ_i, ν_i, 2, age_i+1]
#                                         if temp > V_best
#                                             V_best = temp
#                                             best_a_p_i = a_p_i
#                                             best_x_i = x_i
#                                             best_l_i = l_i
#                                         end
#                                         Sol_all_i += 1
#                                     end
#                                 end
#                             end
#                             @inbounds begin
#                                 variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_best
#                                 variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_a_p_i
#                                 variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_x_i
#                                 variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_l_i
#                             end
#                         end
#                     end
#                 end
#             end
#         end
#     end
#     nothing
# end
# function sol_9!(variables::Mutable_Variables, parameters::NamedTuple)
#     # Unpack parameters
#     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
#     @unpack ν_size, ν_grid, ν_Γ = parameters
#     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
#     @unpack n_size, n_grid, n_Γ, n_max = parameters
#     @unpack a_size, a_grid = parameters
#     @unpack l_size, l_grid, x_size, x_grid = parameters
#     @unpack inf_size, inf_grid = parameters
#     @unpack h_grid = parameters
#     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

#     age_i = 48
#     age = age_grid[age_i]
#     age_ret_i = findfirst(parameters.age_grid .== parameters.age_ret)
#     h = ifelse(age > age_ret, h_grid[age_ret_i], h_grid[age_i])

#     if age == age_ret # At retirement age
#         for ν_i in 1:ν_size
#             @views ν = ν_grid[ν_i]
#             for ϵ_i in 1:ϵ_size
#                 @views ϵ = ϵ_grid[ϵ_i]
#                 w = exp(h + ϵ + ν)
#                 for n_i in 1:n_size
#                     @views n = n_grid[n_i]
#                     for a_i in 1:a_size
#                         @views a = a_grid[a_i]
#                         if n == 0
#                             @views EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
#                             V_all = similar(EV)  # Ensure V_all is correctly initialized
#                             @inbounds @. V_all = utility_function((1.0 + r) * a + w - a_grid, n, 0.0, γ, ψ, κ, q_bar) + β * EV
#                             V_max_i = argmax(V_all)
#                             @inbounds begin
#                                 variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
#                                 variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_max_i
#                             end
#                         else
#                             V_best = -Inf
#                             best_a_p_i, best_x_i, best_l_i = 0, 0, 0
#                             for a_p_i in 1:a_size
#                                 a_p = a_grid[a_p_i]
#                                 for x_i in 1:x_size
#                                     x = x_grid[x_i]
#                                     for l_i in 1:l_size
#                                         l = l_grid[l_i]
#                                         q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
#                                         temp = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) +
#                                                β * variables.V[a_p_i, 1, ϵ_i, ν_i, 2, age_i+1]
#                                         if temp > V_best
#                                             V_best = temp
#                                             best_a_p_i = a_p_i
#                                             best_x_i = x_i
#                                             best_l_i = l_i
#                                         end
#                                     end
#                                 end
#                             end
#                             @inbounds begin
#                                 variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_best
#                                 variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_a_p_i
#                                 variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_x_i
#                                 variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_l_i
#                             end
#                         end
#                     end
#                 end
#             end
#         end
#     end
#     nothing
# end

# function sol_10!(variables::Mutable_Variables, parameters::NamedTuple)
#     # Unpack parameters
#     @unpack age_size, age_grid, age_max, age_ret, age_inf = parameters
#     @unpack ν_size, ν_grid, ν_Γ = parameters
#     @unpack ϵ_size, ϵ_grid, ϵ_Γ = parameters
#     @unpack n_size, n_grid, n_Γ, n_max = parameters
#     @unpack a_size, a_grid = parameters
#     @unpack l_size, l_grid, x_size, x_grid = parameters
#     @unpack inf_size, inf_grid = parameters
#     @unpack h_grid = parameters
#     @unpack b, r, γ, ψ, κ, β, μ, θ, ψ_1, ψ_2, q_bar, q_x = parameters

#     age_i = 48
#     age = age_grid[age_i]
#     age_ret_i = findfirst(parameters.age_grid .== parameters.age_ret)
#     h = ifelse(age > age_ret, h_grid[age_ret_i], h_grid[age_i])

#     # Precompute (1.0 + r) * a_grid
#     # r_a_grid = (1.0 + r) .* a_grid

#     # Preallocate array for V_all
#     V_all = Vector{Float64}(undef, a_size)

#     if age == age_ret # At retirement age
#         for ν_i in 1:ν_size
#             @views ν = ν_grid[ν_i]
#             for ϵ_i in 1:ϵ_size
#                 @views ϵ = ϵ_grid[ϵ_i]
#                 w = exp(h + ϵ + ν)
#                 for n_i in 1:n_size
#                     @views n = n_grid[n_i]
#                     for a_i in 1:a_size
#                         @views a = a_grid[a_i]
#                         if n == 0
#                             # Avoid allocation in utility_function
#                             @views EV = variables.V[:, 1, ϵ_i, ν_i, 2, age_i+1]
#                             # @inbounds @. V_all = utility_function(r_a_grid + w - a_grid, n, 0.0, γ, ψ, κ, q_bar) + β * EV
#                             @inbounds V_all .= utility_function.((1.0 + r) .* a .+ w .- a_grid, Ref(n), Ref(0.0), Ref(γ), Ref(ψ), Ref(κ), Ref(q_bar)) .+ β .* EV
#                             V_max_i = argmax(V_all)
#                             @inbounds begin
#                                 variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_all[V_max_i]
#                                 variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_max_i
#                             end
#                         else
#                             V_best = -Inf
#                             best_a_p_i, best_x_i, best_l_i = 0, 0, 0
#                             for a_p_i in 1:a_size
#                                 a_p = a_grid[a_p_i]
#                                 for x_i in 1:x_size
#                                     x = x_grid[x_i]
#                                     for l_i in 1:l_size
#                                         l = l_grid[l_i]
#                                         q = quality_function(x, l, n, μ, θ, ψ_1, ψ_2)
#                                         temp = utility_function((1.0 + r) * a + (1.0 - l) * w - a_p - q_x * x, n, q, γ, ψ, κ, q_bar) +
#                                                β * variables.V[a_p_i, 1, ϵ_i, ν_i, 2, age_i+1]
#                                         if temp > V_best
#                                             V_best = temp
#                                             best_a_p_i = a_p_i
#                                             best_x_i = x_i
#                                             best_l_i = l_i
#                                         end
#                                     end
#                                 end
#                             end
#                             @inbounds begin
#                                 variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i] = V_best
#                                 variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_a_p_i
#                                 variables.policy_x[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_x_i
#                                 variables.policy_l[a_i, n_i, ϵ_i, ν_i, 2, age_i] = best_l_i
#                             end
#                         end
#                     end
#                 end
#             end
#         end
#     end
#     nothing
# end

# variables_1 = variables_function(parameters);
# # variables_2 = variables_function(parameters);
# # variables_3 = variables_function(parameters);
# # variables_4 = variables_function(parameters);
# variables_5 = variables_function(parameters);
# # variables_6 = variables_function(parameters);
# variables_7 = variables_function(parameters);
# variables_8 = variables_function(parameters);
# variables_9 = variables_function(parameters);
# variables_10 = variables_function(parameters);

# sol_1!(variables_1, parameters)
# # sol_2!(variables_2, parameters)
# # sol_3!(variables_3, parameters)
# # sol_4!(variables_4, parameters)
# sol_5!(variables_5, parameters)
# # sol_6!(variables_6, parameters)
# sol_7!(variables_7, parameters)
# sol_8!(variables_8, parameters)
# sol_9!(variables_9, parameters)
# sol_10!(variables_10, parameters)

# # @assert variables_1.V[:,:,:,:,2,48] == variables_2.V[:,:,:,:,2,48] == variables_3.V[:,:,:,:,2,48] == variables_4.V[:,:,:,:,2,48]
# # @assert variables_1.policy_a_p == variables_2.policy_a_p == variables_3.policy_a_p == variables_4.policy_a_p
# # @assert variables_1.policy_x == variables_2.policy_x == variables_3.policy_x == variables_4.policy_x
# # @assert variables_1.policy_l == variables_2.policy_l == variables_3.policy_l == variables_4.policy_l

# @assert variables_1.V[:,:,:,:,2,48] == variables_5.V[:,:,:,:,2,48] == variables_7.V[:,:,:,:,2,48] == variables_8.V[:,:,:,:,2,48] == variables_9.V[:,:,:,:,2,48] == variables_10.V[:,:,:,:,2,48]
# @assert variables_1.policy_a_p == variables_5.policy_a_p == variables_7.policy_a_p == variables_8.policy_a_p == variables_9.policy_a_p == variables_10.policy_a_p
# @assert variables_1.policy_x == variables_5.policy_x == variables_7.policy_x == variables_8.policy_x == variables_9.policy_x == variables_10.policy_x
# @assert variables_1.policy_l == variables_5.policy_l == variables_7.policy_l == variables_8.policy_l == variables_9.policy_l == variables_10.policy_l

# @btime sol_1!($variables_1, $parameters)
# # @btime sol_2!($variables_2, $parameters)
# # @btime sol_3!($variables_3, $parameters)
# # @btime sol_4!($variables_4, $parameters)
# @btime sol_5!($variables_5, $parameters)
# # @btime sol_6!($variables_6, $parameters)
# @btime sol_7!($variables_7, $parameters)
# @btime sol_8!($variables_8, $parameters)
# @btime sol_9!($variables_9, $parameters)
# @btime sol_10!($variables_10, $parameters)







