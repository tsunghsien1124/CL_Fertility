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
    a = reshape(a_grid, (a_size, 1, 1, 1, 1))
    b = reshape(b_grid, (1, b_size, 1, 1, 1))
    c = reshape(c_grid, (1, 1, c_size, 1, 1))
    d = reshape(d_grid, (1, 1, 1, d_size, 1))
    e = reshape(e_grid, (1, 1, 1, 1, e_size))
    return a .+ b .- c .* d .- e
end

function sol_loop(a_size, a_grid, b_size, b_grid, c_size, c_grid, d_size, d_grid, e_size, e_grid)
    res = zeros(a_size,b_size,c_size,d_size,e_size)
    # Threads.@threads for (a_i, b_i, c_i) in collect(Iterators.product(1:a_size, 1:b_size, 1:c_size))
    for a_i = 1:a_size, b_i = 1:b_size, c_i = 1:c_size, d_i = 1:d_size, e_i = 1:e_size
        res[a_i,b_i,c_i,d_i,e_i] = a_grid[a_i] + b_grid[b_i] - c_grid[c_i] * d_grid[d_i] - e_grid[e_i]
    end
    return res
end

function sol_loop_inline(a_size, a_grid, b_size, b_grid, c_size, c_grid, d_size, d_grid, e_size, e_grid)
    return [a_grid[a_i] + b_grid[b_i] - c_grid[c_i] * d_grid[d_i] - e_grid[e_i] for a_i = 1:a_size, b_i = 1:b_size, c_i = 1:c_size, d_i = 1:d_size, e_i = 1:e_size]
    # return [a_grid[a_i] + b_grid[b_i] - c_grid[c_i] for (a_i, b_i, c_i) in collect(Iterators.product(1:a_size, 1:b_size, 1:c_size))]
end

sol_1 = sol_reshape(a_size, a_grid, b_size, b_grid, c_size, c_grid, d_size, d_grid, e_size, e_grid);
sol_2 = sol_loop(a_size, a_grid, b_size, b_grid, c_size, c_grid, d_size, d_grid, e_size, e_grid);
sol_3 = sol_loop_inline(a_size, a_grid, b_size, b_grid, c_size, c_grid, d_size, d_grid, e_size, e_grid)

@btime sol_reshape($a_size, $a_grid, $b_size, $b_grid, $c_size, $c_grid, $d_size, $d_grid, $e_size, $e_grid);
@btime sol_loop($a_size, $a_grid, $b_size, $b_grid, $c_size, $c_grid, $d_size, $d_grid, $e_size, $e_grid);
@btime sol_loop_inline($a_size, $a_grid, $b_size, $b_grid, $c_size, $c_grid, $d_size, $d_grid, $e_size, $e_grid);