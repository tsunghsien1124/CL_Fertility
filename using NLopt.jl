using NLopt

# example 1
function myfunc(x::Vector, grad::Vector)
    # if length(grad) > 0
    #     grad[1] = 0
    #     grad[2] = 0.5/sqrt(x[2])
    # end
    return sqrt(x[2])
end

function myconstraint(x::Vector, grad::Vector)
    # if length(grad) > 0
    #     grad[1] = 3a * (a*x[1] + b)^2
    #     grad[2] = -1
    # end
    return (a*x[1] + b)^3 - x[2]
end

a = 2
b = 0

opt = Opt(:LN_COBYLA, 2)
# opt = Opt(:LD_MMA, 2)
opt.lower_bounds = [-Inf, 0.]
opt.xtol_rel = 1e-4

opt.min_objective = myfunc
inequality_constraint!(opt, myconstraint, 1e-8)
# inequality_constraint!(opt, (x,g) -> myconstraint(x,g,-1,1), 1e-8)

(minf,minx,ret) = optimize(opt, [5.234, 5.678])
numevals = opt.numevals # the number of function evaluations
println("got $minf at $minx after $numevals iterations (returned $ret)")

# example 2
function myfunc(x::Vector, grad::Vector)
    return x[1]^2
end
opt = Opt(:LN_COBYLA, 1)
opt.lower_bounds = [0.0]
opt.xtol_rel = 1e-4
opt.min_objective = myfunc
(minf,minx,ret) = optimize(opt, [10^10])
numevals = opt.numevals
println("got $minf at $minx after $numevals iterations (returned $ret)")

