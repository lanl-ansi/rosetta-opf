# https://github.com/SciML/Optimization.jl/blob/master/lib/OptimizationMOI/test/runtests.jl

using Optimization
using OptimizationMOI
using ForwardDiff
using Ipopt

rosenbrock(x,p) =  (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
p  = Dict(1 => 1.0, 2=> 100.0)

function con2_c(x,p)
    [x[1]^2 + x[2]^2, x[2]*sin(x[1])-x[1]]
end

optprob = OptimizationFunction(rosenbrock, Optimization.AutoForwardDiff(); cons=con2_c)
prob = OptimizationProblem(optprob, x0, p, lcons=[-Inf,-Inf], ucons=[Inf,Inf])

sol = solve(prob, Ipopt.Optimizer())

println(sol.minimum)
println(sol.u)
