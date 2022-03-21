using GalacticOptim


rosenbrock(x,p) =  (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
_p  = [1.0, 100.0]
l1 = rosenbrock(x0, _p)

function con2_c(x,p)
    [x[1]^2 + x[2]^2, x[2]*sin(x[1])-x[1]]
end

optprob = OptimizationFunction(rosenbrock, GalacticOptim.AutoForwardDiff(); cons= con2_c)
prob = OptimizationProblem(optprob, x0, _p, lcons = [-Inf,-Inf], ucons = [Inf,Inf])
sol = solve(prob, Ipopt.Optimizer())



using GalacticOptim
using Optim
using ForwardDiff

rosenbrock(x,p) =  (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
p  = [1.0,100.0]

f = OptimizationFunction(rosenbrock, GalacticOptim.AutoForwardDiff())
prob = OptimizationProblem(f, x0, p)
sol = solve(prob, BFGS())
