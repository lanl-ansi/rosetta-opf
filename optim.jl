# https://github.com/JuliaNLSolvers/Optim.jl

# apparently supports optimization with manifold constraints
# unable to find an example

using Optim

rosenbrock(x) =  (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

result = optimize(rosenbrock, zeros(2), BFGS())
