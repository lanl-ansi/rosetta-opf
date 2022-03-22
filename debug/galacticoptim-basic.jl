time_start = time()

using GalacticOptim
using ForwardDiff
using Ipopt

pkg_load_time = time() - time_start


time_start = time()

data_load_time = time() - time_start


time_start = time()

# Works

rosenbrock(x,p) =  (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
_p  = Dict(1 => 1.0, 2=> 100.0)

function con2_c(x,p)
    [x[1]^2 + x[2]^2, x[2]*sin(x[1])-x[1]]
end

optprob = OptimizationFunction(rosenbrock, GalacticOptim.AutoForwardDiff(); cons= con2_c)
prob = OptimizationProblem(optprob, x0, _p, lcons = [-Inf,-Inf], ucons = [Inf,Inf])


# Broken?
using GalacticOptim, ForwardDiff, Ipopt

function objective(x,p)
    cost = 0.0
    for (i,v) in p
        var = x[i]
        cost += (v-var)^2
    end
    return cost
end

x0 = [0.0, 0.0, 0.0, 0.0]
parameters  = Dict(2 => 2.0, 4 => 4.0)

optprob = OptimizationFunction(objective, GalacticOptim.AutoForwardDiff())
prob = OptimizationProblem(optprob, x0, parameters)
sol = solve(prob, Ipopt.Optimizer())


model_build_time = time() - time_start


time_start = time()

sol = solve(prob, Ipopt.Optimizer())
cost = sol.minimum

solve_time = time() - time_start

println("")
println("\033[1mSummary\033[0m")
println("   case......: $(file_name)")
println("   cost......: $(round(Int, cost))")
println("   pkg time..: $(pkg_load_time)")
println("   data time.: $(data_load_time)")
println("   build time: $(model_build_time)")
println("   solve time: $(solve_time)")
println("")

