# See also
# https://juliamanifolds.github.io/ManoptExamples.jl/stable/examples/Rosenbrock/
# and
# https://manoptjl.org/stable/tutorials/ConstrainedOptimization/

#
# Trying to reproduce the Optim.jl MWE
#

using Manopt, Manifolds

# Rosenbrock cost
f(M,x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

# the gradient can be improved by using an in-place one
# one could also use a scheme with AD like
# https://manoptjl.org/stable/tutorials/AutomaticDifferentiation/
# or really just AD, since we are Euclidean
function grad_f(M,x)
    return [
        -2*(1.0-x[1]) - 400*(x[2] - x[1]^2)*x[1],
        200.0*(x[2]-x[1]^2)
    ]
end

x0 = [0.25, 0.25]
M = Euclidean(2)

# Then one can for example use
x1 = quasi_Newton(M, f, grad_f, x0)
# Or look at some output as well
quasi_Newton(M, f, grad_f, x0;
    debug=[:Iteration, :Cost, " ", :GradientNorm, "\n", 10, :Stop],
    return_state=true
)
# Or even change the stopping criterion
quasi_Newton(M, f, grad_f, x0;
    debug=[:Iteration, :Cost, " ", :GradientNorm, "\n", 10, :Stop],
    stopping_criterion = StopAfterIteration(1000) | StopWhenGradientNormLess(1e-9),
    return_state=true
)


# Constraints
# I do not understand what the constraints in the Optim.jl MWE mean,
# TwiceDifferentiableConstraints is no documented,
# * so I _guess_ lx and ux are upper and lower bounds for x?
# * what are lc and uc? Upper and lower bounds ... for...?

# So for now just the two constraints that I do understand
# at least I hope that they are meant <= 0 ?
# If we need upper and lower bounds we would have to hardcode them here,
# i.e. two lines each
function g(M,x)
    return [
        x[1]^2 + x[2]^2,
        x[2]*sin(x[1])-x[1]
    ]
end

# alternatively one could also implement them separately
# Then the ALM we have does need the gradient of G of course

augmented_Lagrangian_method(M, f, grad_f; g=g, grad_g=???)