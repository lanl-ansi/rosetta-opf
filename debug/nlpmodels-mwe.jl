# https://juliasmoothoptimizers.github.io/ADNLPModels.jl/stable/tutorial/

using ADNLPModels
using NLPModelsIpopt

f(x) = (x[1] - 1)^2 + 100*(x[2] - x[1]^2)^2

x0 = [-1.2; 1.0]
uvar = [10.0; 10.0]
lvar = [-10.0; -10.0]

c(x) = [x[1]^2 + x[2]^2; x[1]*x[2]]
ucon = [0.75; 0.3]
lcon = [0.0; 0.0]

nlp = ADNLPModel(f, x0, lvar, uvar, c, lcon, ucon)

#output = ipopt(nlp, print_level=0)
output = ipopt(nlp)

println(output.objective)
println(output.solution)
