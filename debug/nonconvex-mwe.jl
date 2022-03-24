# https://julianonconvex.github.io/Nonconvex.jl/stable/
# concerning that Model and DictModel find different solutions

using Nonconvex
Nonconvex.@load Ipopt

f(x) = sqrt(x[2])
g(x, a, b) = (a*x[1] + b)^3 - x[2]

model = Model(f)

addvar!(model, 0.0, 10.0, init=1.0)
addvar!(model, 0.0, 10.0, init=1.0)

add_ineq_constraint!(model, x -> g(x, 2, 0))
add_ineq_constraint!(model, x -> g(x, -1, 1))

x0 = NonconvexCore.getinit(model)
r = optimize(model, IpoptAlg(), x0)

println(r.minimum)
println(r.minimizer)



f1(x) = sqrt(x["2"])
g1(x, a, b) = (a*x["1"] + b)^3 - x["2"]

model = DictModel(f1)

addvar!(model, "1", 0.0, 10.0, init=1.0)
addvar!(model, "2", 0.0, 10.0, init=1.0)

add_ineq_constraint!(model, x -> g1(x, 2, 0))
add_ineq_constraint!(model, x -> g1(x, -1, 1))

x0 = NonconvexCore.getinit(model)
r = optimize(model, IpoptAlg(), x0)

println(r.minimum)
println(r.minimizer)
