# https://julianlsolvers.github.io/Optim.jl/stable/#examples/generated/ipnewton_basics/

using Optim

fun(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

function con2_c!(c, x)
    c[1] = x[1]^2 + x[2]^2     ## First constraint
    c[2] = x[2]*sin(x[1])-x[1] ## Second constraint
    return c
end

x0 = [0.25, 0.25]
df = TwiceDifferentiable(fun, x0)

lc = [-Inf, 0.0]; uc = [0.5^2, 0.0]
lx = [-0.5, -0.5]; ux = [0.5, 0.5]
dfc = TwiceDifferentiableConstraints(con2_c!, lx, ux, lc, uc)

res = optimize(df, dfc, x0, IPNewton(), Optim.Options(show_trace=true))
display(res)

println(res.minimum)
println(res.minimizer)
