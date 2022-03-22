# https://github.com/JuliaNLSolvers/Optim.jl

# WIP

# apparently supports optimization with manifold constraints
# unable to find an example...
# https://julianlsolvers.github.io/Optim.jl/stable/#examples/generated/ipnewton_basics/

time_start = time()

using Optim

pkg_load_time = time() - time_start


time_start = time()

file_name = "../data/pglib_opf_case5_pjm.m"
#file_name = "../data/pglib_opf_case118_ieee.m"

data_load_time = time() - time_start


time_start = time()

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

model_build_time = time() - time_start

time_start = time()

res = optimize(df, dfc, x0, IPNewton(), Optim.Options(show_trace=true))
display(res)

sol = res.minimizer
cost = res.minimum

println(sol)

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
