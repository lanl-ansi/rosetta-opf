time_start = time()

using PowerModels
using NLPModels
using ADNLPModels
using NLPModelsIpopt

pkg_load_time = time() - time_start


time_start = time()

file_name = "../data/pglib_opf_case5_pjm.m"

data = PowerModels.parse_file(file_name)
PowerModels.standardize_cost_terms!(data, order=2)
PowerModels.calc_thermal_limits!(data)
ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

data_load_time = time() - time_start


time_start = time()

f(x) = (x[1] - 1)^2 + 100*(x[2] - x[1]^2)^2

x0 = [-1.2; 1.0]
uvar = [10.0; 10.0]
lvar = [-10.0; -10.0]

c(x) = [x[1]^2 + x[2]^2; x[1]*x[2]]
ucon = [0.75; 0.3]
lcon = [0.0; 0.0]

nlp = ADNLPModel(f, x0, lvar, uvar, c, lcon, ucon)
#ADNLPModel(f, x0, lvar, uvar)


model_build_time = time() - time_start


time_start = time()

#output = ipopt(nlp, print_level=0)
output = ipopt(nlp)
cost = output.objective
println(output.solution)

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

