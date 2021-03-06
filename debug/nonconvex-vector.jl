time_start = time()

using PowerModels
using Nonconvex
Nonconvex.@load Ipopt

pkg_load_time = time() - time_start


time_start = time()

file_name = "../data/pglib_opf_case5_pjm.m"
#file_name = "../data/pglib_opf_case118_ieee.m"

data = PowerModels.parse_file(file_name)
PowerModels.standardize_cost_terms!(data, order=2)
PowerModels.calc_thermal_limits!(data)
ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

data_load_time = time() - time_start


time_start = time()

model = Nonconvex.Model()


bus_pd = Dict(i => 0.0 for (i,bus) in ref[:bus])
bus_qd = Dict(i => 0.0 for (i,bus) in ref[:bus])

bus_gs = Dict(i => 0.0 for (i,bus) in ref[:bus])
bus_bs = Dict(i => 0.0 for (i,bus) in ref[:bus])

for (i,bus) in ref[:bus]
    if length(ref[:bus_loads][i]) > 0
        bus_pd[i] = sum(ref[:load][l]["pd"] for l in ref[:bus_loads][i])
        bus_qd[i] = sum(ref[:load][l]["qd"] for l in ref[:bus_loads][i])
    end

    if length(ref[:bus_shunts][i]) > 0
        bus_gs[i] = sum(ref[:shunt][s]["gs"] for s in ref[:bus_shunts][i])
        bus_bs[i] = sum(ref[:shunt][s]["bs"] for s in ref[:bus_shunts][i])
    end
end


br_g = Dict(i => 0.0 for (i,branch) in ref[:branch])
br_b = Dict(i => 0.0 for (i,branch) in ref[:branch])

br_tr = Dict(i => 0.0 for (i,branch) in ref[:branch])
br_ti = Dict(i => 0.0 for (i,branch) in ref[:branch])
br_tm = Dict(i => 0.0 for (i,branch) in ref[:branch])

br_g_fr = Dict(i => 0.0 for (i,branch) in ref[:branch])
br_b_fr = Dict(i => 0.0 for (i,branch) in ref[:branch])
br_g_to = Dict(i => 0.0 for (i,branch) in ref[:branch])
br_b_to = Dict(i => 0.0 for (i,branch) in ref[:branch])

for (i,branch) in ref[:branch]
    g, b = PowerModels.calc_branch_y(branch)
    tr, ti = PowerModels.calc_branch_t(branch)

    br_g[i] = g
    br_b[i] = b

    br_tr[i] = tr
    br_ti[i] = ti
    br_tm[i] = tr^2 + ti^2

    br_g_fr[i] = branch["g_fr"]
    br_b_fr[i] = branch["b_fr"]
    br_g_to[i] = branch["g_to"]
    br_b_to[i] = branch["b_to"]
end


var_lookup = Dict{String,Int}()

var_idx = 1
for (i,bus) in ref[:bus]
    addvar!(model, -Inf, Inf) #va
    var_lookup["va_$(i)"] = var_idx
    global var_idx += 1

    addvar!(model, bus["vmin"], bus["vmax"], init=1.0) #vm
    var_lookup["vm_$(i)"] = var_idx
    global var_idx += 1
end

for (i,gen) in ref[:gen]
    addvar!(model, gen["pmin"], gen["pmax"]) #pg
    var_lookup["pg_$(i)"] = var_idx
    global var_idx += 1

    addvar!(model, gen["qmin"], gen["qmax"]) #qg
    var_lookup["qg_$(i)"] = var_idx
    global var_idx += 1
end

for (l,i,j) in ref[:arcs]
    branch = ref[:branch][l]

    addvar!(model, -branch["rate_a"], branch["rate_a"]) #p
    var_lookup["p_$(l)_$(i)_$(j)"] = var_idx
    global var_idx += 1

    addvar!(model, -branch["rate_a"], branch["rate_a"]) #q
    var_lookup["q_$(l)_$(i)_$(j)"] = var_idx
    global var_idx += 1
end

@assert var_idx == length(var_lookup)+1


function opf_objective(x::Vector)
    cost = 0.0
    for (i,gen) in ref[:gen]
        pg = x[var_lookup["pg_$(i)"]]
        cost += gen["cost"][1]*pg^2 + gen["cost"][2]*pg + gen["cost"][3]
    end
    return cost
end

set_objective!(model, opf_objective)

#println(model)

model_build_time = time() - time_start


time_start = time()

alg = IpoptAlg()
options = IpoptOptions(print_level = 0)
x0 = NonconvexCore.getinit(model)
#r = optimize(model, alg, x0, options = options)

# Nonconvex v1.0.2
# ERROR: LoadError: MethodError: no method matching getindex(::Dict{Any, Any})
# something to do with :Zygote it seems
r = optimize(model, alg, x0)

solve_time = time() - time_start


println("")
println("\033[1mSummary\033[0m")
println("   case......: $(file_name)")
#println("   cost......: $(round(Int, cost))")
println("   pkg time..: $(pkg_load_time)")
println("   data time.: $(data_load_time)")
println("   build time: $(model_build_time)")
println("   solve time: $(solve_time)")
println("")

