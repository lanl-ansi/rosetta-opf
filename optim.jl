# https://github.com/JuliaNLSolvers/Optim.jl
# https://julianlsolvers.github.io/Optim.jl/stable/#examples/generated/ipnewton_basics/

time_start = time()

using PowerModels
using Optim

pkg_load_time = time() - time_start


time_start = time()

file_name = "data/pglib_opf_case5_pjm.m"
#file_name = "data/pglib_opf_case118_ieee.m"

data = PowerModels.parse_file(file_name)
PowerModels.standardize_cost_terms!(data, order=2)
PowerModels.calc_thermal_limits!(data)
ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

data_load_time = time() - time_start


time_start = time()


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

var_init = Float64[]
var_lb = Float64[]
var_ub = Float64[]

var_idx = 1
for (i,bus) in ref[:bus]
    push!(var_init, 0.0) #va
    push!(var_lb, -Inf)
    push!(var_ub, Inf)
    var_lookup["va_$(i)"] = var_idx
    global var_idx += 1

    push!(var_init, 1.0) #vm
    push!(var_lb, bus["vmin"])
    push!(var_ub, bus["vmax"])
    var_lookup["vm_$(i)"] = var_idx
    global var_idx += 1
end

for (i,gen) in ref[:gen]
    #push!(var_init, 0.0) #pg
    push!(var_init, (gen["pmax"]+gen["pmin"])/2) # non-standard start
    push!(var_lb, gen["pmin"])
    push!(var_ub, gen["pmax"])
    var_lookup["pg_$(i)"] = var_idx
    global var_idx += 1

    #push!(var_init, 0.0) #qg
    push!(var_init, (gen["qmax"]+gen["qmin"])/2) # non-standard start
    push!(var_lb, gen["qmin"])
    push!(var_ub, gen["qmax"])
    var_lookup["qg_$(i)"] = var_idx
    global var_idx += 1
end

for (l,i,j) in ref[:arcs]
    branch = ref[:branch][l]

    push!(var_init, 0.0) #p
    push!(var_lb, -branch["rate_a"])
    push!(var_ub,  branch["rate_a"])
    var_lookup["p_$(l)_$(i)_$(j)"] = var_idx
    global var_idx += 1

    push!(var_init, 0.0) #q
    push!(var_lb, -branch["rate_a"])
    push!(var_ub,  branch["rate_a"])
    var_lookup["q_$(l)_$(i)_$(j)"] = var_idx
    global var_idx += 1
end

@assert var_idx == length(var_init)+1

function opf_objective(x)
    cost = 0.0
    for (i,gen) in ref[:gen]
        pg = x[var_lookup["pg_$(i)"]]
        cost += gen["cost"][1]*pg^2 + gen["cost"][2]*pg + gen["cost"][3]
    end
    return cost
end

function opf_constraints(c,x)
    va = Dict(i => x[var_lookup["va_$(i)"]] for (i,bus) in ref[:bus])
    vm = Dict(i => x[var_lookup["vm_$(i)"]] for (i,bus) in ref[:bus])

    pg = Dict(i => x[var_lookup["pg_$(i)"]] for (i,gen) in ref[:gen])
    qg = Dict(i => x[var_lookup["qg_$(i)"]] for (i,gen) in ref[:gen])

    p = Dict((l,i,j) => x[var_lookup["p_$(l)_$(i)_$(j)"]] for (l,i,j) in ref[:arcs])
    q = Dict((l,i,j) => x[var_lookup["q_$(l)_$(i)_$(j)"]] for (l,i,j) in ref[:arcs])

    vm_fr = Dict(l => vm[branch["f_bus"]] for (l,branch) in ref[:branch])
    vm_to = Dict(l => vm[branch["t_bus"]] for (l,branch) in ref[:branch])
    va_fr = Dict(l => va[branch["f_bus"]] for (l,branch) in ref[:branch])
    va_to = Dict(l => va[branch["t_bus"]] for (l,branch) in ref[:branch])


    va_con = [va[i] for (i,bus) in ref[:ref_buses]]

    #     @constraint(model,
    #         sum(p[a] for a in ref[:bus_arcs][i]) ==
    #         sum(pg[g] for g in ref[:bus_gens][i]) -
    #         sum(load["pd"] for load in bus_loads) -
    #         sum(shunt["gs"] for shunt in bus_shunts)*vm[i]^2
    #     )
    power_balance_p_con = [
       sum(pg[j] for j in ref[:bus_gens][i]; init=0.0) -
       bus_pd[i] -
       bus_gs[i]*vm[i]^2 -
       sum(p[a] for a in ref[:bus_arcs][i])
       for (i,bus) in ref[:bus]
    ]

    #     @constraint(model,
    #         sum(q[a] for a in ref[:bus_arcs][i]) ==
    #         sum(qg[g] for g in ref[:bus_gens][i]) -
    #         sum(load["qd"] for load in bus_loads) +
    #         sum(shunt["bs"] for shunt in bus_shunts)*vm[i]^2
    #     )
    power_balance_q_con = [
       sum(qg[j] for j in ref[:bus_gens][i]; init=0.0) -
       bus_qd[i] +
       bus_bs[i]*vm[i]^2 -
       sum(q[a] for a in ref[:bus_arcs][i])
       for (i,bus) in ref[:bus]
    ]


    # @NLconstraint(model, p_fr ==  (g+g_fr)/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
    power_flow_p_from_con = [
       (br_g[l]+br_g_fr[l])/br_tm[l]*vm_fr[l]^2 +
       (-br_g[l]*br_tr[l]+br_b[l]*br_ti[l])/br_tm[l]*(vm_fr[l]*vm_to[l]*cos(va_fr[l]-va_to[l])) +
       (-br_b[l]*br_tr[l]-br_g[l]*br_ti[l])/br_tm[l]*(vm_fr[l]*vm_to[l]*sin(va_fr[l]-va_to[l])) -
       p[(l,i,j)]
       for (l,i,j) in ref[:arcs_from]
    ]

    # @NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
    power_flow_p_to_con = [
       (br_g[l]+br_g_to[l])*vm_to[l]^2 +
       (-br_g[l]*br_tr[l]-br_b[l]*br_ti[l])/br_tm[l]*(vm_to[l]*vm_fr[l]*cos(va_to[l]-va_fr[l])) +
       (-br_b[l]*br_tr[l]+br_g[l]*br_ti[l])/br_tm[l]*(vm_to[l]*vm_fr[l]*sin(va_to[l]-va_fr[l])) -
       p[(l,i,j)]
       for (l,i,j) in ref[:arcs_to]
    ]

    # @NLconstraint(model, q_fr == -(b+b_fr)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
    power_flow_q_from_con = [
       -(br_b[l]+br_b_fr[l])/br_tm[l]*vm_fr[l]^2 -
       (-br_b[l]*br_tr[l]-br_g[l]*br_ti[l])/br_tm[l]*(vm_fr[l]*vm_to[l]*cos(va_fr[l]-va_to[l])) +
       (-br_g[l]*br_tr[l]+br_b[l]*br_ti[l])/br_tm[l]*(vm_fr[l]*vm_to[l]*sin(va_fr[l]-va_to[l])) - 
       q[(l,i,j)]
       for (l,i,j) in ref[:arcs_from]
    ]

    # @NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
    power_flow_q_to_con = [
       -(br_b[l]+br_b_to[l])*vm_to[l]^2 -
       (-br_b[l]*br_tr[l]+br_g[l]*br_ti[l])/br_tm[l]*(vm_to[l]*vm_fr[l]*cos(va_fr[l]-va_to[l])) +
       (-br_g[l]*br_tr[l]-br_b[l]*br_ti[l])/br_tm[l]*(vm_to[l]*vm_fr[l]*sin(va_to[l]-va_fr[l])) -
       q[(l,i,j)]
       for (l,i,j) in ref[:arcs_to]
    ]

    # @constraint(model, va_fr - va_to <= branch["angmax"])
    # @constraint(model, va_fr - va_to >= branch["angmin"])
    power_flow_vad_con = [
       va_fr[l] - va_to[l]
       for (l,i,j) in ref[:arcs_from]
    ]

    # @constraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
    power_flow_mva_from_con = [
       p[(l,i,j)]^2 + q[(l,i,j)]^2
       for (l,i,j) in ref[:arcs_from]
    ]

    # @constraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2)
    power_flow_mva_to_con = [
       p[(l,i,j)]^2 + q[(l,i,j)]^2
       for (l,i,j) in ref[:arcs_to]
    ]

    c .= [
        va_con...,
        power_balance_p_con...,
        power_balance_q_con...,
        power_flow_p_from_con...,
        power_flow_p_to_con...,
        power_flow_q_from_con...,
        power_flow_q_to_con...,
        power_flow_vad_con...,
        power_flow_mva_from_con...,
        power_flow_mva_to_con...,
    ]

    return c
end

con_lbs = Float64[]
con_ubs = Float64[]

#@constraint(model, va[i] == 0)
for (i,bus) in ref[:ref_buses]
    push!(con_lbs, 0.0)
    push!(con_ubs, 0.0)
end


#power_balance_p_con
for (i,bus) in ref[:bus]
    push!(con_lbs, 0.0)
    push!(con_ubs, 0.0)
    #push!(con_lbs, -Inf)
    #push!(con_ubs, Inf)
end

#power_balance_q_con
for (i,bus) in ref[:bus]
    push!(con_lbs, 0.0)
    push!(con_ubs, 0.0)
    #push!(con_lbs, -Inf)
    #push!(con_ubs, Inf)
end


#power_flow_p_from_con
for (l,i,j) in ref[:arcs_from]
    push!(con_lbs, 0.0)
    push!(con_ubs, 0.0)
end

#power_flow_p_to_con
for (l,i,j) in ref[:arcs_to]
    push!(con_lbs, 0.0)
    push!(con_ubs, 0.0)
end

#power_flow_q_from_con
for (l,i,j) in ref[:arcs_from]
    push!(con_lbs, 0.0)
    push!(con_ubs, 0.0)
end

#power_flow_q_to_con
for (l,i,j) in ref[:arcs_to]
    push!(con_lbs, 0.0)
    push!(con_ubs, 0.0)
end

#power_flow_vad_con
for (l,i,j) in ref[:arcs_from]
    branch = ref[:branch][l]
    push!(con_lbs, branch["angmin"])
    push!(con_ubs, branch["angmax"])
end

#power_flow_mva_from_con
for (l,i,j) in ref[:arcs_from]
    branch = ref[:branch][l]
    push!(con_lbs, -Inf)
    push!(con_ubs, branch["rate_a"]^2)
end

#power_flow_mva_to_con
for (l,i,j) in ref[:arcs_to]
    branch = ref[:branch][l]
    push!(con_lbs, -Inf)
    push!(con_ubs, branch["rate_a"]^2)
end

println("variables: $(length(var_init)), $(length(var_lb)), $(length(var_ub))")
println("constraints: $(length(opf_constraints(zeros(length(con_lbs)), var_init))), $(length(con_lbs)), $(length(con_ubs))")

df = TwiceDifferentiable(opf_objective, var_init)
dfc = TwiceDifferentiableConstraints(opf_constraints, var_lb, var_ub, con_lbs, con_ubs)

model_build_time = time() - time_start

time_start = time()

options = Optim.Options(show_trace=true)

# NOTE: had to change initial guess to be an interior point, otherwise getting NaN values
res = optimize(df, dfc, var_init, IPNewton(), options)
#res = optimize(df, dfc, var_init, LBFGS(), options) #  StackOverflowError:
#res = optimize(df, dfc, var_init, NelderMead(), options) #  StackOverflowError:
display(res)

sol = res.minimizer
cost = res.minimum

#println(cost)
#println(sol)

# NOTE: confirmed these constraint violations can be eliminated
# if a better starting point is used
sol_eval = opf_constraints(zeros(length(con_lbs)), sol)
vio_lb = [max(v,0) for v in (con_lbs .- sol_eval)]
vio_ub = [max(v,0) for v in (sol_eval .- con_ubs)]
const_vio = vio_lb .+ vio_ub
#println(const_vio)
println("total constraint violation: $(sum(const_vio))")

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

if sum(const_vio) > 1.0
    error("optimize failed to satify the problem constraints")
end

