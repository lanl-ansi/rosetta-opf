time_start = time()

using PowerModels
using Ipopt
using JuMP

pkg_load_time = time() - time_start


time_start = time()

file_name = "../data/pglib_opf_case5_pjm.m"
#file_name = "data/pglib_opf_case118_ieee.m"

data = PowerModels.parse_file(file_name)
PowerModels.standardize_cost_terms!(data, order=2)
PowerModels.calc_thermal_limits!(data)
ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

data_load_time = time() - time_start



time_start = time()

model = Model(Ipopt.Optimizer)
#set_optimizer_attribute(model, "print_level", 0)

@variable(model, va[i in keys(ref[:bus])])
@variable(model, -ref[:bus][i]["vmax"] <= vr[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)
@variable(model, -ref[:bus][i]["vmax"] <= vi[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)

@variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
@variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

@variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
@variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

@objective(model, Min, sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]))

for (i,bus) in ref[:ref_buses]
    @constraint(model, vi[i] == 0)
end

for (i,bus) in ref[:bus]
    bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
    bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

    @constraint(model, bus["vmin"]^2 <= vr[i]^2 + vi[i]^2)
    @constraint(model, bus["vmax"]^2 >= vr[i]^2 + vi[i]^2)

    @constraint(model,
        sum(p[a] for a in ref[:bus_arcs][i]) ==
        sum(pg[g] for g in ref[:bus_gens][i]) -
        sum(load["pd"] for load in bus_loads) -
        sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2 + vi[i]^2)
    )

    @constraint(model,
        sum(q[a] for a in ref[:bus_arcs][i]) ==
        sum(qg[g] for g in ref[:bus_gens][i]) -
        sum(load["qd"] for load in bus_loads) +
        sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2 + vi[i]^2)
    )

end

# Branch power flow physics and limit constraints
for (i,branch) in ref[:branch]
    f_idx = (i, branch["f_bus"], branch["t_bus"])
    t_idx = (i, branch["t_bus"], branch["f_bus"])

    p_fr = p[f_idx]                     # p_fr is a reference to the optimization variable p[f_idx]
    q_fr = q[f_idx]                     # q_fr is a reference to the optimization variable q[f_idx]
    p_to = p[t_idx]                     # p_to is a reference to the optimization variable p[t_idx]
    q_to = q[t_idx]                     # q_to is a reference to the optimization variable q[t_idx]

    vr_fr = vr[branch["f_bus"]]         # vm_fr is a reference to the optimization variable vm on the from side of the branch
    vr_to = vr[branch["t_bus"]]         # vm_to is a reference to the optimization variable vm on the to side of the branch
    vi_fr = vi[branch["f_bus"]]         # va_fr is a reference to the optimization variable va on the from side of the branch
    vi_to = vi[branch["t_bus"]]         # va_fr is a reference to the optimization variable va on the to side of the branch

    g, b = PowerModels.calc_branch_y(branch)
    tr, ti = PowerModels.calc_branch_t(branch)
    tm = tr^2 + ti^2
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    g_to = branch["g_to"]
    b_to = branch["b_to"]

    # From side of the branch flow
    @constraint(model, p_fr ==  (g+g_fr)/tm*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm*(vi_fr*vr_to - vr_fr*vi_to) )
    @constraint(model, q_fr == -(b+b_fr)/tm*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm*(vi_fr*vr_to - vr_fr*vi_to) )

    # To side of the branch flow
    @constraint(model, p_to ==  (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm*(-(vi_fr*vr_to - vr_fr*vi_to)) )
    @constraint(model, q_to == -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm*(-(vi_fr*vr_to - vr_fr*vi_to)) )

    # Voltage angle difference limit
    @constraint(model, (vi_fr*vr_to - vr_fr*vi_to) <= tan(branch["angmax"])*(vr_fr*vr_to + vi_fr*vi_to))
    @constraint(model, (vi_fr*vr_to - vr_fr*vi_to) >= tan(branch["angmin"])*(vr_fr*vr_to + vi_fr*vi_to))

    # Apparent power limit, from side and to side
    @constraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
    @constraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2)
end

model_build_time = time() - time_start


time_start = time()

optimize!(model)
cost = objective_value(model)

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
