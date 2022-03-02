time_start = time()

using PowerModels
using Ipopt
using JuMP
using ComplexOptInterface

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

model = Model(Ipopt.Optimizer)
#set_optimizer_attribute(model, "print_level", 0)

ComplexOptInterface.add_all_bridges(model)


@variable(model, -ref[:bus][i]["vmax"] <= vr[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)
@variable(model, -ref[:bus][i]["vmax"] <= vi[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=0.0)
V = Dict(i => vr[i] + vi[i]*im for i in keys(ref[:bus]))

@variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
@variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])
G = Dict(i => pg[i] + qg[i]*im for i in keys(ref[:gen]))

@variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
@variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
S = Dict(k => p[k] + q[k]*im for k in ref[:arcs])


@objective(model, Min, sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]))

for (i,bus) in ref[:ref_buses]
    @constraint(model, vi[i] == 0)
end

for (i,bus) in ref[:bus]
    bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
    bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

    @constraint(model, abs2(V[i]) >= bus["vmin"]^2)
    @constraint(model, abs2(V[i]) <= bus["vmax"]^2)

    @constraint(model,
        sum(S[a] for a in ref[:bus_arcs][i]) ==
        sum(G[g] for g in ref[:bus_gens][i]) -
        sum(load["pd"] + load["qd"]*im for load in bus_loads) -
        sum(shunt["gs"] - shunt["bs"]*im for shunt in bus_shunts)*abs2(V[i])
    )
end

# Branch power flow physics and limit constraints
for (i,branch) in ref[:branch]
    f_idx = (i, branch["f_bus"], branch["t_bus"])
    t_idx = (i, branch["t_bus"], branch["f_bus"])

    S_fr = S[f_idx]
    S_to = S[t_idx]

    V_fr = V[branch["f_bus"]]
    V_to = V[branch["t_bus"]]

    g, b = PowerModels.calc_branch_y(branch)
    Y = g + b*im

    tr, ti = PowerModels.calc_branch_t(branch)
    T = tr + ti*im

    Y_fr = branch["g_fr"] + branch["b_fr"]*im
    Y_to = branch["g_to"] + branch["b_to"]*im

    # Branch Flow
    @constraint(model, S_fr == conj(Y+Y_fr)/(T*conj(T))*V_fr*conj(V_fr) - conj(Y)/T*V_fr*conj(V_to))
    @constraint(model, S_to == conj(Y+Y_to)*V_to*conj(V_to) - conj(Y)/conj(T)*conj(V_fr)*V_to)

    # Voltage angle difference limit
    W = V_fr*conj(V_to)
    @constraint(model, imag(W) <= tan(branch["angmax"])*real(W))
    @constraint(model, imag(W) >= tan(branch["angmin"])*real(W))

    # Apparent power limit, from side and to side
    @constraint(model, abs2(S_fr) <= branch["rate_a"]^2)
    @constraint(model, abs2(S_to) <= branch["rate_a"]^2)
end

model_build_time = time() - time_start

#println("model building done")
#println(model)

time_start = time()

println("try to optimize")
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
