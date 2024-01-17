#!/usr/bin/env julia
###### AC-OPF using JuMP ######
#
# implementation reference: https://github.com/lanl-ansi/PowerModelsAnnex.jl/blob/master/src/model/ac-opf.jl
#

import PowerModels
import Ipopt
import JuMP

function solve_opf(file_name)
    time_data_start = time()

    data = PowerModels.parse_file(file_name)
    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    data_load_time = time() - time_data_start



    time_model_start = time()

    model = JuMP.Model(Ipopt.Optimizer)
    #set_optimizer_attribute(model, "print_level", 0)

    JuMP.@variable(model, va[i in keys(ref[:bus])])
    JuMP.@variable(model, -ref[:bus][i]["vmax"] <= vr[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)
    JuMP.@variable(model, -ref[:bus][i]["vmax"] <= vi[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=0.0)

    JuMP.@variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    JuMP.@variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

    JuMP.@objective(model, Min, sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]))

    for (i,bus) in ref[:ref_buses]
        JuMP.@constraint(model, vi[i] == 0)
    end

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        JuMP.@constraint(model, bus["vmin"]^2 <= vr[i]^2 + vi[i]^2)
        JuMP.@constraint(model, bus["vmax"]^2 >= vr[i]^2 + vi[i]^2)

        JuMP.@constraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2 + vi[i]^2)
        )

        JuMP.@constraint(model,
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
        ttm = tr^2 + ti^2
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]

        # From side of the branch flow
        JuMP.@constraint(model, p_fr ==  (g+g_fr)/ttm*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/ttm*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/ttm*(vi_fr*vr_to - vr_fr*vi_to) )
        JuMP.@constraint(model, q_fr == -(b+b_fr)/ttm*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/ttm*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/ttm*(vi_fr*vr_to - vr_fr*vi_to) )

        # To side of the branch flow
        JuMP.@constraint(model, p_to ==  (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/ttm*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/ttm*(-(vi_fr*vr_to - vr_fr*vi_to)) )
        JuMP.@constraint(model, q_to == -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/ttm*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/ttm*(-(vi_fr*vr_to - vr_fr*vi_to)) )

        # Voltage angle difference limit
        JuMP.@constraint(model, (vi_fr*vr_to - vr_fr*vi_to) <= tan(branch["angmax"])*(vr_fr*vr_to + vi_fr*vi_to))
        JuMP.@constraint(model, (vi_fr*vr_to - vr_fr*vi_to) >= tan(branch["angmin"])*(vr_fr*vr_to + vi_fr*vi_to))

        # Apparent power limit, from side and to side
        JuMP.@constraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
        JuMP.@constraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2)
    end

    model_variables = JuMP.num_variables(model)

    # for consistency with other solvers, skip the variable bounds in the constraint count
    non_nl_constraints = sum(JuMP.num_constraints(model, ft, st) for (ft, st) in JuMP.list_of_constraint_types(model) if ft != JuMP.VariableRef)
    model_constraints = JuMP.num_nonlinear_constraints(model) + non_nl_constraints

    model_build_time = time() - time_model_start

    time_solve_start = time()

    JuMP.optimize!(model)
    cost = JuMP.objective_value(model)
    feasible = (JuMP.termination_status(model) == JuMP.LOCALLY_SOLVED)

    solve_time = time() - time_solve_start
    total_time = time() - time_data_start

    println("")
    println("\033[1mSummary\033[0m")
    println("   case........: $(file_name)")
    println("   variables...: $(model_variables)")
    println("   constraints.: $(model_constraints)")
    println("   feasible....: $(feasible)")
    println("   cost........: $(round(Int, cost))")
    println("   total time..: $(total_time)")
    println("     data time.: $(data_load_time)")
    println("     build time: $(model_build_time)")
    println("     solve time: $(solve_time)")
    println("")

    return Dict(
        "case" => file_name,
        "variables" => model_variables,
        "constraints" => model_constraints,
        "feasible" => feasible,
        "cost" => cost,
        "time_total" => total_time,
        "time_data" => data_load_time,
        "time_build" => model_build_time,
        "time_solve" => solve_time,
    )
end

if isinteractive() == false
    solve_opf("$(@__DIR__)/../data/opf_warmup.m")
end