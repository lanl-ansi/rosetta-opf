#!/usr/bin/env julia
###### AC-OPF using JuMP ######
#
# implementation reference: https://github.com/lanl-ansi/PowerModelsAnnex.jl/blob/master/src/model/ac-opf.jl
# requires JuMP v0.23.1, ComplexOptInterface v0.1.1
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

    JuMP.@variable(model,
        V[i in keys(ref[:bus])] in JuMP.ComplexPlane(),
        lower_bound = -ref[:bus][i]["vmax"] - ref[:bus][i]["vmax"] * im,
        upper_bound =  ref[:bus][i]["vmax"] + ref[:bus][i]["vmax"] * im,
        start=1.0 + 0.0im
    )

    JuMP.@variable(model,
        G[i in keys(ref[:gen])] in JuMP.ComplexPlane(),
        lower_bound = ref[:gen][i]["pmin"] + ref[:gen][i]["qmin"] * im,
        upper_bound = ref[:gen][i]["pmax"] + ref[:gen][i]["qmax"] * im,
    )

    JuMP.@variable(model,
        S[(l,i,j) in ref[:arcs]] in JuMP.ComplexPlane(),
        lower_bound = -ref[:branch][l]["rate_a"] - ref[:branch][l]["rate_a"] * im,
        upper_bound =  ref[:branch][l]["rate_a"] + ref[:branch][l]["rate_a"] * im,
    )

    JuMP.@objective(model, Min, sum(gen["cost"][1]*real(G[i])^2 + gen["cost"][2]*real(G[i]) + gen["cost"][3] for (i,gen) in ref[:gen]))

    for (i,bus) in ref[:ref_buses]
        JuMP.@constraint(model, imag(V[i]) == 0)
    end

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        JuMP.@constraint(model, abs2(V[i]) >= bus["vmin"]^2)
        JuMP.@constraint(model, abs2(V[i]) <= bus["vmax"]^2)

        JuMP.@constraint(model,
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
        JuMP.@constraint(model, S_fr == conj(Y+Y_fr)/(T*conj(T))*V_fr*conj(V_fr) - conj(Y)/T*V_fr*conj(V_to))
        JuMP.@constraint(model, S_to == conj(Y+Y_to)*V_to*conj(V_to) - conj(Y)/conj(T)*conj(V_fr)*V_to)

        # Voltage angle difference limit
        W = V_fr*conj(V_to)
        JuMP.@constraint(model, imag(W) <= tan(branch["angmax"])*real(W))
        JuMP.@constraint(model, imag(W) >= tan(branch["angmin"])*real(W))

        # Apparent power limit, from side and to side
        JuMP.@constraint(model, abs2(S_fr) <= branch["rate_a"]^2)
        JuMP.@constraint(model, abs2(S_to) <= branch["rate_a"]^2)
    end

    model_variables = JuMP.num_variables(model)

    # for consistency with other solvers, skip the variable bounds in the constraint count
    non_nl_constraints = sum(JuMP.num_constraints(model, ft, st) for (ft, st) in JuMP.list_of_constraint_types(model) if ft != JuMP.VariableRef)
    model_constraints = JuMP.num_nonlinear_constraints(model) + non_nl_constraints

    model_build_time = time() - time_model_start

    #println("model building done")
    #println(model)

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
    solve_opf("$(@__DIR__)/../data/pglib_opf_case5_pjm.m")
end
