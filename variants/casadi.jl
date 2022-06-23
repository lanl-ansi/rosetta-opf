#!/usr/bin/env julia
###### AC-OPF using CasADi via PythonCall ######
#
# implementation reference: https://github.com/lanl-ansi/PowerModelsAnnex.jl/blob/master/src/model/ac-opf.jl
#

import PowerModels
import PythonCall

function solve_opf(file_name)
    time_data_start = time()

    data = PowerModels.parse_file(file_name)
    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    data_load_time = time() - time_data_start

    time_model_start = time()

    casadi = PythonCall.pyimport("casadi")

    x, x0, lbx, ubx, cons, lbg, ubg = [], [], [], [], [], [], []
    va, vm = Dict{Int,Any}(), Dict{Int,Any}()
    for (k, _) in ref[:bus]
        va[k] = casadi.SX.sym("va$k")
        push!(x, va[k])
        push!(x0, 0.0)
        push!(lbx, -casadi.inf)
        push!(ubx, casadi.inf)
        vm[k] = casadi.SX.sym("vm$k")
        push!(x, vm[k])
        push!(x0, 1.0)
        push!(lbx, ref[:bus][k]["vmin"])
        push!(ubx, ref[:bus][k]["vmax"])
    end
    pg, qg = Dict{Int,Any}(), Dict{Int,Any}()
    for (k, ) in ref[:gen]
        pg[k] = casadi.SX.sym("pg$k")
        push!(x, pg[k])
        push!(x0, 0.0)
        push!(lbx, ref[:gen][k]["pmin"])
        push!(ubx, ref[:gen][k]["pmax"])
        qg[k] = casadi.SX.sym("qg$k")
        push!(x, qg[k])
        push!(x0, 0.0)
        push!(lbx, ref[:gen][k]["qmin"])
        push!(ubx, ref[:gen][k]["qmax"])
    end
    p, q = Dict{NTuple{3,Int},Any}(), Dict{NTuple{3,Int},Any}()
    for k in ref[:arcs]
        a = ref[:branch][k[1]]["rate_a"]
        p[k] = casadi.SX.sym("p$k")
        push!(x, p[k])
        push!(x0, 0.0)
        push!(lbx, -a)
        push!(ubx, a)
        q[k] = casadi.SX.sym("q$k")
        push!(x, q[k])
        push!(x0, 0.0)
        push!(lbx, -a)
        push!(ubx, a)
    end
    f = sum(
        cost["cost"][1] * pg[k]^2 +
        cost["cost"][2] * pg[k] +
        cost["cost"][3] for (k, cost) in ref[:gen]
    )
    for (k, _) in ref[:ref_buses]
        push!(cons, va[k])
        push!(lbg, 0)
        push!(ubg, 0)
    end
    for (i, _) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        push!(
            cons,
            sum(p[k] for k in ref[:bus_arcs][i]) -
            sum(pg[g] for g in ref[:bus_gens][i]; init = 0) +
            sum(load["pd"] for load in bus_loads; init = 0) +
            sum(shunt["gs"] for shunt in bus_shunts; init = 0) * vm[i]^2
        )
        push!(lbg, 0)
        push!(ubg, 0)
        push!(
            cons,
            sum(q[k] for k in ref[:bus_arcs][i]) -
            sum(qg[g] for g in ref[:bus_gens][i]; init = 0) +
            sum(load["qd"] for load in bus_loads; init = 0) -
            sum(shunt["bs"] for shunt in bus_shunts; init = 0) * vm[i]^2
        )
        push!(lbg, 0)
        push!(ubg, 0)
    end
    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])
        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]
        vm_fr = vm[branch["f_bus"]]
        vm_to = vm[branch["t_bus"]]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        ttm = tr^2 + ti^2
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        push!(
            cons,
            (g+g_fr)/ttm*vm_fr^2 + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*casadi.cos(va_fr-va_to)) + (-b*tr-g*ti)/ttm*(vm_fr*vm_to*casadi.sin(va_fr-va_to)) - p_fr
        )
        push!(
            cons,
            -(b+b_fr)/ttm*vm_fr^2 - (-b*tr-g*ti)/ttm*(vm_fr*vm_to*casadi.cos(va_fr-va_to)) + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*casadi.sin(va_fr-va_to)) - q_fr
        )
        push!(
            cons,
            (g+g_to)*vm_to^2 + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*casadi.cos(va_to-va_fr)) + (-b*tr+g*ti)/ttm*(vm_to*vm_fr*casadi.sin(va_to-va_fr)) - p_to
        )
        push!(
            cons,
            -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/ttm*(vm_to*vm_fr*casadi.cos(va_to-va_fr)) + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*casadi.sin(va_to-va_fr)) - q_to
        )
        for i in 1:4
            push!(lbg, 0)
            push!(ubg, 0)
        end
        push!(cons, va_fr - va_to)
        push!(lbg, branch["angmin"])
        push!(ubg, branch["angmax"])
        push!(cons, p_fr^2 + q_fr^2)
        push!(lbg, -casadi.inf)
        push!(ubg, branch["rate_a"]^2)
        push!(cons, p_to^2 + q_to^2)
        push!(lbg, -casadi.inf)
        push!(ubg, branch["rate_a"]^2)
    end

    nlp = Dict("x" => casadi.vcat(x), "f" => f, "g" => casadi.vcat(cons))
    options = PythonCall.pydict(Dict("error_on_fail" => true, "ipopt" => PythonCall.pydict(Dict("linear_solver" => "ma27"))))
    model = casadi.nlpsol("model", "ipopt", PythonCall.pydict(nlp), options)

    model_variables = length(x)
    model_constraints = length(lbg)

    model_build_time = time() - time_model_start

    time_solve_start = time()
    solution = model(; lbx = lbx, ubx = ubx, lbg = lbg, ubg = ubg, x0 = x0)
    cost = PythonCall.pyconvert(Float64, (PythonCall.pyfloat(solution["f"])))
    feasible = true # error if not feasible

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
    # println("      callbacks: $(total_callback_time)")
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
        #"time_callbacks" => TBD,
    )
end

if isinteractive() == false
    solve_opf("$(@__DIR__)/../data/pglib_opf_case5_pjm.m")
end
