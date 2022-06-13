#!/usr/bin/env julia
###### AC-OPF using Nonconvex ######
#
# implementation reference: https://julianonconvex.github.io/Nonconvex.jl/stable/problem/
# currently does not converge due to an upstream issue with the AD backend Zygote: https://github.com/JuliaNonconvex/Nonconvex.jl/issues/130
#


import PowerModels
import Nonconvex
Nonconvex.@load Ipopt


function solve_opf(file_name)
    time_data_start = time()

    data = PowerModels.parse_file(file_name)
    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    data_load_time = time() - time_data_start


    time_model_start = time()

    model = Nonconvex.DictModel()

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
    br_ttm = Dict(i => 0.0 for (i,branch) in ref[:branch])

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
        br_ttm[i] = tr^2 + ti^2

        br_g_fr[i] = branch["g_fr"]
        br_b_fr[i] = branch["b_fr"]
        br_g_to[i] = branch["g_to"]
        br_b_to[i] = branch["b_to"]
    end

    for (i,bus) in ref[:bus]
        addvar!(model, "va_$(i)", -Inf, Inf) #va
        addvar!(model, "vm_$(i)", bus["vmin"], bus["vmax"], init=1.0) #vm
    end

    for (i,gen) in ref[:gen]
        addvar!(model, "pg_$(i)", gen["pmin"], gen["pmax"]) #pg
        addvar!(model, "qg_$(i)", gen["qmin"], gen["qmax"]) #qg
    end

    for (l,i,j) in ref[:arcs]
        branch = ref[:branch][l]
        addvar!(model, "p_$(l)_$(i)_$(j)", -branch["rate_a"], branch["rate_a"]) #p
        addvar!(model, "q_$(l)_$(i)_$(j)", -branch["rate_a"], branch["rate_a"]) #q
    end


    #total_callback_time = 0.0
    function opf_objective(x::OrderedDict)
        #start = time()
        cost = 0.0
        for (i,gen) in ref[:gen]
            pg = x["pg_$(i)"]
            cost += gen["cost"][1]*pg^2 + gen["cost"][2]*pg + gen["cost"][3]
        end
        #total_callback_time += time() - start
        return cost
    end
    Nonconvex.set_objective!(model, opf_objective)


    function const_ref_bus(x::OrderedDict, i)
        return x["va_$(i)"]
    end
    for (i,bus) in ref[:ref_buses]
        add_eq_constraint!(model, x -> const_ref_bus(x,i))
    end


    # TBD once AD is resolved
    #add_ineq_constraint!(model, f)
    #add_eq_constraint!(model, f)

    #=
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


    # @NLconstraint(model, p_fr ==  (g+g_fr)/ttm*vm_fr^2 + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )
    power_flow_p_from_con = [
       (br_g[l]+br_g_fr[l])/br_ttm[l]*vm_fr[l]^2 +
       (-br_g[l]*br_tr[l]+br_b[l]*br_ti[l])/br_ttm[l]*(vm_fr[l]*vm_to[l]*cos(va_fr[l]-va_to[l])) +
       (-br_b[l]*br_tr[l]-br_g[l]*br_ti[l])/br_ttm[l]*(vm_fr[l]*vm_to[l]*sin(va_fr[l]-va_to[l])) -
       p[(l,i,j)]
       for (l,i,j) in ref[:arcs_from]
    ]

    # @NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )
    power_flow_p_to_con = [
       (br_g[l]+br_g_to[l])*vm_to[l]^2 +
       (-br_g[l]*br_tr[l]-br_b[l]*br_ti[l])/br_ttm[l]*(vm_to[l]*vm_fr[l]*cos(va_to[l]-va_fr[l])) +
       (-br_b[l]*br_tr[l]+br_g[l]*br_ti[l])/br_ttm[l]*(vm_to[l]*vm_fr[l]*sin(va_to[l]-va_fr[l])) -
       p[(l,i,j)]
       for (l,i,j) in ref[:arcs_to]
    ]

    # @NLconstraint(model, q_fr == -(b+b_fr)/ttm*vm_fr^2 - (-b*tr-g*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )
    power_flow_q_from_con = [
       -(br_b[l]+br_b_fr[l])/br_ttm[l]*vm_fr[l]^2 -
       (-br_b[l]*br_tr[l]-br_g[l]*br_ti[l])/br_ttm[l]*(vm_fr[l]*vm_to[l]*cos(va_fr[l]-va_to[l])) +
       (-br_g[l]*br_tr[l]+br_b[l]*br_ti[l])/br_ttm[l]*(vm_fr[l]*vm_to[l]*sin(va_fr[l]-va_to[l])) -
       q[(l,i,j)]
       for (l,i,j) in ref[:arcs_from]
    ]

    # @NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )
    power_flow_q_to_con = [
       -(br_b[l]+br_b_to[l])*vm_to[l]^2 -
       (-br_b[l]*br_tr[l]+br_g[l]*br_ti[l])/br_ttm[l]*(vm_to[l]*vm_fr[l]*cos(va_to[l]-va_fr[l])) +
       (-br_g[l]*br_tr[l]-br_b[l]*br_ti[l])/br_ttm[l]*(vm_to[l]*vm_fr[l]*sin(va_to[l]-va_fr[l])) -
       q[(l,i,j)]
       for (l,i,j) in ref[:arcs_to]
    ]
    =#

    model_variables = Nonconvex.NonconvexCore.getnvars(model)
    model_constraints = Nonconvex.NonconvexCore.getnconstraints(model)


    model_build_time = time() - time_model_start


    time_solve_start = time()

    x0 = NonconvexCore.getinit(model)
    first_order = true
    options = IpoptOptions(; first_order)
    sym_model = Nonconvex.symbolify(model, hessian=!first_order, sparse=true, simplify=true)
    result = Nonconvex.optimize(sym_model, IpoptAlg(), x0; options)

    println(result.minimizer)
    cost = result.minimum
    feasible = result.status == 0 # just guessing this is correct

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
    solve_opf("$(@__DIR__)/data/pglib_opf_case5_pjm.m")
end
