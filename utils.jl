import PowerModels

function create_common_datastructures(file_name)
    time_data_start = time()
    data = PowerModels.parse_file(file_name)
    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    data_load_time = time() - time_data_start
    time_model_start = time()
    bus_pd = Dict(i => 0.0 for (i, _) in ref[:bus])
    bus_qd = Dict(i => 0.0 for (i, _) in ref[:bus])
    bus_gs = Dict(i => 0.0 for (i, _) in ref[:bus])
    bus_bs = Dict(i => 0.0 for (i, _) in ref[:bus])
    for (i, _) in ref[:bus]
        if length(ref[:bus_loads][i]) > 0
            bus_pd[i] = sum(ref[:load][l]["pd"] for l in ref[:bus_loads][i])
            bus_qd[i] = sum(ref[:load][l]["qd"] for l in ref[:bus_loads][i])
        end
        if length(ref[:bus_shunts][i]) > 0
            bus_gs[i] = sum(ref[:shunt][s]["gs"] for s in ref[:bus_shunts][i])
            bus_bs[i] = sum(ref[:shunt][s]["bs"] for s in ref[:bus_shunts][i])
        end
    end
    br_g = Dict(i => 0.0 for (i, _) in ref[:branch])
    br_b = Dict(i => 0.0 for (i, _) in ref[:branch])
    br_tr = Dict(i => 0.0 for (i, _) in ref[:branch])
    br_ti = Dict(i => 0.0 for (i, _) in ref[:branch])
    br_tm = Dict(i => 0.0 for (i, _) in ref[:branch])
    br_g_fr = Dict(i => 0.0 for (i, _) in ref[:branch])
    br_b_fr = Dict(i => 0.0 for (i, _) in ref[:branch])
    br_g_to = Dict(i => 0.0 for (i, _) in ref[:branch])
    br_b_to = Dict(i => 0.0 for (i, _) in ref[:branch])
    for (i, branch) in ref[:branch]
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
    for (i, bus) in ref[:bus]
        push!(var_init, 0.0) #va
        push!(var_lb, -Inf)
        push!(var_ub, Inf)
        var_lookup["va_$(i)"] = var_idx
        var_idx += 1
        push!(var_init, 1.0) #vm
        push!(var_lb, bus["vmin"])
        push!(var_ub, bus["vmax"])
        var_lookup["vm_$(i)"] = var_idx
        var_idx += 1
    end
    for (i, gen) in ref[:gen]
        push!(var_init, 0.0) #pg
        push!(var_lb, gen["pmin"])
        push!(var_ub, gen["pmax"])
        var_lookup["pg_$(i)"] = var_idx
        var_idx += 1
        push!(var_init, 0.0) #qg
        push!(var_lb, gen["qmin"])
        push!(var_ub, gen["qmax"])
        var_lookup["qg_$(i)"] = var_idx
        var_idx += 1
    end
    for (l, i, j) in ref[:arcs]
        branch = ref[:branch][l]
        push!(var_init, 0.0) #p
        push!(var_lb, -branch["rate_a"])
        push!(var_ub,  branch["rate_a"])
        var_lookup["p_$(l)_$(i)_$(j)"] = var_idx
        var_idx += 1
        push!(var_init, 0.0) #q
        push!(var_lb, -branch["rate_a"])
        push!(var_ub,  branch["rate_a"])
        var_lookup["q_$(l)_$(i)_$(j)"] = var_idx
        var_idx += 1
    end
    @assert var_idx == length(var_init)+1
    # param argument is for GalacticOptim.
    function opf_objective(x, param = nothing)
        cost = 0.0
        for (i, gen) in ref[:gen]
            pg = x[var_lookup["pg_$(i)"]]
            cost += gen["cost"][1]*pg^2 + gen["cost"][2]*pg + gen["cost"][3]
        end
        return cost
    end

    # param argument is for GalacticOptim.
    function opf_constraints(x, param = nothing)
        va = Dict(i => x[var_lookup["va_$(i)"]] for (i, _) in ref[:bus])
        vm = Dict(i => x[var_lookup["vm_$(i)"]] for (i, _) in ref[:bus])
        pg = Dict(i => x[var_lookup["pg_$(i)"]] for (i, _) in ref[:gen])
        qg = Dict(i => x[var_lookup["qg_$(i)"]] for (i, _) in ref[:gen])
        p = Dict(
            (l, i, j) => x[var_lookup["p_$(l)_$(i)_$(j)"]] for
            (l, i, j) in ref[:arcs]
        )
        q = Dict(
            (l, i, j) => x[var_lookup["q_$(l)_$(i)_$(j)"]] for
            (l, i, j) in ref[:arcs]
        )
        vm_fr = Dict(l => vm[branch["f_bus"]] for (l, branch) in ref[:branch])
        vm_to = Dict(l => vm[branch["t_bus"]] for (l, branch) in ref[:branch])
        va_fr = Dict(l => va[branch["f_bus"]] for (l, branch) in ref[:branch])
        va_to = Dict(l => va[branch["t_bus"]] for (l, branch) in ref[:branch])
        va_con = [va[i] for (i, _) in ref[:ref_buses]]
        # @constraint(model,
        #     sum(p[a] for a in ref[:bus_arcs][i]) ==
        #     sum(pg[g] for g in ref[:bus_gens][i]) -
        #     sum(load["pd"] for load in bus_loads) -
        #     sum(shunt["gs"] for shunt in bus_shunts)*vm[i]^2
        # )
        power_balance_p_con = [
           sum(pg[j] for j in ref[:bus_gens][i]; init = 0.0) -
           bus_pd[i] -
           bus_gs[i] * vm[i]^2 -
           sum(p[a] for a in ref[:bus_arcs][i])
           for (i, _) in ref[:bus]
        ]
        # @constraint(model,
        #     sum(q[a] for a in ref[:bus_arcs][i]) ==
        #     sum(qg[g] for g in ref[:bus_gens][i]) -
        #     sum(load["qd"] for load in bus_loads) +
        #     sum(shunt["bs"] for shunt in bus_shunts)*vm[i]^2
        # )
        power_balance_q_con = [
           sum(qg[j] for j in ref[:bus_gens][i]; init = 0.0) -
           bus_qd[i] +
           bus_bs[i] * vm[i]^2 -
           sum(q[a] for a in ref[:bus_arcs][i])
           for (i, _) in ref[:bus]
        ]
        # @NLconstraint(model, p_fr ==  (g+g_fr)/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        power_flow_p_from_con = [
           (br_g[l] + br_g_fr[l]) / br_tm[l] * vm_fr[l]^2 +
           (-br_g[l] * br_tr[l] + br_b[l] * br_ti[l]) / br_tm[l] * (vm_fr[l] * vm_to[l] * cos(va_fr[l] - va_to[l])) +
           (-br_b[l] * br_tr[l] - br_g[l] * br_ti[l]) / br_tm[l] * (vm_fr[l] * vm_to[l] * sin(va_fr[l] - va_to[l])) -
           p[(l, i, j)]
           for (l, i, j) in ref[:arcs_from]
        ]
        # @NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        power_flow_p_to_con = [
           (br_g[l]+br_g_to[l])*vm_to[l]^2 +
           (-br_g[l]*br_tr[l]-br_b[l]*br_ti[l])/br_tm[l]*(vm_to[l]*vm_fr[l]*cos(va_to[l]-va_fr[l])) +
           (-br_b[l]*br_tr[l]+br_g[l]*br_ti[l])/br_tm[l]*(vm_to[l]*vm_fr[l]*sin(va_to[l]-va_fr[l])) -
           p[(l, i, j)]
           for (l, i, j) in ref[:arcs_to]
        ]
        # @NLconstraint(model, q_fr == -(b+b_fr)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        power_flow_q_from_con = [
           -(br_b[l]+br_b_fr[l])/br_tm[l]*vm_fr[l]^2 -
           (-br_b[l]*br_tr[l]-br_g[l]*br_ti[l])/br_tm[l]*(vm_fr[l]*vm_to[l]*cos(va_fr[l]-va_to[l])) +
           (-br_g[l]*br_tr[l]+br_b[l]*br_ti[l])/br_tm[l]*(vm_fr[l]*vm_to[l]*sin(va_fr[l]-va_to[l])) -
           q[(l, i, j)]
           for (l, i, j) in ref[:arcs_from]
        ]
        # @NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        power_flow_q_to_con = [
           -(br_b[l] + br_b_to[l]) * vm_to[l]^2 -
           (-br_b[l] * br_tr[l] + br_g[l] * br_ti[l]) / br_tm[l] * (vm_to[l] * vm_fr[l] * cos(va_fr[l] - va_to[l])) +
           (-br_g[l] * br_tr[l] - br_b[l] * br_ti[l]) / br_tm[l] * (vm_to[l] * vm_fr[l] * sin(va_to[l] - va_fr[l])) -
           q[(l, i, j)]
           for (l, i, j) in ref[:arcs_to]
        ]
        # @constraint(model, va_fr - va_to <= branch["angmax"])
        # @constraint(model, va_fr - va_to >= branch["angmin"])
        power_flow_vad_con = [
           va_fr[l] - va_to[l]
           for (l, _, _) in ref[:arcs_from]
        ]
        # @constraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
        power_flow_mva_from_con = [
           p[(l, i, j)]^2 + q[(l, i, j)]^2
           for (l, i, j) in ref[:arcs_from]
        ]
        # @constraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2)
        power_flow_mva_to_con = [
           p[(l, i, j)]^2 + q[(l, i, j)]^2
           for (l, i, j) in ref[:arcs_to]
        ]
        return [
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
    end
    con_lbs = Float64[]
    con_ubs = Float64[]
    #@constraint(model, va[i] == 0)
    for (i, bus) in ref[:ref_buses]
        push!(con_lbs, 0.0)
        push!(con_ubs, 0.0)
    end
    #power_balance_p_con
    for (i, bus) in ref[:bus]
        push!(con_lbs, 0.0)
        push!(con_ubs, 0.0)
    end
    #power_balance_q_con
    for (i, bus) in ref[:bus]
        push!(con_lbs, 0.0)
        push!(con_ubs, 0.0)
    end
    #power_flow_p_from_con
    for (l, i, j) in ref[:arcs_from]
        push!(con_lbs, 0.0)
        push!(con_ubs, 0.0)
    end
    #power_flow_p_to_con
    for (l, i, j) in ref[:arcs_to]
        push!(con_lbs, 0.0)
        push!(con_ubs, 0.0)
    end
    #power_flow_q_from_con
    for (l, i, j) in ref[:arcs_from]
        push!(con_lbs, 0.0)
        push!(con_ubs, 0.0)
    end
    #power_flow_q_to_con
    for (l, i, j) in ref[:arcs_to]
        push!(con_lbs, 0.0)
        push!(con_ubs, 0.0)
    end
    #power_flow_vad_con
    for (l, i, j) in ref[:arcs_from]
        branch = ref[:branch][l]
        push!(con_lbs, branch["angmin"])
        push!(con_ubs, branch["angmax"])
    end
    #power_flow_mva_from_con
    for (l, i, j) in ref[:arcs_from]
        branch = ref[:branch][l]
        push!(con_lbs, -Inf)
        push!(con_ubs, branch["rate_a"]^2)
    end
    #power_flow_mva_to_con
    for (l, i, j) in ref[:arcs_to]
        branch = ref[:branch][l]
        push!(con_lbs, -Inf)
        push!(con_ubs, branch["rate_a"]^2)
    end
    return (
        var_init = var_init,
        var_lb = var_lb,
        var_ub = var_ub,
        con_lb = con_lbs,
        con_ub = con_ubs,
        cons = opf_constraints,
        obj = opf_objective,
        ref = ref,
        data_load_time = data_load_time,
        time_data_start = time_data_start,
        time_model_start = time_model_start,
    )
end

function solve_opf(build_model, solve, file_name)
    data = create_common_datastructures(file_name)
    model = build_model(data)
    model_build_time = time() - data.time_model_start
    time_solve_start = time()
    solution = solve(model)
    solve_time = time() - time_solve_start
    total_time = time() - data.time_data_start
    println("")
    println("\033[1mSummary\033[0m")
    println("   case........: $(file_name)")
    println("   feasible....: $(solution.feasible)")
    println("   cost........: $(round(Int, solution.cost))")
    println("   total time..: $(total_time)")
    println("     data time.: $(data.data_load_time)")
    println("     build time: $(model_build_time)")
    println("     solve time: $(solve_time)")
    println("")
    return Dict(
        "case" => file_name,
        "feasible" => solution.feasible,
        "cost" => solution.cost,
        "time_total" => total_time,
        "time_data" => data.data_load_time,
        "time_build" => model_build_time,
        "time_solve" => solve_time,
    )
end
