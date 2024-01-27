#!/usr/bin/env julia
###### AC-OPF using ExaModels ######
#
# implementation reference: https://exanauts.github.io/ExaModels.jl/stable/guide/
# only the built-in AD library is supported
#

import PowerModels
import ExaModels
import NLPModelsIpopt
import LinearAlgebra

function solve_opf(file_name)
    time_data_start = time()

    data = PowerModels.parse_file(file_name)
    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    
    arcdict = Dict(a => k for (k, a) in enumerate(ref[:arcs]))
    busdict = Dict(k => i for (i, (k, v)) in enumerate(ref[:bus]))
    gendict = Dict(k => i for (i, (k, v)) in enumerate(ref[:gen]))
    branchdict = Dict(k => i for (i, (k, v)) in enumerate(ref[:branch]))

    data = (
        bus = [
            begin
                bus_loads = [ref[:load][l] for l in ref[:bus_loads][k]]
                bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][k]]
                pd = sum(load["pd"] for load in bus_loads; init = 0.0)
                gs = sum(shunt["gs"] for shunt in bus_shunts; init = 0.0)
                qd = sum(load["qd"] for load in bus_loads; init = 0.0)
                bs = sum(shunt["bs"] for shunt in bus_shunts; init = 0.0)
                (i = busdict[k], pd = pd, gs = gs, qd = qd, bs = bs)
            end for (k, v) in ref[:bus]
        ],
        gen = [
            (
                i = gendict[k],
                cost1 = v["cost"][1],
                cost2 = v["cost"][2],
                cost3 = v["cost"][3],
                bus = busdict[v["gen_bus"]],
            ) for (k, v) in ref[:gen]
        ],
        arc = [
            (i = k, rate_a = ref[:branch][l]["rate_a"], bus = busdict[i]) for
            (k, (l, i, j)) in enumerate(ref[:arcs])
        ],
        branch = [
            begin
                f_idx = arcdict[i, branch["f_bus"], branch["t_bus"]]
                t_idx = arcdict[i, branch["t_bus"], branch["f_bus"]]
                g, b = PowerModels.calc_branch_y(branch)
                tr, ti = PowerModels.calc_branch_t(branch)
                ttm = tr^2 + ti^2
                g_fr = branch["g_fr"]
                b_fr = branch["b_fr"]
                g_to = branch["g_to"]
                b_to = branch["b_to"]
                c1 = (-g * tr - b * ti) / ttm
                c2 = (-b * tr + g * ti) / ttm
                c3 = (-g * tr + b * ti) / ttm
                c4 = (-b * tr - g * ti) / ttm
                c5 = (g + g_fr) / ttm
                c6 = (b + b_fr) / ttm
                c7 = (g + g_to)
                c8 = (b + b_to)
                (
                    i = branchdict[i],
                    j = 1,
                    f_idx = f_idx,
                    t_idx = t_idx,
                    f_bus = busdict[branch["f_bus"]],
                    t_bus = busdict[branch["t_bus"]],
                    c1 = c1,
                    c2 = c2,
                    c3 = c3,
                    c4 = c4,
                    c5 = c5,
                    c6 = c6,
                    c7 = c7,
                    c8 = c8,
                    rate_a_sq = branch["rate_a"]^2,
                )
            end for (i, branch) in ref[:branch]
        ],
        ref_buses = [busdict[i] for (i, k) in ref[:ref_buses]],
        vmax = [v["vmax"] for (k, v) in ref[:bus]],
        vmin = [v["vmin"] for (k, v) in ref[:bus]],
        pmax = [v["pmax"] for (k, v) in ref[:gen]],
        pmin = [v["pmin"] for (k, v) in ref[:gen]],
        qmax = [v["qmax"] for (k, v) in ref[:gen]],
        qmin = [v["qmin"] for (k, v) in ref[:gen]],
        rate_a = [ref[:branch][l]["rate_a"] for (k, (l, i, j)) in enumerate(ref[:arcs])],
        angmax = [b["angmax"] for (i, b) in ref[:branch]],
        angmin = [b["angmin"] for (i, b) in ref[:branch]],
    )

    data_load_time = time() - time_data_start


    time_model_start = time()

    w = ExaModels.ExaCore()

    va = ExaModels.variable(w, length(data.bus);)

    vm = ExaModels.variable(
        w,
        length(data.bus);
        start = fill!(similar(data.bus, Float64), 1.0),
        lvar = data.vmin,
        uvar = data.vmax,
    )
    pg = ExaModels.variable(w, length(data.gen); lvar = data.pmin, uvar = data.pmax)

    qg = ExaModels.variable(w, length(data.gen); lvar = data.qmin, uvar = data.qmax)

    p = ExaModels.variable(w, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    q = ExaModels.variable(w, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    o = ExaModels.objective(
        w,
        g.cost1 * pg[g.i]^2 + g.cost2 * pg[g.i] + g.cost3 for g in data.gen
    )

    c1 = ExaModels.constraint(w, va[i] for i in data.ref_buses)

    c2 = ExaModels.constraint(
        w,
        p[b.f_idx] - b.c5 * vm[b.f_bus]^2 -
        b.c3 * (vm[b.f_bus] * vm[b.t_bus] * cos(va[b.f_bus] - va[b.t_bus])) -
        b.c4 * (vm[b.f_bus] * vm[b.t_bus] * sin(va[b.f_bus] - va[b.t_bus])) for
        b in data.branch
    )

    c3 = ExaModels.constraint(
        w,
        q[b.f_idx] +
        b.c6 * vm[b.f_bus]^2 +
        b.c4 * (vm[b.f_bus] * vm[b.t_bus] * cos(va[b.f_bus] - va[b.t_bus])) -
        b.c3 * (vm[b.f_bus] * vm[b.t_bus] * sin(va[b.f_bus] - va[b.t_bus])) for
        b in data.branch
    )

    c4 = ExaModels.constraint(
        w,
        p[b.t_idx] - b.c7 * vm[b.t_bus]^2 -
        b.c1 * (vm[b.t_bus] * vm[b.f_bus] * cos(va[b.t_bus] - va[b.f_bus])) -
        b.c2 * (vm[b.t_bus] * vm[b.f_bus] * sin(va[b.t_bus] - va[b.f_bus])) for
        b in data.branch
    )

    c5 = ExaModels.constraint(
        w,
        q[b.t_idx] +
        b.c8 * vm[b.t_bus]^2 +
        b.c2 * (vm[b.t_bus] * vm[b.f_bus] * cos(va[b.t_bus] - va[b.f_bus])) -
        b.c1 * (vm[b.t_bus] * vm[b.f_bus] * sin(va[b.t_bus] - va[b.f_bus])) for
        b in data.branch
    )

    c6 = ExaModels.constraint(
        w,
        va[b.f_bus] - va[b.t_bus] for b in data.branch;
        lcon = data.angmin,
        ucon = data.angmax,
    )
    c7 = ExaModels.constraint(
        w,
        p[b.f_idx]^2 + q[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )
    c8 = ExaModels.constraint(
        w,
        p[b.t_idx]^2 + q[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c9 = ExaModels.constraint(w, b.pd + b.gs * vm[b.i]^2 for b in data.bus)

    c10 = ExaModels.constraint(w, b.qd - b.bs * vm[b.i]^2 for b in data.bus)

    c11 = ExaModels.constraint!(w, c9, a.bus => p[a.i] for a in data.arc)
    c12 = ExaModels.constraint!(w, c10, a.bus => q[a.i] for a in data.arc)

    c13 = ExaModels.constraint!(w, c9, g.bus => -pg[g.i] for g in data.gen)
    c14 = ExaModels.constraint!(w, c10, g.bus => -qg[g.i] for g in data.gen)

    model = ExaModels.ExaModel(w)

    model_build_time = time() - time_model_start


    time_solve_start = time()

    result = NLPModelsIpopt.ipopt(model)

    cost = result.objective

    feasible = result.status == :first_order

    solve_time = time() - time_solve_start
    total_time = time() - time_data_start

    # TODO: once a built-in timer is implemented in ExaModels, report callback timing
    
    # total_callback_time =
    #     nlp_block.evaluator.eval_objective_timer +
    #     nlp_block.evaluator.eval_objective_gradient_timer +
    #     nlp_block.evaluator.eval_constraint_timer +
    #     nlp_block.evaluator.eval_constraint_jacobian_timer +
    #     nlp_block.evaluator.eval_hessian_lagrangian_timer

    println("")
    println("\033[1mSummary\033[0m")
    println("   case........: $(file_name)")
    println("   variables...: $(model.meta.nvar)")
    println("   constraints.: $(model.meta.ncon)")
    println("   feasible....: $(feasible)")
    println("   cost........: $(round(Int, cost))")
    println("   total time..: $(total_time)")
    println("     data time.: $(data_load_time)")
    println("     build time: $(model_build_time)")
    println("     solve time: $(solve_time)")
    # println("      callbacks: $(total_callback_time)")
    println("")
    # println("   callbacks time:")
    # println("   * obj.....: $(nlp_block.evaluator.eval_objective_timer)")
    # println("   * grad....: $(nlp_block.evaluator.eval_objective_gradient_timer)")
    # println("   * cons....: $(nlp_block.evaluator.eval_constraint_timer)")
    # println("   * jac.....: $(nlp_block.evaluator.eval_constraint_jacobian_timer)")
    # println("   * hesslag.: $(nlp_block.evaluator.eval_hessian_lagrangian_timer)")
    # println("")

    va_sol = ExaModels.solution(result, va)
    va_dict = Dict("va_$(i)" => va_sol[b.i] for (i,b) in enumerate(data.bus))
    
    vm_sol = ExaModels.solution(result, vm)
    vm_dict = Dict("vm_$(i)" => vm_sol[b.i] for (i,b) in enumerate(data.bus))

    pg_sol = ExaModels.solution(result, pg)
    pg_dict = Dict("pg_$(i)" => pg_sol[b.i] for (i,b) in enumerate(data.gen))

    qg_sol = ExaModels.solution(result, qg)
    qg_dict = Dict("qg_$(i)" => qg_sol[b.i] for (i,b) in enumerate(data.gen))

    p_sol = ExaModels.solution(result, p)
    p_dict = Dict("p_$(write_out_tuple(ref[:arcs][i]))" => p_sol[i] for (i,b) in enumerate(data.arc))

    q_sol = ExaModels.solution(result, q)
    q_dict = Dict("q_$(write_out_tuple(ref[:arcs][i]))" => q_sol[i] for (i,b) in enumerate(data.arc))

    return Dict(
        "case" => file_name,
        "variables" => model.meta.nvar,
        "constraints" => model.meta.ncon,
        "feasible" => feasible,
        "cost" => cost,
        "time_total" => total_time,
        "time_data" => data_load_time,
        "time_build" => model_build_time,
        "time_solve" => solve_time,
        # "time_callbacks" => total_callback_time,
        "solution" => Dict(
            va_dict...,
            vm_dict...,
            pg_dict...,
            qg_dict...,
            p_dict...,
            q_dict...),
    )
end

write_out_tuple((i,j,k)) = "$(i)_$(j)_$(k)"

if isinteractive() == false
    solve_opf("$(@__DIR__)/data/opf_warmup.m")
end
