#!/usr/bin/env julia
###### AC-OPF using Pyomo via PythonCall ######
#
# implementation reference: https://github.com/lanl-ansi/PowerModelsAnnex.jl/blob/master/src/model/ac-opf.jl
#

import PowerModels
import PythonCall

function solve_opf(file_name)
    pyo = PythonCall.pyimport("pyomo.environ")

    time_data_start = time()

    # User PowerModels to parse matpower file
    data = PowerModels.parse_file(file_name)
    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)

    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    data_load_time = time() - time_data_start

    time_model_start = time()

    m = pyo.ConcreteModel()

    #
    # Declare sets
    #
    buses = [key for key in keys(ref[:bus])]
    # Assuming generators that are not in-service have already been filtered
    generators = [key for key in keys(ref[:gen])]
    branches = [key for key in keys(ref[:branch])]
    arcs_fr = [key for key in ref[:arcs_from]]
    arcs_to = [key for key in ref[:arcs_to]]
    arcs = [key for key in ref[:arcs]]

    m.buses = pyo.Set(initialize=buses)
    m.generators = pyo.Set(initialize=generators)
    m.branches = pyo.Set(initialize=branches)
    m.arcs_fr = pyo.Set(initialize=arcs_fr)
    m.arcs_to = pyo.Set(initialize=arcs_to)
    m.arcs = pyo.Set(initialize=arcs)
    ###

    #
    # Declare variables
    #
    vm_bounds = Dict(i => (ref[:bus][i]["vmin"], ref[:bus][i]["vmax"]) for i in buses)
    pg_bounds = Dict(i => (ref[:gen][i]["pmin"], ref[:gen][i]["pmax"]) for i in generators)
    qg_bounds = Dict(i => (ref[:gen][i]["qmin"], ref[:gen][i]["qmax"]) for i in generators)
    branch_pq_bounds = Dict(
        (l, i, j) => (-ref[:branch][l]["rate_a"], ref[:branch][l]["rate_a"])
        for (l, i, j) in arcs
    )

    m.va = pyo.Var(m.buses, initialize=0.0)
    m.vm = pyo.Var(m.buses, bounds=vm_bounds, initialize=1.0)
    m.pg = pyo.Var(m.generators, bounds=pg_bounds)
    m.qg = pyo.Var(m.generators, bounds=qg_bounds)
    m.p = pyo.Var(m.arcs, bounds=branch_pq_bounds)
    m.q = pyo.Var(m.arcs, bounds=branch_pq_bounds)
    ###

    #
    # Declare objective
    #
    m.obj = pyo.Objective(
        sense=pyo.minimize,
        expr=sum(
            gen["cost"][1]*m.pg[i]^2 + gen["cost"][2]*m.pg[i] + gen["cost"][3]
            for (i, gen) in ref[:gen]
        )
    )
    ###

    #
    # Declare constraints
    #
    for (i, bus) in ref[:ref_buses]
        m.va[i].fix(0.0)
    end

    # Power balance constraints
    bus_loads = Dict(i => [ref[:load][l] for l in ref[:bus_loads][i]] for i in buses)
    bus_shunts = Dict(i => [ref[:shunt][s] for s in ref[:bus_shunts][i]] for i in buses)

    _p_balance_rule = Dict(
        i => (
            sum(m.p[a] for a in ref[:bus_arcs][i]; init = 0)
            == sum(m.pg[g] for g in ref[:bus_gens][i]; init = 0)
            - sum(load["pd"] for load in bus_loads[i]; init = 0)
            - m.vm[i]^2 * sum(shunt["gs"] for shunt in bus_shunts[i]; init = 0)
        )
        for i in buses
    )
    m.p_balance = pyo.Constraint(m.buses, rule=_p_balance_rule)

    _q_balance_rule = Dict(
        i => (
            sum(m.q[a] for a in ref[:bus_arcs][i]; init = 0)
            == sum(m.qg[g] for g in ref[:bus_gens][i]; init = 0)
            - sum(load["qd"] for load in bus_loads[i]; init = 0)
            + m.vm[i]^2 * sum(shunt["bs"] for shunt in bus_shunts[i]; init = 0)
        )
        for i in buses
    )
    m.q_balance = pyo.Constraint(m.buses, rule=_q_balance_rule)

    # Branch constraints
    branch_fr_idx = Dict(l => (l, i, j) for (l, i, j) in arcs_fr)
    branch_to_idx = Dict(l => (l, i, j) for (l, i, j) in arcs_to)
    branch_vm_fr = Dict(i => m.vm[branch["f_bus"]] for (i, branch) in ref[:branch])
    branch_vm_to = Dict(i => m.vm[branch["t_bus"]] for (i, branch) in ref[:branch])
    branch_va_fr = Dict(i => m.va[branch["f_bus"]] for (i, branch) in ref[:branch])
    branch_va_to = Dict(i => m.va[branch["t_bus"]] for (i, branch) in ref[:branch])

    branch_gb = Dict(i => PowerModels.calc_branch_y(branch) for (i, branch) in ref[:branch])
    branch_g = Dict(i => g for (i, (g, _)) in branch_gb)
    branch_b = Dict(i => b for (i, (_, b)) in branch_gb)
    branch_g_fr = Dict(i => branch["g_fr"] for (i, branch) in ref[:branch])
    branch_g_to = Dict(i => branch["g_to"] for (i, branch) in ref[:branch])
    branch_b_fr = Dict(i => branch["b_fr"] for (i, branch) in ref[:branch])
    branch_b_to = Dict(i => branch["b_to"] for (i, branch) in ref[:branch])

    branch_trti = Dict(i => PowerModels.calc_branch_t(branch) for (i, branch) in ref[:branch])
    branch_tr = Dict(i => tr for (i, (tr, _)) in branch_trti)
    branch_ti = Dict(i => ti for (i, (_, ti)) in branch_trti)
    branch_ttm = Dict(i => branch_tr[i]^2 + branch_ti[i]^2 for i in branches)

    m.p_flow_eq_from = pyo.Constraint(m.branches)
    m.q_flow_eq_from = pyo.Constraint(m.branches)
    m.p_flow_eq_to = pyo.Constraint(m.branches)
    m.q_flow_eq_to = pyo.Constraint(m.branches)
    m.va_difference_limit = pyo.Constraint(m.branches)
    m.apparent_power_limit_from = pyo.Constraint(m.branches)
    m.apparent_power_limit_to = pyo.Constraint(m.branches)

    for (i, branch) in ref[:branch]

        f_idx = branch_fr_idx[i]
        t_idx = branch_to_idx[i]

        p_fr = m.p[f_idx]
        q_fr = m.q[f_idx]
        p_to = m.p[t_idx]
        q_to = m.q[t_idx]

        vm_fr = branch_vm_fr[i]
        vm_to = branch_vm_to[i]
        va_fr = branch_va_fr[i]
        va_to = branch_va_to[i]

        g = branch_g[i]
        b = branch_b[i]
        tr = branch_tr[i]
        ti = branch_ti[i]
        ttm = branch_ttm[i]
        g_fr = branch_g_fr[i]
        b_fr = branch_b_fr[i]
        g_to = branch_g_to[i]
        b_to = branch_b_to[i]

        m.p_flow_eq_from[i] = (
            p_fr
            == (g + g_fr)/ttm * vm_fr^2
            + (-g*tr + b*ti)/ttm * (vm_fr*vm_to*pyo.cos(va_fr - va_to))
            + (-b*tr - g*ti)/ttm * (vm_fr*vm_to*pyo.sin(va_fr - va_to))
        )
        m.q_flow_eq_from[i] = (
            q_fr
            == - (b + b_fr)/ttm * vm_fr^2
            - (-b*tr - g*ti)/ttm * (vm_fr*vm_to*pyo.cos(va_fr - va_to))
            + (-g*tr + b*ti)/ttm * (vm_fr*vm_to*pyo.sin(va_fr - va_to))
        )
        m.p_flow_eq_to[i] = (
            p_to
            == (g + g_to)*vm_to^2
            + (-g*tr - b*ti)/ttm * (vm_to*vm_fr*pyo.cos(va_to - va_fr))
            + (-b*tr + g*ti)/ttm * (vm_to*vm_fr*pyo.sin(va_to - va_fr))
        )
        m.q_flow_eq_to[i] = (
            q_to
            == - (b + b_to)*vm_to^2
            - (-b*tr + g*ti)/ttm * (vm_to*vm_fr*pyo.cos(va_to - va_fr))
            + (-g*tr - b*ti)/ttm * (vm_to*vm_fr*pyo.sin(va_to - va_fr))
        )
        m.va_difference_limit[i] = (branch["angmin"], va_fr - va_to, branch["angmax"])
        m.apparent_power_limit_from[i] = (p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
        
    end

    model_build_time = time() - time_model_start

    #
    # Solve model
    #
    solve_time_start = time()

    solver = pyo.SolverFactory("ipopt")
    solver.options["linear_solver"] = "ma27"
    #solver.options["print_timing_statistics"] = "yes"
    results = solver.solve(m, tee=true)

    solve_time = time() - solve_time_start
    ###

    total_time = time() - time_data_start

    n_variables = PythonCall.pyconvert(Int, results.problem.number_of_variables)
    n_constraints = PythonCall.pyconvert(Int, results.problem.number_of_constraints)
    status = results.solver.termination_condition
    feasible_termination_conditions = Set([
        pyo.TerminationCondition.optimal,
        pyo.TerminationCondition.feasible,
        pyo.TerminationCondition.locallyOptimal,
        pyo.TerminationCondition.globallyOptimal,
    ])
    feasible = (status in feasible_termination_conditions)

    cost = pyo.value(m.obj)
    cost = PythonCall.pyconvert(Float64, cost)

    println("")
    println("\033[1mSummary\033[0m")
    println("   case........: $(file_name)")
    println("   variables...: $(n_variables)")
    println("   constraints.: $(n_constraints)")
    println("   feasible....: $(feasible)")
    println("   cost........: $(round(Int, cost))")
    println("   total time..: $(total_time)")
    println("     data time.: $(data_load_time)")
    println("     build time: $(model_build_time)")
    println("     solve time: $(solve_time)")
    # TODO: Get callback times by writing ipopt log to a file and reading
    # the file to get the detailed timing statistics.
    #println("      callbacks: $(total_callback_time)")
    #println("")
    #println("   callbacks time:")
    #println("   * obj.....: $()")
    #println("   * grad....: $()")
    #println("   * cons....: $()")
    #println("   * jac.....: $()")
    #println("   * hesslag.: $()")
    println("")

    return Dict(
        "case" => file_name,
        "variables" => n_variables,
        "constraints" => n_constraints,
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
