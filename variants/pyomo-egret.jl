#!/usr/bin/env julia
###### AC-OPF using Pyomo and Egret via PythonCall ######
#
# implementation reference: https://github.com/lanl-ansi/PowerModelsAnnex.jl/blob/master/src/model/ac-opf.jl
#

import PythonCall

function solve_opf(file_name)
    matpower_parser = PythonCall.pyimport("egret.parsers.matpower_parser")
    egret_acopf = PythonCall.pyimport("egret.models.acopf")
    pyo = PythonCall.pyimport("pyomo.environ")

    time_data_start = time()
    model_data = matpower_parser.create_ModelData(file_name)
    data_load_time = time() - time_data_start

    time_model_start = time()
    model, _ = egret_acopf.create_psv_acopf_model(model_data)
    model_build_time = time() - time_model_start

    time_solve_start = time()
    solver = pyo.SolverFactory("ipopt")
    solver.options["print_timing_statistics"] = "yes"
    results = solver.solve(model, tee=true)
    solve_time = time() - time_solve_start

    total_time = time() - time_data_start

    n_variables = results.problem.number_of_variables
    n_constraints = results.problem.number_of_constraints
    status = results.solver.termination_condition
    feasible_termination_conditions = Set([
        pyo.TerminationCondition.optimal,
        pyo.TerminationCondition.feasible,
        pyo.TerminationCondition.locallyOptimal,
        pyo.TerminationCondition.globallyOptimal,
    ])
    feasible = (status in feasible_termination_conditions)
    objective = nothing
    for obj in model.component_data_objects(pyo.Objective, active=true)
        if objective !== nothing
            # More than one objective?!
            throw(Exception)
        end
        objective = obj
    end
    cost = pyo.value(model.obj)
    cost = PythonCall.pyconvert(Float64, cost)

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
