#!/usr/bin/env julia
###### AC-OPF using GalacticOptim ######
#
# Constraint optimization implementation reference:
#    https://github.com/SciML/GalacticOptim.jl/blob/master/test/rosenbrock.jl#L30
# Other AD libraries can be considered:
#    https://galacticoptim.sciml.ai/stable/API/optimization_function/#Defining-Optimization-Functions-Via-AD
# However, ForwardDiff is the only one that is compatible with constraint
# functions.

import ForwardDiff  # GalacticOptim loads with Requires.jl
import GalacticOptim
import Ipopt

include("utils.jl")

function solve_opf(file_name)
    function build_model(data)
        optprob = GalacticOptim.OptimizationFunction(
            data.obj,
            GalacticOptim.AutoForwardDiff();
            cons = data.cons,
        )
        return GalacticOptim.OptimizationProblem(
            optprob,
            data.var_init,
            data.ref,
            lb = data.var_lb,
            ub = data.var_ub,
            lcons = data.con_lb,
            ucons = data.con_ub,
        )
    end
    function solve(model)
        sol = GalacticOptim.solve(model, Ipopt.Optimizer())
        return (cost = sol.minimum, feasible = (sol.retcode == :LOCALLY_SOLVED))
    end
    return solve_opf(build_model, solve, file_name)
end

if isinteractive() == false
    solve_opf("$(@__DIR__)/data/pglib_opf_case5_pjm.m")
end
