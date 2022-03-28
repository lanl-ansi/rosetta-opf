#!/usr/bin/env julia
###### AC-OPF using ADNLPModels ######
#
# Implementation reference:
#   https://juliasmoothoptimizers.github.io/ADNLPModels.jl/stable/tutorial/
# Other AD libraries can be considered:
#   https://juliasmoothoptimizers.github.io/ADNLPModels.jl/stable/

import ADNLPModels
import NLPModelsIpopt

include("utils.jl")

function solve_opf(file_name)
    function build_model(data)
        return ADNLPModels.ADNLPModel(
            data.obj,
            data.var_init,
            data.var_lb,
            data.var_ub,
            data.cons,
            data.con_lb,
            data.con_ub,
        )
    end
    function solve(model)
        output = NLPModelsIpopt.ipopt(model)
        return (cost = output.objective, feasible = output.primal_feas <= 1e-6)
    end
    return solve_opf(build_model, solve, file_name)
end

if isinteractive() == false
    solve_opf("$(@__DIR__)/data/pglib_opf_case5_pjm.m")
end
