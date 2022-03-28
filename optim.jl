#!/usr/bin/env julia
###### AC-OPF using Optim ######
#
# implementation reference: https://julianlsolvers.github.io/Optim.jl/stable/#examples/generated/ipnewton_basics/
# currently does not converge to a feasible point, root cause in unclear
# `debug/optim-debug.jl` can be used to confirm it will converge if given a suitable starting point
#

import Optim

include("utils.jl")

function solve_opf(file_name)
    function build_model(data)
        # NOTE: had to change initial guess to be an interior point, otherwise
        #       we get NaN values
        for i in 1:length(data.var_init)
            if data.var_init[i] == data.var_lb[i]
                data.var_init[i] = 0.5 * (data.var_lb[i] + data.var_ub[i])
            end
        end
        df = Optim.TwiceDifferentiable(data.obj, data.var_init)
        dfc = Optim.TwiceDifferentiableConstraints(
            (c, x) -> (c .= data.cons(x)),
            data.var_lb,
            data.var_ub,
            data.con_lb,
            data.con_ub,
        )
        return (data = data, df = df, dfc = dfc)
    end
    function solve(model)
        options = Optim.Options(show_trace=true)
        algorithm = Optim.IPNewton()
        # algorithm = Optim.LBFGS()       # StackOverflowError
        # algorithm = Optim.NelderMead()  # StackOverflowError
        res = Optim.optimize(
            model.df,
            model.dfc,
            model.data.var_init,
            algorithm,
            options,
        )
        display(res)
        # NOTE: confirmed these constraint violations can be eliminated if a
        # better starting point is used.
        sol_eval = model.data.cons(res.minimizer)
        vio_lb = max.(0, model.data.con_lb .- sol_eval)
        vio_ub = max.(0, sol_eval .- model.data.con_ub)
        const_vio = vio_lb .+ vio_ub
        println("total constraint violation: $(sum(const_vio))")
        feasible = sum(const_vio) <= 1e-6
        if !feasible
            @warn "Optim optimize failed to satify the problem constraints"
        end
        return (cost = res.minimum, feasible = feasible)
    end
    return solve_opf(build_model, solve, file_name)
end

if isinteractive() == false
    solve_opf("$(@__DIR__)/data/pglib_opf_case5_pjm.m")
end
