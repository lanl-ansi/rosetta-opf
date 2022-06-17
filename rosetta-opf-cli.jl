#!/usr/bin/env julia

using ArgParse

s = ArgParseSettings()

@add_arg_table s begin
    "--solver", "-s"
        help = "solver to load"
        required = true
    "--file", "-f"
        help = "the opf data file (.m)"
        required = true
end

args = parse_args(s)


if args["solver"] == "jump"
    include("jump.jl")
elseif args["solver"] == "nlpmodels"
    include("nlpmodels.jl")
elseif args["solver"] == "nonconvex"
    include("nonconvex.jl")
elseif args["solver"] == "optim"
    include("optim.jl")
elseif args["solver"] == "optimization"
    include("optimization.jl")
elseif args["solver"] == "jump-nl"
    include("variants/jump-nl.jl")
elseif args["solver"] == "jump-sad"
    include("variants/jump-symbolic-ad.jl")
else
    error("unknwon solver type $(args["solver"])")
end


result = solve_opf(args["file"])


data = [
    "DATA",
    result["case"],
    result["variables"],
    result["constraints"],
    result["feasible"],
    result["cost"],
    result["time_total"],
    result["time_data"],
    result["time_build"],
    result["time_solve"],
    get(result, "time_callbacks", 0)
]
println(join(data, ", "))

