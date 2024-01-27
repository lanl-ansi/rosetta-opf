using Test

include("validator.jl")

@testset "Rosetta OPF" begin
    @testset "$framework" for framework in [
        "examodels",
        "jump",
        "nlpmodels",
        "nonconvex",
        # "optim", # does not converge to feasible solution
        "optimization",
    ]
        include(joinpath(dirname(@__DIR__), "$framework.jl"))
        @testset "$case" for case in [
            "opf_warmup.m",
            "pglib_opf_case5_pjm.m",
            "pglib_opf_case14_ieee.m",
            "pglib_opf_case24_ieee_rts.m",
        ]
            test_case = joinpath(dirname(@__DIR__), "data", case)
            result = solve_opf(test_case)
            validate_result(test_case, result)
        end
    end
end
