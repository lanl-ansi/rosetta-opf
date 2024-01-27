using Test

include(joinpath(dirname(@__DIR__), "test", "validator.jl"))

@testset "Rosetta OPF" begin
    framework = ENV["MODELING_FRAMEWORK"]
    include(joinpath(@__DIR__, "$framework.jl"))
    @testset "$case" for case in [
        "case5_pjm",
        "case14_ieee",
        "case24_ieee_rts",
    ]
        test_case = joinpath(dirname(@__DIR__), "data/pglib_opf_$case.m")
        result = solve_opf(test_case)
        validate_result(test_case, result)
    end
end
