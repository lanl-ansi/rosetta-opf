using Test

test_case = "../data/pglib_opf_case5_pjm.m"
test_case_cost = 17551.891

@testset "Rosetta OPF" begin

    @testset "GalacticOptim" begin
        include("../galacticoptim.jl")
        result = solve_opf(test_case)

        @test result["feasible"]
        @test isapprox(result["cost"], test_case_cost)
    end

    @testset "JuMP" begin
        include("../jump.jl")
        result = solve_opf(test_case)

        @test result["feasible"]
        @test isapprox(result["cost"], test_case_cost)
    end

    @testset "NLPModels" begin
        include("../nlpmodels.jl")
        result = solve_opf(test_case)

        @test result["feasible"]
        @test isapprox(result["cost"], test_case_cost)
    end

    #=
    # currently blocked by https://github.com/JuliaNonconvex/Nonconvex.jl/issues/130
    @testset "NonConvex" begin
        include("../nonconvex.jl")
        result = solve_opf(test_case)

        @test result["feasible"]
        @test isapprox(result["cost"], test_case_cost)
    end
    =#

    # does not converge to a feasible solution
    @testset "Optim" begin
        include("../optim.jl")
        result = solve_opf(test_case)

        @test !result["feasible"]
        #@test isapprox(result["cost"], test_case_cost)
    end

end