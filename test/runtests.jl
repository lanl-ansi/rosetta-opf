using Test

test_case = "../data/pglib_opf_case5_pjm.m"
test_case_cost = 17551.891

result_keys_required = ["case", "variables", "constraints",
    "feasible", "cost", "time_total", "time_data", "time_build", "time_solve"]

function check_result_keys(result::Dict)
    for k in result_keys_required
        if !haskey(result, k)
            @warn "result dict missing key \"$(k)\""
            return false
        end
    end
    return true
end


@testset "Rosetta OPF" begin

    @testset "Optimization" begin
        include("../optimization.jl")
        result = solve_opf(test_case)

        @test check_result_keys(result)
        @test result["feasible"]
        @test isapprox(result["cost"], test_case_cost)
    end

    @testset "JuMP" begin
        include("../jump.jl")
        result = solve_opf(test_case)

        @test check_result_keys(result)
        @test result["feasible"]
        @test isapprox(result["cost"], test_case_cost)
    end

    @testset "NLPModels" begin
        include("../nlpmodels.jl")
        result = solve_opf(test_case)

        @test check_result_keys(result)
        @test result["feasible"]
        @test isapprox(result["cost"], test_case_cost)
    end

    @testset "NonConvex" begin
        include("../nonconvex.jl")
        result = solve_opf(test_case)

        @test check_result_keys(result)
        @test result["feasible"]
        @test isapprox(result["cost"], test_case_cost)
    end

    # does not converge to a feasible solution
    @testset "Optim" begin
        include("../optim.jl")
        result = solve_opf(test_case)

        @test check_result_keys(result)
        @test !result["feasible"]
        #@test isapprox(result["cost"], test_case_cost)
    end

end
