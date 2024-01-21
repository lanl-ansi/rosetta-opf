using Test

function validate_result(file_name::String, result)
    for k in [
        "case",
        "variables",
        "constraints",
        "feasible",
        "cost",
        "time_total",
        "time_data",
        "time_build",
        "time_solve",
        "solution",
    ]
        @test haskey(result, k)
    end
    @test result["feasible"] === true
    validate_solution(file_name, result["solution"], result["cost"])
    return
end

function validate_solution(
    file_name::String,
    x::Dict{String,Float64},
    cost::Float64;
    atol = 1e-6,
    rtol = 1e-6,
)
    data = PowerModels.parse_file(file_name)
    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    va = Dict{Int,Float64}(i => x["va_$i"] for i in keys(ref[:bus]))
    vm = Dict{Int,Float64}()
    for i in keys(ref[:bus])
        vm[i] = x["vm_$i"]
        @test vm[i] <= ref[:bus][i]["vmax"] + atol
        @test ref[:bus][i]["vmin"] <= vm[i] + atol
    end
    pg = Dict{Int,Float64}()
    for i in keys(ref[:gen])
        pg[i] = x["pg_$i"]
        @test pg[i] <= ref[:gen][i]["pmax"] + atol
        @test ref[:gen][i]["pmin"] <= pg[i] + atol
    end
    qg = Dict{Int,Float64}()
    for i in keys(ref[:gen])
        qg[i] = x["qg_$i"]
        @test qg[i] <= ref[:gen][i]["qmax"] + atol
        @test ref[:gen][i]["qmin"] <= qg[i] + atol
    end
    p = Dict{NTuple{3,Int},Float64}()
    for (l,i,j) in ref[:arcs]
        p[(l, i, j)] = x["p_$(l)_$(i)_$(j)"]
        @test p[(l, i, j)] <= ref[:branch][l]["rate_a"] + atol
        @test -ref[:branch][l]["rate_a"] <= p[(l, i, j)] + atol
    end
    q = Dict{NTuple{3,Int},Float64}()
    for (l,i,j) in ref[:arcs]
        q[(l, i, j)] = x["q_$(l)_$(i)_$(j)"]
        @test q[(l, i, j)] <= ref[:branch][l]["rate_a"] + atol
        @test -ref[:branch][l]["rate_a"] <= q[(l, i, j)] + atol
    end
    # Test objective function
    @test isapprox(
        sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]),
        cost;
        atol,
        rtol,
    )
    for (i,bus) in ref[:ref_buses]
        @test isapprox(va[i], 0; atol = atol)
    end

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        @test isapprox(
            sum(p[a] for a in ref[:bus_arcs][i]; init = 0.0),
            sum(pg[g] for g in ref[:bus_gens][i]; init = 0.0) -
                sum(load["pd"] for load in bus_loads; init = 0.0) -
                sum(shunt["gs"] for shunt in bus_shunts; init = 0.0)*vm[i]^2;
            atol,
            rtol,
        )

        @test isapprox(
            sum(q[a] for a in ref[:bus_arcs][i]; init = 0.0),
            sum(qg[g] for g in ref[:bus_gens][i]; init = 0.0) -
                sum(load["qd"] for load in bus_loads; init = 0.0) +
                sum(shunt["bs"] for shunt in bus_shunts; init = 0.0)*vm[i]^2;
            atol,
            rtol,
        )
    end

    # Branch power flow physics and limit constraints
    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        vm_fr = vm[branch["f_bus"]]
        vm_to = vm[branch["t_bus"]]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        ttm = tr^2 + ti^2
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]

        # From side of the branch flow
        @test isapprox(
            p_fr,
            (g+g_fr)/ttm*vm_fr^2 + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to));
            atol,
            rtol,
        )
        @test isapprox(
            q_fr,
            -(b+b_fr)/ttm*vm_fr^2 - (-b*tr-g*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to));
            atol,
            rtol,
        )

        # To side of the branch flow
        @test isapprox(
            p_to,
            (g+g_to)*vm_to^2 + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr));
            atol,
            rtol,
        )
        @test isapprox(
            q_to,
            -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr));
            atol,
            rtol,
        )

        # Voltage angle difference limit
        @test branch["angmin"] <= va_fr - va_to + atol
        @test va_fr - va_to <= branch["angmax"] + atol

        # Apparent power limit, from side and to side
        @test p_fr^2 + q_fr^2 <= branch["rate_a"]^2 + atol
        @test p_to^2 + q_to^2 <= branch["rate_a"]^2 + atol
    end
    return
end

@testset "Rosetta OPF" begin
    @testset "$framework" for framework in [
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
