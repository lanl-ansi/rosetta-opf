#=
    Implement the vectorized OPF model proposed in:

    Feasible Path Identification in Optimal Power Flow with Sequential Convex Restriction
    Dongchan Lee, Konstantin Turitsyn, Daniel K. Molzahn, Line A. Roald

    Preprint: https://arxiv.org/abs/1906.09483
    Code:     https://github.com/dclee131/PowerFlowCVXRS

    The vectorized version of the OPF leads to a compact
    expression tree, with all the nonlinearities factorized inside
    a single nonlinear basis ψ. Comparing to the scalar version,
    the vectorized formulation allows to streamline the evaluation
    of the derivatives.

    The expression tree is given directly as:

                                    /------ (power flow equations)
                                   /
    (vm, va, pg) --->  ψ(vm, va) -/-------- (reactive power generations qg)
                                  \
                                   \------- (line flow)

    In addition, this model considers only a subset of the power
    flow equations, and lead to a more compact KKT system than
    the formulation adopted in PowerModels.jl.

=#

time_start = time()

using LinearAlgebra
using SparseArrays
using JuMP, PowerModels

pkg_load_time = time() - time_start

time_start = time()

file_name = "../data/pglib_opf_case5_pjm.m"
# file_name = "../data/pglib_opf_case118_ieee.m"

#=
    Import data with PowerModels
=#
data = PowerModels.parse_file(file_name)
PowerModels.standardize_cost_terms!(data, order=2)
PowerModels.calc_thermal_limits!(data)
ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

data_load_time = time() - time_start

gen = data["gen"]
bus = data["bus"]
branch = data["branch"]
shunt = data["shunt"]
load = data["load"]
nbus = length(data["bus"])
ngen = length(data["gen"])
nlines = length(data["branch"])
nloads = length(data["load"])
nshunts = length(data["shunt"])
## GENERATORS
pg0 = Float64[gen["$i"]["pg"] for i in 1:ngen]
pgmin = Float64[gen["$i"]["pmin"] for i in 1:ngen]
pgmax = Float64[gen["$i"]["pmax"] for i in 1:ngen]
qgmin = Float64[gen["$i"]["qmin"] for i in 1:ngen]
qgmax = Float64[gen["$i"]["qmax"] for i in 1:ngen]
gen2bus = Int[gen["$i"]["gen_bus"] for i in 1:ngen]
## BUSES
busid = Int[bus["$i"]["bus_i"] for i in 1:nbus]
bustype = Int[bus["$i"]["bus_type"] for i in 1:nbus]
vm0 = Float64[bus["$i"]["vm"] for i in 1:nbus]
va0 = Float64[bus["$i"]["va"] for i in 1:nbus]
vmin = Float64[bus["$i"]["vmin"] for i in 1:nbus]
vmax = Float64[bus["$i"]["vmax"] for i in 1:nbus]
## SHUNTS
b_sh = Int[shunt["$i"]["shunt_bus"] for i in 1:nshunts]
gs_sh = Float64[shunt["$i"]["gs"] for i in 1:nshunts]
bs_sh = Float64[shunt["$i"]["bs"] for i in 1:nshunts]
shunts = gs_sh .+ 1im .* bs_sh
## LOADS
load2bus = Int[load["$i"]["load_bus"] for i in 1:nloads]
pd = Float64[load["$i"]["pd"] for i in 1:nloads]
qd = Float64[load["$i"]["qd"] for i in 1:nloads]
## BRANCHES
f_bus = Int[branch["$i"]["f_bus"] for i in 1:nlines]
t_bus = Int[branch["$i"]["t_bus"] for i in 1:nlines]
tap = Float64[branch["$i"]["tap"] for i in 1:nlines]
stat = Float64[branch["$i"]["br_status"] for i in 1:nlines]
br_r = Float64[branch["$i"]["br_r"] for i in 1:nlines]
br_x = Float64[branch["$i"]["br_x"] for i in 1:nlines]
b_fr = Float64[branch["$i"]["b_fr"] for i in 1:nlines]
b_to = Float64[branch["$i"]["b_to"] for i in 1:nlines]
br_b = b_fr + b_to
shift = Float64[branch["$i"]["shift"] for i in 1:nlines]
sline_max = Float64[branch["$i"]["rate_a"] for i in 1:nlines]
## COSTS
c2 = Float64[gen["$i"]["cost"][1] for i in 1:ngen]
c1 = Float64[gen["$i"]["cost"][2] for i in 1:ngen]
c0 = Float64[gen["$i"]["cost"][3] for i in 1:ngen]

#=
    Buses' classification
=#
pv = findall(isequal(2), bustype)
pq = findall(isequal(1), bustype)
ref = findall(isequal(3), bustype)
refgen = findfirst(isequal(ref[1]), gen2bus)
npv = length(pv)
npq = length(pq)
nref = length(ref)

#=
    Build incidence matrices encoding the problem's topology
=#
Cg = sparse(gen2bus, 1:ngen, ones(ngen), nbus, ngen)
Cl = sparse(load2bus, 1:nloads, ones(nloads), nbus, nloads)
Cf = sparse(f_bus, 1:nlines, ones(nlines), nbus, nlines)
Ct = sparse(t_bus, 1:nlines, ones(nlines), nbus, nlines)
E = Cf - Ct
Ysh = sparse(b_sh, b_sh, shunts, nbus, nbus)

#=
    Build admittance matrices
    (take expressions from MATPOWER)
=#
Ys = stat ./ (br_r + 1im * br_x)       # series admittance
Bc = stat .* br_b                      # line charging susceptance
tap = tap .* exp.(1im*pi/180 .* shift) # add phase shifters
Ytt = Ys .+ 1im .* Bc ./ 2.0
Yff = Ytt ./ (tap .* conj.(tap))
Yft = - Ys ./ conj.(tap)
Ytf = - Ys ./ tap

i = [1:nlines; 1:nlines]
Yf = sparse(i, [f_bus; t_bus], [Yff; Yft], nlines, nbus)
Yt = sparse(i, [f_bus; t_bus], [Ytf; Ytt], nlines, nbus)
Ybus = Cf * Yf + Ct * Yt + Ysh

Yc = Cf * Diagonal(Yft) + Ct * Diagonal(Ytf)
Ys = Cf * Diagonal(Yft) - Ct * Diagonal(Ytf)
Yd = Cf * Diagonal(Yff) * Cf' + Ct * Diagonal(Ytt) * Ct' + Ysh

yff = Diagonal(Yff)
yft = Diagonal(Yft)
ytf = Diagonal(Ytf)
ytt = Diagonal(Ytt)

#=
    Compact matrices
=#
# Eq (5): power flow equations
M = [ real(Yc)  imag(Ys)  real(Yd);
     -imag(Yc)  real(Ys) -imag(Yd)]
# Eqs (9-10): line flow constraints
Lfp = [real(yft)  imag(yft)  real(yff) * Cf']
Lfq = [-imag(yft) real(yft) -imag(yff) * Cf']
Ltp = [real(ytf)  -imag(ytf)  real(ytt) * Ct']
Ltq = [-imag(ytf) -real(ytf) -imag(ytt) * Ct']

# Eq (7): Power flow (equality)
M_eq = M[[pv; pq; nbus .+ pq], :]
Cg_eq = [Cg[pv, :] ; spzeros(2 * npq, ngen)]
τ_eq = [Cl[[pv; pq], :] * pd; Cl[pq, :] * qd]
# Eq (8): Power flow (inequality)
M_ineq = M[[ref; nbus .+ ref; nbus .+ pv], :]
τ_ineq = [Cl[ref, :] * pd; Cl[[ref; pv], :] * qd]
pmin = [Cg[ref, :] * pgmin; Cg[[ref;pv], :] *qgmin]
pmax = [Cg[ref, :] * pgmax; Cg[[ref;pv], :] *qgmax]

#=
    Build JuMP model
=#
time_start = time()

acopf = Model(Ipopt.Optimizer)

## 1. Variables

@variable(acopf, vmin[i] <= vm[i=1:nbus] <= vmax[i], start=vm0[i])
@variable(acopf, va[i=1:nbus], start=va0[i])
@variable(acopf, pgmin[i] <= pg[i=1:ngen] <= pgmax[i], start=pg0[i])
# nonlinear basis (encodes all the nonlinearities in the problem)
@variable(acopf, ψsin[i=1:nlines])
@variable(acopf, ψcos[i=1:nlines])
@variable(acopf, ψq[i=1:nbus])
# line flow
@variable(acopf, sfp[i = 1:nlines])
@variable(acopf, stp[i = 1:nlines])
@variable(acopf, sfq[i = 1:nlines])
@variable(acopf, stq[i = 1:nlines])

## 2. Constraints

@constraint(acopf, va[ref] .== 0)
# Angle difference
@expression(acopf, φ, E' * va)
# The nonlinear basis encodes all the nonlinearities in the problem.
@NLconstraint(acopf, [i=1:nlines], ψsin[i] == vm[f_bus[i]] * vm[t_bus[i]] * sin(φ[i]))
@NLconstraint(acopf, [i=1:nlines], ψcos[i] == vm[f_bus[i]] * vm[t_bus[i]] * cos(φ[i]))
@NLconstraint(acopf, [i=1:nbus],  ψq[i] == vm[i]^2)

# eq(11b): recover power flow equations with sparse operations
@constraint(acopf, τ_eq - Cg_eq * pg + M_eq * [ψcos; ψsin; ψq] .== 0)

# eq(11c) - eq(20): active power bounds on slack + reactive power bounds
@constraint(acopf, pmin .<= M_ineq * [ψcos; ψsin; ψq] .+ τ_ineq .<= pmax)

# Line-flow constraints
# eq(11e) - eq(20)
@constraint(acopf, sfp .== Lfp * [ψcos; ψsin; ψq])
@constraint(acopf, sfq .== Lfq * [ψcos; ψsin; ψq])
@constraint(acopf, stp .== Ltp * [ψcos; ψsin; ψq])
@constraint(acopf, stq .== Ltq * [ψcos; ψsin; ψq])
@constraint(acopf, sfp.^2 .+ sfq.^2  .<= sline_max.^2)
@constraint(acopf, stp.^2 .+ stq.^2  .<= sline_max.^2)

# Recover active power generation at slack node
@constraint(acopf, pg[refgen] .== M_ineq[1:1, :] * [ψcos; ψsin; ψq] .+ τ_ineq[1:1])

## 3. Objective

@objective(acopf, Min, dot(c2, pg.^2) + dot(c1, pg) + sum(c0))

model_build_time = time() - time_start

#=
    Resolution
=#

time_start = time()

optimize!(acopf)

solve_time = time() - time_start

#=
    Analysis
=#

nlp_block = MOI.get(acopf, MOI.NLPBlock())

println("")
println("\033[1mSummary\033[0m")
println("   case......: $(file_name)")
println("   cost......: $(round(Int, JuMP.objective_value(acopf)))")
println("   pkg time..: $(pkg_load_time)")
println("   data time.: $(data_load_time)")
println("   build time: $(model_build_time)")
println("   solve time: $(solve_time)")
println("   callbacks time:")
println("   * obj.....: $(nlp_block.evaluator.eval_objective_timer)")
println("   * grad....: $(nlp_block.evaluator.eval_objective_gradient_timer)")
println("   * cons....: $(nlp_block.evaluator.eval_constraint_timer)")
println("   * jac.....: $(nlp_block.evaluator.eval_constraint_jacobian_timer)")
println("   * hesslag.: $(nlp_block.evaluator.eval_hessian_lagrangian_timer)")
println("")
