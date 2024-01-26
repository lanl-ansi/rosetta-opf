# Rosetta OPF

The AC Optimal Power Flow problem (AC-OPF) is one of the most foundational optimization problems that arises in the design and operations of power networks.
Mathematically the AC-OPF is a large-scale, sparse, non-convex nonlinear continuous optimization problem. 
The models make heavy use of scalar indexing in the definition of the nonlinear functions, and most operations are simple floating point operations (`+`, `*`, `cos`, and `sin`) applied to scalar data.
In practice AC-OPF is most often solved to local optimality conditions using interior point methods.
This project proposes AC-OPF as _proxy-application_ for testing the viability of different nonlinear optimization frameworks, as performant solutions to AC-OPF has proven to be a necessary (but not always sufficient) condition for solving a wide range of industrial network optimization tasks.

### Objectives
- Communicate the technical requirements for solving real-world continuous non-convex mathematical optimization problems.
- Highlight scalability requirements for the problem sizes that occur in practice.
- Provide a consistent implementation for solving AC-OPF in different NLP modeling frameworks.

### Use-case Assumptions
- All decision variables are continuous
- The objective function may be non-convex
- The constraints include a system of equality and inequality functions, which can be non-convex (the equality constraints usually cannot be expressed explicitly as a manifold)
- The constraints can take the form of polynomial and transcendental functions (e.g. `x^2*y^3`, `sin(x)*cos(y)`)
- Derivative computations should be handled by the modeling layer (e.g., via Automatic Differentiation). The user lacks the time or expertise to _hard-code_ Jacobian and Hessian oracles.

### Scope
At present this project is focused on comparing NLP modeling layers that are available in the Julia programming langue.  However, other modeling layers may be considered in the future.

This work is not intended for comparing different nonlinear optimization algorithms, which are often independent of the NLP modeling layer. Consequently, the [Ipopt](https://github.com/jump-dev/Ipopt.jl) solver is used as a standard NLP algorithm whenever it is accessible from the modeling layer.

### Comparison with SciMLBenchmarks

All framework implementations use PowerModels.jl to read input `.m` data files into non-concrete dictionaries, from which the objective and constraint functions are defined.
Reading the data into non-concrete dictionaries simulates a very common non-expert user of NLP tooling, but biases the performance results in favor of symbolic front-ends which perform additional optimizations on the objective and constraint functions before performing the optimization.

Moreover, the `total time` reported by `rosetta-opf` includes Julia's compilation time.
This biases performance against frameworks which generate large amounts of new code that requires compilation.
We do not attempt to remove the compilation overhead, for example, by using a tool like [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl).

For a variation of `rosetta-opf` that removes these factors to focus on core implementation differences of the numerical backends, see the [SciMLBenchmarks adaptation](https://docs.sciml.ai/SciMLBenchmarksOutput/stable/OptimizationFrameworks/optimal_powerflow/).

## Mathematical and Data Models
This work adopts the mathematical model and data format that is used in the IEEE PES benchmark library for AC-OPF, [PGLib-OPF](https://github.com/power-grid-lib/pglib-opf). The Julia package [PowerModels](https://github.com/lanl-ansi/PowerModels.jl) is used for parsing the problem data files and making standard data transformations.

&nbsp;
![The Mathematical Model of the Optimal Power Flow Problem](MODEL.png?raw=true "AC Optimal Power Flow")
&nbsp;

In this formulation the equations encode the following properties: (1) minimization of generator fuel costs; (2) a voltage phase reference angle; (3) power balance (i.e. energy conservation); (4,5) Ohm's Law for the flow power; (6) power flow limits; and (7) angle difference limits.
AC-OPF is naturally modeled as continuous nonlinear optimization problem over complex data and variables. However, the implementations in this repository use the projection into real numbers to support the broadest possible set of optimization modeling frameworks.
In this formulation the complex voltage terms expand into the following expressions, `|Vᵢ|^2 = (vᵢ)^2`, `∠Vᵢ = θᵢ`, `Vᵢ * conj(Vⱼ) = vᵢ*vⱼ*cos(θᵢ-θⱼ) + im*vᵢ*vⱼ*sin(θᵢ-θⱼ)`, which is the source of the transcendental functions in the implementations.

## Code Overview

### AC-OPF Implementations
Each of these files is designed to be _stand-alone_ and can be tested with minimal package dependencies.
Consequently, there is some code replication between implementations.
- `jump.jl`: implementation using [JuMP](https://github.com/jump-dev/JuMP.jl)
- `nlpmodels.jl`: implementation using [ADNLPModels](https://github.com/JuliaSmoothOptimizers/ADNLPModels.jl)
- `nonconvex.jl`: implementation using [Nonconvex](https://github.com/JuliaNonconvex/Nonconvex.jl)
- `optim.jl`: implementation using [Optim](https://github.com/JuliaNLSolvers/Optim.jl)
- `optimization.jl`: implementation using [Optimization](https://github.com/SciML/Optimization.jl)
- `examodels.jl`: implemenation using [ExaModels](https://github.com/exanauts/ExaModels.jl)

### Other Files
- `data/*`: small example datasets
- `test/*`: basic tests for the primary AC-OPF implementations
- `debug/*`: scripts for debugging NLP modeling layers
- `variants/*`: additional variants of the AC-OPF problem

### Running an AC-OPF Case File

By convention each AC-OPF stand-alone file implements a function called `solve_opf` that requires one argument, `file_name`. The argument should be a path to an AC-OPF case file that is compatible with PowerModels' `parse_file` function (usually a [MATPOWER](https://matpower.org/) case file). The `solve_opf` function parses the case file, builds the suitable AC-OPF optimization problem, solves it with Ipopt (or the next-best available algorithm) and returns a dictionary of basic information about the solution process. This dictionary includes items such as,
- `case`: the name of the file that was solved
- `feasible`: if the modeling layer determined the solution at the competition of the solve process satisfies the model constraints 
- `cost`: evaluation of the objective function at the competition of the solve process
- `time_total`: the total wall-clock time of the `solve_opf` function in seconds

The standard steps for testing the `solve_opf` function on a specific AC-OPF data file are,
```
julia --project=.                               # Start Julia with the project environment provided with this repo
julia> include("jump.jl")                       # Load the solve_opf function from one of the example files
julia> solve_opf("data/pglib_opf_case5_pjm.m")  # Run the solve_opf function on a specific AC-OPF case file
```
Note that due to Julia's JIT, it is very likely that the first call to `solve_opf` will take significantly more time than the second call.


## License

This code is provided under a BSD license as part of the Grid Optimization Competition Solvers project, C19076.
