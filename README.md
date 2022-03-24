# Rosetta OPF

The AC Optimal Power Flow problem (AC-OPF) is one of the most foundational optimization problems that arises in the design and operations of power networks.
Mathematically the AC-OPF is a large-scale, sparse, non-convex nonlinear continuous optimization problem. 
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
- The constraints can take the form of polynomial and transcendental functions (e.g. `x^2*y^3`, `sin(x)*cos(x)`)
- Derivative computations should be handled by the modeling layer (e.g., via Automatic Differentiation). The user lacks the time or expertise to _hard-code_ Jacobian and Hessian oracles.

### Scope
At present this project is focused on comparing NLP modeling layers that are available in the Julia programming langue.  However, other modeling layers may be considered in the future.

This work is not intended for comparing different nonlinear optimization algorithms, which are often independent of the NLP modeling layer. Consequently, the [Ipopt](https://github.com/jump-dev/Ipopt.jl) solver is used as a standard NLP algorithm whenever it is accessible from the modeling layer.

## Mathematical and Data Models
This work adopts the mathematical model and data format that is used in [PGLib-OPF](https://github.com/power-grid-lib/pglib-opf) in the polar voltage form, which is the IEEE PES benchmark library for AC-OPF. The Julia package [PowerModels](https://github.com/lanl-ansi/PowerModels.jl) is used for parsing the problem data files.

## Code Overview

### AC-OPF Implementations
Each of these files is design to be _stand-alone_ and can be tested with minimal package dependencies.
Consequently, there is some code replication between implementations.
- `galacticoptim.jl`: implementation using [GalacticOptim](https://github.com/SciML/GalacticOptim.jl)
- `jump.jl`: implementation using [JuMP](https://github.com/jump-dev/JuMP.jl)
- `nlpmodels.jl`: implementation using [ADNLPModels](https://github.com/JuliaSmoothOptimizers/ADNLPModels.jl)
- `nonconvex.jl`: implementation using [Nonconvex](https://github.com/JuliaNonconvex/Nonconvex.jl)
- `optim.jl`: implementation using [Optim](https://github.com/JuliaNLSolvers/Optim.jl)

### Other Files
- `data/*`: small example datasets
- `test/*`: basic tests for the primary AC-OPF implementations
- `debug/*`: scripts for debugging NLP modeling layers
- `variants/*`: additional variants of the AC-OPF problem

## License

This code is provided under a BSD license as part of the Grid Optimization Competition Solvers project, C19076.
