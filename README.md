# Rosetta OPF

Various AC OPF Implementations.

Objectives:
- communicate the technical requirements for solving real-world contious non-convex mathematical optimization problems.
- test NLP modeling frameworks performance and scalablity
Assumes that derivative computations should be handeled by the modeling layer, the user lacks the time or expertise to _hard-code_ Jacobian and Hessian oracles.

Models that pass derivative computations on to the underlying solvers (e.g. QCQP, SDP) are not generally of interest to this effort (e.g. Convex.jl) as this only tests the modeling layer's ability to encode the problem for the underlying solver and not an iterative loop between the modeling layer (for computing Jacobian and Hessian matricies) and the solver, as is the case for most NLP solvers supporting non-convex functions.


## License

This code is provided under a BSD license as part of the Grid Optimization Competition Solvers project, C19076.
