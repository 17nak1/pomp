# Subplex
## Minimization Of A Function By The Subplex Algorithm
Subplx uses the subplex method to solve unconstrained
optimization problems.  The method is well suited for
optimizing objective functions that are noisy or are
discontinuous at the solution.
### Usage
**subplex.f** = The function to be minimized. Its first argument must be the vector of parameters to be optimized over. It should return a scalar result.

**subplex.x0** = Object of an initial guess of the parameters to be optimized over.

**subplex.tol** = Tolerance. The good option can be 2.220446e-16.

**subplex.run()**;

### Value
*subplex.run()* returns an array containing the following:

**x:** computed optimum.

**fx:** value of f at x.

**nfe:** number of function evaluations

**msg:** message if there is a warning.